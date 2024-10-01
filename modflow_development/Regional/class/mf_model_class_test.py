import os
from os.path import join, exists, dirname, expanduser
import sys
import glob
from importlib import reload

import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import geopandas as gpd

class CosumnesWS:
    def __init__(self):
        self.usr_dir = expanduser('~')
        self.doc_dir = join(self.usr_dir, 'Documents')
        self.gwfm_dir = join(self.usr_dir, 'Box','research_cosumnes','GWFlowModel')
        self.sfr_dir = join(self.gwfm_dir, 'SFR_data/')
        self.model_ws = None

    def add_path(self, fxn_dir):
        """ Insert fxn directory into first position on path so local functions supercede the global"""
        if fxn_dir not in sys.path:
            sys.path.insert(0, fxn_dir)

    def setup_environment(self):
        self.add_path(join(self.doc_dir, 'GitHub/flopy'))
        import flopy

        self.add_path(join(self.doc_dir, 'GitHub/CosumnesRiverRecharge/python_utilities'))
        self.add_path(join(self.doc_dir, 'GitHub', 'CosumnesRiverRecharge', 'modflow_development'))
        
        import DIS.f_dis_utilities
        reload(DIS.f_dis_utilities)
        from DIS.f_dis_utilities import set_layer_botm
        
        from mf_utility import get_layer_from_elev, param_load, find_isolated_active
        
class CosumnesModel(CosumnesWS):
    def __init__(self, parent):
        self.parent = parent
        
    def load_model_domain(self):
        m_domain = gpd.read_file(join(self.gwfm_dir, 'DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp'))
        xul, yul = list(m_domain.geometry.values[0].exterior.coords)[1]
        return xul, yul

    def define_model_attributes(self, ss_bool, ss_strt, strt_date, end_date):
        self.ss_bool = ss_bool
        self.ss_strt = pd.to_datetime(ss_strt)
        self.strt_date = pd.to_datetime(strt_date)
        self.end_date = pd.to_datetime(end_date)
        dates = pd.date_range(self.strt_date, self.end_date)
        self.nper = len(dates)
        if self.ss_bool:
            self.nper += 1

        self.perlen = np.ones(self.nper)
        if self.ss_bool:
            self.perlen[0] = 1/86400
        self.steady = np.zeros(self.nper)
        if self.ss_bool:
            self.steady[0] = 1
        self.steady = self.steady.astype('bool').tolist()
        self.nstp = np.ones(self.nper)

    def configure_paths(self, base_dir='C:/WRDAPP', model_nam = 'historical_simple_geology', scenario='reconnection'):

        loadpth = join(base_dir, 'GWFlowModel/Cosumnes/Regional/')
        model_ws = join(loadpth, '')
        if scenario == 'reconnection':
            model_ws += '_' + scenario
        # add dates to work space to mark difference in runs
        model_ws += '_' + str(self.strt_date.year) + '_' + str(self.end_date.year)
        if self.nper <= 2:
            model_ws = join(loadpth, 'steadystate')
        self.model_ws = model_ws
        os.makedirs(join(self.model_ws, 'input_data'), exist_ok=True)
    
class GroundwaterModel:
    def __init__(self, parent):
        self.parent = parent
        self.mf = None

    def create_model(self, modelname='MF', exe_name='mf-owhm.exe', version='mfnwt'):
        self.mf = flopy.modflow.Modflow(modelname=modelname, exe_name=exe_name, version=version, model_ws=self.parent.model_ws)

    def configure_discretization(self, nrow, ncol, nlay, delr, delc, rotation, xul, yul, proj4_str):
        self.dis = flopy.modflow.ModflowDis(
            model=self.mf, nrow=nrow, ncol=ncol, nlay=nlay, delr=delr, delc=delc,
            lenuni=2, itmuni=4, xul=xul, yul=yul, rotation=rotation, proj4_str=proj4_str,
            nper=self.parent.nper, perlen=self.parent.perlen, nstp=self.parent.nstp, steady=self.parent.steady,
            start_datetime=self.parent.strt_date
        )

    def load_elevation_data(self, gwfm_dir):
        dem_data = np.loadtxt(join(gwfm_dir, 'DIS_data/dem_52_9_200m_mean.tsv'))
        np.savetxt(join(self.parent.model_ws, 'input_data', 'dem_data.txt'), dem_data)
        return dem_data

    def configure_ibound(self, top_botm, dem_data):
        ibound = define_deep_aquifer_layering(top_botm, dem_data, cutoff_elev=56)
        # copy the ibound array to alter the geology array to set these cells as low permeability formations
        # either marine or volcanic based
        deep_geology = np.invert(ibound[:,:,:].astype(bool))
        np.savetxt(model_ws+'/input_data/deep_geology.tsv', np.reshape(deep_geology, (nlay*nrow,ncol)), delimiter ='\t')
        return deep_geology
    
    def set_ibound(self, botm, dem_data, hd_strt):
        strt = np.ones((self.dis.nlay, self.dis.nrow, self.dis.ncol), dtype=np.float32)
        if self.parent.ss_bool:
            strt[:, :, :] = self.dis.top[:, :]
        else:
            strt[:, :, :] = hd_strt

        self.ibound = np.ones((self.dis.nlay, self.dis.nrow, self.dis.ncol))
        self.ibound[botm > dem_data] = 0

        for k in range(self.dis.nlay):
            self.ibound[k, find_isolated_active(self.ibound[k])] = 0

        self.bas = flopy.modflow.ModflowBas(model=self.mf, ibound=self.ibound, strt=strt)

    def load_ghb_boundaries(self, gwfm_dir, strt_date, end_date, nrow):
        ghb_dir = join(gwfm_dir, 'GHB_data')
        strtyear = strt_date.year
        endyear = end_date.year + 1

        kriged_fall = np.zeros((int(endyear - strtyear), nrow, self.dis.ncol))
        kriged_spring = np.zeros((int(endyear - strtyear), nrow, self.dis.ncol))

        for t, year in enumerate(np.arange(strtyear, endyear)):
            spring_filename = glob.glob(join(ghb_dir, f'final_WSEL_arrays/spring{year}_kriged_WSEL.tsv'))[0]
            kriged_spring[t, :, :] = np.loadtxt(spring_filename) * 0.3048
            fall_filename = glob.glob(join(ghb_dir, f'final_WSEL_arrays/fall{year}_kriged_WSEL.tsv'))[0]
            kriged_fall[t, :, :] = np.loadtxt(fall_filename) * 0.3048

        return kriged_fall, kriged_spring



