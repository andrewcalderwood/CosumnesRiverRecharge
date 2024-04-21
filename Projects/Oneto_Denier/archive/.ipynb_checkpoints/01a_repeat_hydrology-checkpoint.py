# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
# standard python utilities
import os
import sys
from os.path import basename, dirname, join, exists, expanduser
import glob
import time
import shutil

import pandas as pd
import numpy as np
import numpy.ma as ma
from scipy.stats import gmean

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
# import pyproj # for converting proj4string
import geopandas as gpd
from osgeo import gdal
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %% editable=true slideshow={"slide_type": ""}
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study
proj_dir = gwfm_dir + '/Oneto_Denier/'
dat_dir = proj_dir+'Stream_level_data/'

sfr_dir = gwfm_dir+'/SFR_data/'
uzf_dir = gwfm_dir+'/UZF_data/'


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

# from mf_utility import get_layer_from_elev
from mf_utility import get_layer_from_elev, param_load

from map_cln import gdf_bnds, plt_cln

# %%
tprogs_fxn_dir = doc_dir+'/GitHub/CosumnesRiverRecharge/tprogs_utilities'
add_path(tprogs_fxn_dir)
import tprogs_cleaning as tc

from importlib import reload
reload(tc)
tprogs_info = [80, -80, 320]

# %%
ext_dir = 'F:/WRDAPP'
c_dir = 'C:/WRDAPP'

if os.path.exists(ext_dir):
    loadpth = ext_dir 
elif os.path.exists(c_dir):
    loadpth = c_dir 
    
loadpth += '/GWFlowModel/Cosumnes/Stream_seepage/'
model_nam = 'oneto_denier_upscale4x_2014_2020'
model_nam = 'oneto_denier_upscale4x_2014_2020_no_reconnection'

model_ws = loadpth+ model_nam

ws_out = join(loadpth, 'repeat_hydrology')


# %% [markdown]
# # Script goals
# 1. Load the existing model.
# 2. Define a new model with time DIS multiplying the previous nper by 5 to represent a repeated period.
# 3. Reuse spatial DIS and reuse exactly the packages BAS6, LPF, NWT
# 4. Repeat the inputs of WEL, RCH, SFR/tab, LAK, GHB 5 times
# 5. Extend the OC package

# %%
m0 = flopy.modflow.Modflow.load(join(model_ws, 'MF.nam'))

# %% [markdown]
# New discretization

# %%

# Extract spatial and temporal variables
dis0 = m0.dis
# may not help much
dis_sp_attr = ['nlay', 'nrow','ncol', 'delr', 'delc','top','botm','itmuni', 'lenuni']
dis_sp_val = list(getattr(dis0, a) for a in dis_sp_attr)
dis_tm_attr = ['nper', 'perlen', 'nstp','steady']
dis_tm_val = list(getattr(dis0, a) for a in dis_tm_attr)


# %%
ncycle = 5
nper = dis0.nper*ncycle
perlen = np.repeat(dis0.perlen.array, ncycle) 
nstp = np.repeat(dis0.nstp.array, ncycle)
steady = np.repeat(dis0.steady.array, ncycle)

# %%
strt_date = pd.to_datetime(dis0.start_datetime) +pd.DateOffset(days=dis0.nper)

# %%

m = flopy.modflow.Modflow(modelname = 'MF', exe_name = 'mf-owhm', 
                          version = 'mfnwt', model_ws=ws_out)

#lenuni = 1 is in ft, lenuni = 2 is in meters
# itmuni is time unit 5 = years, 4=days, 3 =hours, 2=minutes, 1=seconds
dis = flopy.modflow.ModflowDis(nrow=dis0.nrow, ncol=dis0.ncol, 
                               nlay=dis0.nlay, delr=dis0.delr, delc=dis0.delc,
                               model=m, lenuni = 2, itmuni = 4,
                               
                               # xul = xul, yul = yul,rotation=rotation, proj4_str=proj4_str,
                              nper = nper, perlen=perlen, nstp=nstp, steady = steady,
                              start_datetime = strt_date)

# %%
dis.top = dis0.top
dis.botm = dis0.botm

# %%
dis.nper

# %% [markdown]
# # Static packages
# Reuse directly the LPF, NWT. Reuse the BAS after updating starting heads

# %%
# get heads from last stress period as starting for next simulation year
hdobj = flopy.utils.HeadFile(model_ws+'/MF.hds')
strt_hds = hdobj.get_data((0,dis0.nper-1))

bas = flopy.modflow.ModflowBas.load(join(model_ws, 'MF.bas'), model = m)
bas.strt = strt_hds

# %%
upw = flopy.modflow.ModflowUpw.load(join(model_ws, 'MF.upw'), model = m)


# %%
nwt = flopy.modflow.ModflowNwt.load(join(model_ws, 'MF.nwt'), model = m)


# %% [markdown]
# # Variable packages
# Repeat the inputs
# - WEL, RCH, SFR/tab, LAK, GHB 5 times

# %% [markdown]
# ## WEL

# %%
# runs fairly quickly (few seconds)
# empty well dictionary
wel_dict = {}
# repeat for the set number of cycles
for c in np.arange(0, ncycle):
    # each cycle is a repeat of the original transient model
    for s in np.arange(0, dis0.nper):
        wel_dict[s+(c*dis0.nper)] = m0.wel.stress_period_data[s]

# %%
wel = flopy.modflow.ModflowWel(m, stress_period_data=wel_dict, ipakcb=55)

# %% [markdown]
# ## RCH

# %%
# need to open first to avoid over allocating 
rech = m0.rch.rech.array[:,0]

# %%
# opening the array by itself runs quickly
# t0 = time.time()
rch_dict = dict()
for c in np.arange(0, ncycle):
    for s in np.arange(0, dis0.nper):
        rch_dict[s+(c*dis0.nper)] = rech[s]

# t1 = time.time()


# %%
rch = flopy.modflow.ModflowRch(m, rech=rch_dict, nrchop = m0.rch.nrchop, ipakcb=55)

# %% [markdown]
# ## EVT

# %%
# need to open first to avoid over allocating 
evtr = m0.evt.evtr.array[:,0]
surf = m0.evt.surf.array[0,0]
exdp = m0.evt.exdp.array[0,0]
ievt = m0.evt.ievt.array[0,0]
nevtop = m0.evt.nevtop

# %%
# opening the array by itself runs quickly
# t0 = time.time()
evtr_dict = dict()
for c in np.arange(0, ncycle):
    for s in np.arange(0, dis0.nper):
        evtr_dict[s+(c*dis0.nper)] = evtr[s]

# t1 = time.time()


# %%

# surf : et surface elevation. evtr: maximum ET flux
# exdp: extinction depth. ievt : layer indicator variable
# nevtop = 3 -> highest active layer
# nevtop = 2 -> layer defined in ievt
evt = flopy.modflow.ModflowEvt(model=m, nevtop = nevtop, ievt = ievt, 
                               evtr = evtr_dict, exdp = exdp,  
                               surf = surf, ipakcb = 55)



# %%
# evt.write_file()

# %% [markdown]
# ## GHB

# %%
# runs fairly quickly (few seconds)
# empty well dictionary
ghb_dict = {}
# repeat for the set number of cycles
for c in np.arange(0, ncycle):
    # each cycle is a repeat of the original transient model
    for s in np.arange(0, dis0.nper):
        ghb_dict[s+(c*dis0.nper)] = m0.ghb.stress_period_data[s]

# %%
ghb = flopy.modflow.ModflowGhb(m, stress_period_data=ghb_dict, ipakcb=55)

# %% [markdown]
# ## SFR
# Floodplain flow threshold ($m^3/d$)  
# Baseline: 9417600.0  
# Restoration: 1987200.0  

# %%
# could reuse SFR as is, but need to update dataset 5
# m0.sfr

# %%
# sfr = type('sfr', m0.sfr.__bases__, dict(m0.sfr.__dict__))
sfr0 = m0.sfr

sfr = flopy.modflow.ModflowSfr2(model = m,
    nstrm = sfr0.nstrm, nss = sfr0.nss, const = sfr0.const, dleak = sfr0.dleak, isfropt = sfr0.isfropt, ipakcb = 55, istcb2 = 54,
    irtflg =sfr0.irtflg, numtim = sfr0.numtim, flwtol = sfr0.flwtol, weight = sfr0.weight,
    reachinput=True, transroute=True, tabfiles=True, tabfiles_dict={1: {'numval': nper, 'inuit': 56}},
    reach_data = sfr0.reach_data, segment_data = sfr0.segment_data[0], channel_geometry_data = sfr0.channel_geometry_data)

# only works on GitHub version with edits
# Add option block at the top of the sfr input file for tabfiles
options_line = ' reachinput transroute tabfiles 1 ' + str(nper) + ' no_reach_layer_change'
tab_option = flopy.utils.OptionBlock(options_line = options_line, package = sfr, block = True)
sfr.options = tab_option

# %%
tab = np.loadtxt(join(model_ws, 'MF.tab'))

# %%
tab1 = np.tile(tab[:,1], ncycle)
tab0 = np.tile(tab[:,0], ncycle)+ np.repeat(np.arange(0, nper, dis0.nper), dis0.nper)
np.savetxt(join(ws_out, 'MF.tab'), np.transpose((tab0, tab1)))

# %% [markdown]
# ## LAK

# %%
# lake stage (elevation), volume, and area (3 numbers per line)
bathtxt = np.loadtxt(model_ws+'/MF.bath',  delimiter = '\t')

stages = bathtxt[:,0].min() +0.01
# (ssmn, ssmx) max and min stage of each lake for steady state solution, there is a stage range for each lake
# so double array is necessary
stage_range = [[bathtxt[:,0].min(), bathtxt[:,0].max()]]

# %%
lakarr =  m0.lak.lakarr.array[0]
bdlknc = m0.lak.bdlknc.array[0]
flux_data = m0.lak.flux_data 

# %%
lak = flopy.modflow.ModflowLak(model = m, lakarr = lakarr, bdlknc = bdlknc,  stages=stages, 
                               stage_range=stage_range, flux_data = flux_data,
                               theta = 1, nssitr = 1000, sscncr = 1E-5, surfdepth = 0.02, # take from Shasta model
                               tabdata= True, tab_files='MF.bath', tab_units=[57],ipakcb=55)

lak.options = ['TABLEINPUT']
# # need to reset tabdata as True before writing output for LAK
lak.tabdata = True

# %% [markdown]
# ## OC

# %%
# get ALL stress periods and time steps list, not just those in the output
kstpkper = []
for n,stps in enumerate(m.dis.nstp.array):
    kstpkper += list(zip(np.arange(0,stps),np.full(stps,n)))

# %%
# Output control
# default unit number for heads is 51, cell by cell is 53 and drawdown is 52

# For later model runs when all the data is needed to be saved
oc_dict = { (sp, ts): ['save head'] for ts, sp in kstpkper}

oc = flopy.modflow.ModflowOc(model = m, stress_period_data = oc_dict, compact = True)

# %% [markdown]
# # Write out files

# %% [markdown]
# Writes out pretty slowly. It could be worth segmenting this into the five year blocks to make more manageable outputs.
#
# GHB might be the slowest because the final size is 3.6 GB while RCH is 1 GB. It might make more sense to use a simple long-term average for the GHB to save memory

# %%
# m.write_input()

# %%
# ghb.write_file()
# m.write_name_file()

# %%
files = ['nam', 'gage', 'bath', 'hob']
for f in files:
    f = 'MF.'+f
    shutil.copy(join(model_ws, f), join(m.model_ws, f))

# %%
# sfr.write_file()

# %%
# m.get_package_list()
