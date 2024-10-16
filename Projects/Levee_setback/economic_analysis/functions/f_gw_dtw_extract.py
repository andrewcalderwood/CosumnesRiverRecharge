# ---
# jupyter:
#   jupytext:
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

# %% [markdown]
# # Extract groundwater data

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import glob
import sys
import time
from importlib import reload
import h5py

import pandas as pd
import numpy as np

# standard geospatial python utilities
import shapely
import geopandas as gpd

import matplotlib.pyplot as plt

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

# # flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# # other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)
from mf_utility import get_layer_from_elev, param_load
from report_cln import base_round

# %%
# resampled ground surface elevation
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')

# %%
year = int(2015)
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Regional/'
model_ws = loadpth+'historical_simple_geology_reconnection'
model_ws = loadpth+'strhc1_scale'


# %%
# year = int(2015)

def get_dtw(year, model_ws):
# %%

    load_only = ['DIS','BAS6']
    m = flopy.modflow.Modflow.load('MF.nam', model_ws=model_ws, 
                                    exe_name='mf-owhm.exe', version='mfnwt', load_only=load_only)
    botm = np.copy(m.dis.botm.array)

    # %%
    # also need shapefile of pumping well locations for each parcel
    parcel_wells = gpd.read_file(join(gwfm_dir, 'WEL_data', 'parcels_to_wells', 'parcels_to_wells.shp'))
    frow = parcel_wells.row-1
    fcol = parcel_wells.column-1
    # # parcel_wells layers (make 1-based
    parcel_wells['layer'] = get_layer_from_elev(dem_data[frow,fcol] - parcel_wells.depth_m*0.9, botm[:, frow,fcol], m.dis.nlay) + 1
    # get elevation
    well_dem = dem_data[parcel_wells.row-1, parcel_wells.column-1]
    
    nfield = parcel_wells.shape[0]

    # %%
    # sample heads
    # model_ws_last = loadpth + 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date())
    # hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')
    hdobj = flopy.utils.HeadFile(model_ws + '/MF.hds')

    # %%
    # the latest season date is 12/31 and earliest is 11/15 the year prior
    # easy fix would be to sample head until 12/31 assuming we can let model run until then
    # m_strt = pd.to_datetime(str(year-1)+'-10-1')
    m_strt = pd.to_datetime(str(year)+'-1-1')
    
    m_end = m_strt + pd.DateOffset(years=1) - pd.DateOffset(days=1)#+pd.DateOffset(months=3)
    # m_end = all_run_dates.iloc[m_per+1].date
    nper = (m_end - m_strt).days+1

    # could sample heads at daily or monthly + linear interpolation steps
    steps = pd.date_range(m_strt, m_end, freq='MS')
    print(m_strt, m_end)

    # %%
    # it's going to be very slow to sample the time series for all wells for all stress periods and time steps
    # it does seem more efficient to sample on a monthly scale like a farmer would (irrigation dates) and interpolate
    well_kij = list(zip(parcel_wells.layer-1, parcel_wells.row-1, parcel_wells.column-1))

    # %%
    dtw_days = (steps - m_strt).days.values # may need to compensate for steady state if added (better to adjust start heads)
    dtw_days = np.hstack((dtw_days, [nper-1])) # add end
    dtw = np.zeros((nper, nfield))

    # dtw_days = np.hstack((irr_days,[nper-1])) + spd_strt # irr_days alternate
    for n, kper in enumerate(dtw_days):
    # for kper in np.arange(0,nper): # daily sampling
        head = hdobj.get_data(kstpkper = (0,kper))
        dtw[kper,:] = well_dem - head[parcel_wells.layer-1, parcel_wells.row-1, parcel_wells.column-1]
        if n >=1: # not necessary except for plotting
            # fill in between as farmer might expect linear change
            dtw[dtw_days[n-1]:kper+1,:] = np.linspace(dtw[dtw_days[n-1]], dtw[kper], kper-dtw_days[n-1]+1)
    # convert from meters to feet
    dtw /= 0.3048

    # %%
    # some wells have relatively static dtw and some show a decline
    # import matplotlib.pyplot as plt
    # plt.plot(dtw[:,2080]);
    plt.plot(dtw[:,::100]);
    plt.ylim(0,200)
    # dtw.shape

    # %%
    # save output dtw for the year
    print('Saving output')
    # np.savetxt(join(model_ws, 'field_SWB', 'dtw_ft_parcels_'+str(year)+'.txt'), dtw)

    # convert to dataframe
    dtw_df = pd.DataFrame(dtw, columns = parcel_wells.UniqueID.astype(int))
    dtw_df.index = pd.date_range(m_strt, m_end)
    dtw_df.index.name = 'dt'
    return dtw_df
    # dtw_df.to_csv(join(model_ws, 'field_SWB', 'dtw_ft_parcels_'+str(year)+'.csv'))


# %%

# %%
def sample_dtw(end_heads, botm):
    """
    For a simple function we will insert the heads 
    that have already been sampled
    INPUT:
    end_heads: array of nlay, nrow, ncol with head values for the end
        of a simulatoin period
    botm: modflow layer bottom elevations
    OUTPUT:
    well_dtw: depth to water for each well for each parcel for
        the given heads input
    """
    nlay, nrow,ncol = botm.shape


    # %%
    # also need shapefile of pumping well locations for each parcel
    parcel_wells = gpd.read_file(join(gwfm_dir, 'WEL_data', 'parcels_to_wells', 'parcels_to_wells.shp'))
    frow = parcel_wells.row-1
    fcol = parcel_wells.column-1
    # # parcel_wells layers (make 1-based
    parcel_wells['layer'] = get_layer_from_elev(dem_data[frow,fcol] - parcel_wells.depth_m*0.9, botm[:, frow,fcol], nlay) + 1
    # get elevation
    well_dem = dem_data[parcel_wells.row-1, parcel_wells.column-1]
    
    nfield = parcel_wells.shape[0]

    # %%
    # sample the groundwater elevation where it is needed for wells for irrigation
    well_head = end_heads[parcel_wells.layer-1, parcel_wells.row-1, parcel_wells.column-1]
    # calculate the depth to water at each well    
    well_dtw = well_dem - well_head
    # convert from meters to feet
    well_dtw /= 0.3048

    # return output as a dataframe
    parcel_wells_out =  parcel_wells[['UniqueID','row','column']].copy()
    parcel_wells_out['dtw_ft'] = well_dtw

    # %%
    return(parcel_wells_out)


# %%
def avg_heads(sp_last, hdobj, m):
    """
    Given an array of stress periods, sample the heads
    then take the average
    """
    nspd = sp_last.shape[0]
    # allocate an array for the heads
    end_heads = np.zeros((nspd, m.dis.nlay,m.dis.nrow,m.dis.ncol))
    # sample the heads for each stress period and insert in the array
    for n in np.arange(nspd):
        end_heads[n] = hdobj.get_data(sp_last[n])
    # take the average head value for the previous 7 days
    end_heads = end_heads.mean(axis=0)
    return(end_heads)


# %%
## simple representative DTW for linear steps from dtw_min to dtw_max
## with a 5 ft decline from June to December based on observed data

# get min and max dtw to reference
def calc_simple_dtw(well_dtw, year, 
                    dtw_step = 10, decline_total = 5):
    """
    Given the minimum and maximum DTW in the domain
    Create a series of DTW profiles to run to be representative
    """
    # old version with heads from continuous model run
    # dtw_min, dtw_max = dtw_df.mean(axis=0).quantile([0,1]).apply(lambda x: base_round(x, dtw_step))
    # new version with heads from auto-updated
    dtw_min, dtw_max = well_dtw.dtw_ft.quantile([0,1]).apply(lambda x: base_round(x, dtw_step))
    
    # calculate dates for DTW decline, shouldn't really go longer than a year since it 
    # will be run at least once per year
    dtw_avg = pd.DataFrame(pd.date_range(str(year-1)+'-11-1', str(year)+'-12-31'), columns=['date'])
    dtw_avg = dtw_avg.assign(decline = 0).set_index('date')
    # dates where a decline date is specified
    decline_dates = dtw_avg.index[dtw_avg.index >=str(year)+'-6-1']
    decline = np.cumsum(np.full(len(decline_dates), decline_total/len(decline_dates)))
    dtw_avg.loc[decline_dates, 'decline'] = decline
    dtw_simple = np.repeat(np.reshape(np.arange(dtw_min, dtw_max, dtw_step), (1,-1)), len(dtw_avg), axis=0)
    dtw_simple = dtw_simple + np.reshape(dtw_avg.decline, (-1,1))
    dtw_simple_df = pd.DataFrame(dtw_simple, dtw_avg.index)
    return(dtw_simple_df)

