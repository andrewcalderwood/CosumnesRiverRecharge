# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Goal: take already copied/split up modflow model and copy to the scenario directory

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename
import sys
import glob
from importlib import reload
import shutil

import pandas as pd
import numpy as np
from scipy.stats import hmean, gmean

# import calendar
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import geopandas as gpd

# # mapping utilities
# import contextily as ctx
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as fm

# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
gwfm_dir


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

from mf_utility import get_dates, get_layer_from_elev, clean_wb, clean_hob
from map_cln import gdf_bnds
from report_cln import nse, base_round


# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
run_dir = 'C:/WRDAPP/GWFlowModel'
run_dir = 'F://WRDAPP/GWFlowModel'
# run_dir = 'D://WRDAPP/GWFlowModel'

# loadpth = run_dir +'/Cosumnes/levee_setback/streamflow/'
# # model_nam = 'setback_streamflow'
# model_nam = 'historical_streamflow'

loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology_reconnection'
# model_nam = 'input_write_2014_2020'
model_nam = 'input_write_2014_2022'

base_model_ws = loadpth+model_nam

# %%
load_only = ['DIS','BAS6','SFR']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only,
                              )

nrow,ncol,nlay,delr,delc = (m.dis.nrow, m.dis.ncol, m.dis.nlay, m.dis.delr, m.dis.delc)

# if 'LPF' in m.get_package_list():
#     gel_nam = 'LPF'
# else:
#     gel_nam = 'UPW'
# gel = m.__getattr__(gel_nam)

# %%
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')


# %%
# I've only ran the 2014-2020 model now, switch to 2014-2022 later

wb, out_cols, in_cols = clean_wb(base_model_ws, dt_ref)
# # clean_wb?

# %%
def overview_plt(wb):
    fig,ax= plt.subplots(3,1, sharex=True)
    ax2 = ax[0].twinx()
    ax2.set_ylim(0,5)
    wb.plot(y='PERCENT_ERROR', ax=ax2)
    ax2.set_ylabel('Percent Error')
    wb.plot(y=['dSTORAGE_sum'], legend=False, ax=ax[0], color='black')
    ax[0].set_ylabel('Cumulative Storage \nChange ($m^3$)')
    
    wb.plot(y=out_cols, ax=ax[1], legend=True)
    wb.plot(y=in_cols, ax=ax[2], legend=True)



# %%
overview_plt(wb)


# %% [markdown]
# # Iterate over different economic runs
# There isn't output for the flow budget from OWHM or hob out because I turned them off thinking I wouldn't need to check them.  
# - turns out it was just because the 2014-2022 base run didn't have the manual edit to the bas file. Should try to add this for future runs.
#
# The models take too long to to turn back on and re-run so I should plan to just run zonebudget with the cbc file for periods available and check heads at select locations. Maybe start by sampling heads because this might be faster if I choose just a few representative points. -> look at code that extracts simulated time series for HOB wells in oneto denier stream seepage

# %%
loadpth = run_dir +'/Cosumnes/Economic/'
model_nam = 'input_write_2014_2022'
scen = 'R204'
# it is probably better to create a slightly different file name then to copy these over for a set scenario
econ_model_ws = join(loadpth, model_nam+'_'+scen, 'crop_modflow')

all_run_dates = pd.read_csv(join(econ_model_ws, 'all_run_dates.csv'))



# %%
# load_only = ['DIS']
# m_per = flopy.modflow.Modflow.load('MF.nam', model_ws= join(econ_model_ws,d), 
#                                 exe_name='mf-owhm', version='mfnwt',
#                               load_only = load_only,
#                               )
# it looks like the start date time isn't correct in the economic model ws so can't really use this
# strt_date_per, end_date_per, dt_ref_per2 = get_dates(m_per.dis, ref='strt')


# %%
wb_per_all = pd.DataFrame()
for n, d in enumerate(all_run_dates.date[:-3]):
    # create locally referenced dt_ref file to avoid reloading dis
    dt_ref_per = dt_ref[(dt_ref.dt<=all_run_dates.date[n+1])&(dt_ref.dt>=d)].copy()
    dt_ref_per.kstpkper = list(zip(np.zeros(len(dt_ref_per),dtype=int), np.arange(len(dt_ref_per))))

    econ_ws_yr =join(econ_model_ws, d)
    wb_per, out_cols_per, in_cols_per  = clean_wb(econ_ws_yr, dt_ref_per)
    wb_per_all = pd.concat((wb_per_all, wb_per))
# econ_ws_yr

# %%
overview_plt(wb_per_all)


# %%
wb_per_all_ann = wb_per_all.resample('YS-OCT').sum()
wb_ann = wb.resample('YS-OCT').sum()


# %%
# overview_plt(wb_per_all_ann)
# overview_plt(wb_ann)

# should create overlay plots of major components, WEL, GHB_IN, RCH_IN, SFR_IN
# high level overview shows order of magnitude is similar after accounting for difference in recharge due to scenario

# %% [markdown]
# ## sample head time series
#
# change in storage in baseline shows some gradual decline with upswings from wet periods. need to cross reference with heads between scenarios to see if the scenario has the same pattern. Want to see how much changes in DTW are causing predicted crop to vary from observed.
#
# Take the hob package from the original run then sample the continuous time series for select wells fromt he original run and baseline model connect run.

# %%
def nse(targets,predictions):
    return 1-(np.sum((targets-predictions)**2)/np.sum((targets-np.mean(predictions))**2))

def clean_hob(model_ws):
    hobout = pd.read_csv(join(model_ws,'MF.hob.out'),delimiter=r'\s+', header = 0,names = ['sim_val','obs_val','obs_nam'],
                         dtype = {'sim_val':float,'obs_val':float,'obs_nam':object})
    hobout[['Sensor', 'spd']] = hobout.obs_nam.str.split('p',n=2, expand=True)
    hobout['kstpkper'] = list(zip(np.full(len(hobout),0), hobout.spd.astype(int)))
    hobout = hobout.join(dt_ref.set_index('kstpkper'), on='kstpkper')
    hobout.loc[hobout.sim_val.isin([-1e30, -999.99, -9999]), 'sim_val'] = np.nan
    hobout = hobout.dropna(subset='sim_val')
    hobout['error'] = hobout.obs_val - hobout.sim_val
    hobout['sq_error'] = hobout.error**2
    
    return(hobout)


# %%
obs_df = pd.read_csv(join(base_model_ws, 'input_data','all_obs_grid_prepared.csv'))

# obs_df.groupby('site_code')['gwe'].count()
obs_sites = obs_df.site_code.unique()
obs_sites.shape

# %%
hobout = clean_hob(base_model_ws, dt_ref, split_c='.')
# hobout

# %%
# # find where lake existed
# lak_lay, lak_row, lak_col = np.where(lakarr[0]==1)
# lak_kij = pd.DataFrame(np.transpose(np.where(m.lak.lakarr.array[0]==1)), columns=['k','i','j'])
# # get first layer below lake cells
# lak_kij = (lak_kij.groupby(['i','j']).max()+1).reset_index()
# # create tuples for sampling
# lak_idx = list(zip(lak_kij.k, lak_kij.i, lak_kij.j))

# %%
def get_lak_head(hdobj, lak_idx):
    """
    Return the spatially averaged head for the maximum head at the input locations (idx)
    hdobj: flopy head object
    idx: list of tuples as (layer, row, column)
    """
    # get heads under the lake
    lak_ts = hdobj.get_ts(lak_idx)
    lak_ts_df = pd.DataFrame(lak_ts, columns=['totim']+lak_idx)
    lak_ts_df = lak_ts_df.set_index('totim')
    lak_ts_df = lak_ts_df.melt(ignore_index=False)
    lak_ts_df[['k','i','j']] = lak_ts_df.variable.tolist()
    lak_ts_df = lak_ts_df.drop(columns='variable') # drop to speed up groupby
    lak_head = lak_ts_df.groupby(['totim','i','j']).max().groupby('totim').mean()
    return lak_head
