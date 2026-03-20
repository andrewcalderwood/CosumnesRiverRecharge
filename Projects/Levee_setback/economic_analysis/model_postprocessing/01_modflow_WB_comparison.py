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
dt_ref['totim'] = (dt_ref.dt-dt_ref.dt.iloc[0]).dt.days+1


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
# overview_plt(wb)
overview_plt(wb.resample('MS').sum())


# %% [markdown]
# # Iterate over different economic runs
#
# Goal: 
# - check change in storage is consistent, i.e., overall mass balance is consistent
# - check that pumping dynamics are consistent on seasonal and annual basis, timing will be different due to daily vs timed irrigation periods
# - recharge may be the most affected due to larger irrigation events and also because the SWB is different between native and non-native and that irrigated lands have a simpler SWB run in the winter. Currently soil moisture is reset to field capacity (I think, maybe wilting point) instead of passing the actual soil moisture from previous run end.
#
#

# %%
loadpth = run_dir +'/Cosumnes/Economic/'
model_nam = 'input_write_2014_2022'
scen = 'R200'
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
for n, d in enumerate(all_run_dates.date[:-1]):
    # create locally referenced dt_ref file to avoid reloading dis
    dt_ref_per = dt_ref[(dt_ref.dt<=all_run_dates.date[n+1])&(dt_ref.dt>=d)].copy()
    dt_ref_per.kstpkper = list(zip(np.zeros(len(dt_ref_per),dtype=int), np.arange(len(dt_ref_per))))

    econ_ws_yr =join(econ_model_ws, d)
    wb_per, out_cols_per, in_cols_per  = clean_wb(econ_ws_yr, dt_ref_per)
    wb_per_all = pd.concat((wb_per_all, wb_per))
# econ_ws_yr

# %%
# overview_plt(wb_per_all)
overview_plt(wb_per_all.resample('MS').sum())


# %%
wb_per_all_ann = wb_per_all.resample('YS-OCT').sum()
wb_ann = wb.resample('YS-OCT').sum()

wb_per_all_mon = wb_per_all.resample('MS').sum()
wb_mon = wb.resample('MS').sum()

# %%
overview_plt(wb_per_all_ann)
overview_plt(wb_ann)

# should create overlay plots of major components, WEL, GHB_IN, RCH_IN, SFR_IN
# high level overview shows order of magnitude is similar after accounting for difference in recharge due to scenario

# %%
fig,ax = plt.subplots(3,1)

for n, p in enumerate(['RCH_IN','WEL_OUT', 'SFR_IN']):
    ax_n = ax[n]
    wb_per_all_mon.plot(y=p, ax=ax_n, label='Connected')
    wb_mon.plot(y=p, ax=ax_n, label='Original')
    ax_n.set_ylabel(p)

# there is a significant drop in recharge estimates in the wet season
# need to evaluate if due to loss of soil moisture passing or from more efficient irrigation

# %% [markdown]
# After re-running with fixed dtw for crop choice pumping is still fairly close even with more fields which is odd but maybe the added fields were lower usage? Recharge is still too low.

# %% [markdown]
# ## Do more direct comparison of pumping and recharge
#
# Want to see if pumping is fairly close and investigate why recharge seems slightly lower. Is this because of the change in irrigation style or because the soil moisture is not passed between runs so it resets and that has to refill before good percolation? Best practice might be to simply add in a way to assign input soil moisture and create simple csv style to save output for each in a subfolder.

# %%

# %%

# %% [markdown]
# # sample head time series
#
# change in storage in baseline shows some gradual decline with upswings from wet periods. need to cross reference with heads between scenarios to see if the scenario has the same pattern. Want to see how much changes in DTW are causing predicted crop to vary from observed.
#
# Take the hob package from the original run then sample the continuous time series for select wells fromt he original run and baseline model connect run.

# %%
# load dataset with pre-processed info
obs_df = pd.read_csv(join(base_model_ws, 'input_data','all_obs_grid_prepared.csv'))
# find information on just sites
wlm_cols =  obs_df.columns[obs_df.columns.str.contains('wlm')].tolist()
sites_df = obs_df.drop_duplicates('site_code').drop(columns=['date','gwe','obs_nam']+wlm_cols)
obs_sites = sites_df.site_code.unique()
sites_df.shape

# 1 based row, columns
# get 0-based index
sites_df[['i','j','k']] = sites_df[['row','column','layer']] - 1
# create tuple for time series sampling
sites_idx = list(zip(sites_df.k, sites_df.i, sites_df.j))

len(sites_idx)

# %%
hobout = clean_hob(base_model_ws, dt_ref, split_c='.')
# hobout

# %%
from hyd_utility import hob_fit_scatter

# %%
# hob_fit_scatter(hobout)

# %%
dt_ref_mon = dt_ref[dt_ref.dt.dt.day==1]
mon_kper = dt_ref_mon.kstpkper.tolist()

# %%
# pretty slow to load (<1 min though)
hdobj = flopy.utils.HeadFile(base_model_ws+'/MF.hds')
spd_stp = hdobj.get_kstpkper()
times = hdobj.get_times()
# cbc = base_model_ws+'/MF.cbc'

# %%
def get_ts_head(hdobj, cell_idx):
    """
    Return the spatially averaged head for the maximum head at the input locations (idx)
    hdobj: flopy head object
    idx: list of tuples as (layer, row, column)
    """
    # get heads under the lake
    hd_ts = hdobj.get_ts(cell_idx)
    hd_ts_df = pd.DataFrame(hd_ts, columns=['totim']+cell_idx)
    hd_ts_df = hd_ts_df.set_index('totim')
    hd_ts_df = hd_ts_df.melt(ignore_index=False)
    hd_ts_df[['k','i','j']] = hd_ts_df.variable.tolist()
    hd_ts_df = hd_ts_df.drop(columns='variable') # drop to speed up groupby
    # lak_head = lak_ts_df.groupby(['totim','i','j']).max().groupby('totim').mean()
    return hd_ts_df


# %%
# might be more effective to sample on a monthly scale if possible

# took 13 minutes
tic = time.time()
obs_ts_df = get_ts_head(hdobj, sites_idx)
toc = time.time()
print('Time was %.2f min' %((toc-tic)/60))
# add info on date
obs_ts_df = obs_ts_df.reset_index().merge(dt_ref)


# %%

obs_ts_df_all = pd.DataFrame()
for n, d in enumerate(all_run_dates.date[:-1]):
    # create locally referenced dt_ref file to avoid reloading dis
    dt_ref_per = dt_ref[(dt_ref.dt<=all_run_dates.date[n+1])&(dt_ref.dt>=d)].copy()
    dt_ref_per.kstpkper = list(zip(np.zeros(len(dt_ref_per),dtype=int), np.arange(len(dt_ref_per))))
    dt_ref_per['totim'] = (dt_ref_per.dt-dt_ref_per.dt.iloc[0]).dt.days+1
    
    econ_ws_yr =join(econ_model_ws, d)
    hdobj_yr = flopy.utils.HeadFile(econ_ws_yr+'/MF.hds')
    obs_ts_df_yr = get_ts_head(hdobj_yr, sites_idx)
    obs_ts_df_yr = obs_ts_df_yr.reset_index().merge(dt_ref_per)

    obs_ts_df_all = pd.concat((obs_ts_df_all, obs_ts_df_yr))

# %% [markdown]
# ## simulated head comparison
#
# Main goal: show connected model is on par with original model.
#
# I made a version with observations as well but this isn't directly helpful, it generally shows magnitude is met in most areas but we undersimulate on the lower Cosumnes floodplain area. Should look into the reason why the local model does this better than regional.

# %%

# %%
obs_comp = obs_ts_df_all.merge(sites_df[['i','j','k','site_code']])

# %%
# concat obs original and model connect baseline for comparison
obs_comp = pd.concat((obs_ts_df_all.assign(name='scenario'), 
           obs_ts_df.assign(name='original')))

obs_comp = obs_comp.merge(sites_df[['i','j','k','site_code']])

# could think about adding in obs heads as wells for comparison
obs_df_join = obs_df.rename(columns={'date':'dt','gwe':'value'})[['dt','value','site_code']]

obs_comp = pd.concat((obs_comp, obs_df_join.assign(name='observed')))

# %%
import seaborn as sns
# obs_comp

# %%
# break up hydrographs by region
sites_df.dem_elev.hist()
high_sites = sites_df[sites_df.dem_elev>=20].site_code
mid_sites = sites_df[(sites_df.dem_elev>=10)&(sites_df.dem_elev<20)].site_code
low_sites = sites_df[sites_df.dem_elev<10].site_code


# %%
sns.relplot(obs_comp[obs_comp.site_code.isin(low_sites)], x='dt', y='value', hue='name',col='site_code',col_wrap=4, kind='line')


# %%
sns.relplot(obs_comp[obs_comp.site_code.isin(mid_sites)], x='dt', y='value', hue='name',col='site_code',col_wrap=4, kind='line')


# %%
sns.relplot(obs_comp[obs_comp.site_code.isin(high_sites)], x='dt', y='value', hue='name',col='site_code',col_wrap=4, kind='line')


# %%
