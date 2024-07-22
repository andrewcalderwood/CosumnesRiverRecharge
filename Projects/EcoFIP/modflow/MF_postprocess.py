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

# %% [markdown]
# Load the modflow model output for all concepts and compare against the baseline scenario. 
# - Compare the water budget to review the impact of changes in farmland in the WEL/RCH 
# - compare the change in floodplain recharge and stream seepage
#

# %%
# standard python utilities
import os
from os.path import basename, dirname, join, exists, expanduser
import sys
from importlib import reload
import glob

import pandas as pd
import numpy as np
import time

# standard python plotting utilities
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.dates as mdates

# standard geospatial python utilities
import geopandas as gpd
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator


# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = join(usr_dir,'Box/research_cosumnes/GWFlowModel')

bas_dir = join(gwfm_dir, 'BAS6')
# proj_dir = join(gwfm_dir,'Projects','EcoFIP')


# %%
lwa_dir = join(usr_dir, 'LWA Dropbox','01_Project-Teams')
proj_dir = join(lwa_dir, '669.03 - DWR Cosumnes Floodplain Recharge')
main_concept_dir = join(proj_dir, 'Concepts')
gis_dir = join(main_concept_dir,'GIS')

fig_dir = join(main_concept_dir,'figures')


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
from map_cln import gdf_bnds, plt_cln
from mf_utility import get_dates


# %%
# load the baseline groundwater model for reference
# each concept will be loaded for just output files

run_dir = 'C://WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'

model_ws = loadpth+model_nam


# %%
load_only = ['DIS','BAS6','UPW','OC','SFR','LAK',
            'RCH', 'WEL'
            ]
m = flopy.modflow.Modflow.load('MF.nam', model_ws= model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only
                              )

if 'LPF' in m.get_package_list():
    gel_nam = 'LPF'
else:
    gel_nam = 'UPW'
gel = m.__getattr__(gel_nam)

# %%
# load datetime reference for the model
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')

# %%
base_model_ws = join(loadpth, 'EcoFIP')

# list the modflow model concept workspaces
concepts = os.listdir(base_model_ws)


# %% [markdown]
# # Water Budget plotting

# %%
import plot_utilities
reload(plot_utilities)

from plot_utilities import plt_wb, plt_wb_diff, sig_fill, get_lak_head

from mf_utility import clean_wb

# %%
wb0, out_cols, in_cols = clean_wb(model_ws, dt_ref)


# %%
n= 0
concept_ws = join(base_model_ws, concepts[n])
wb, out_cols, in_cols = clean_wb(concept_ws, dt_ref)

# %%
label_restoration = concepts[n]
label_baseline = 'Baseline'

# %%
# we can't plot the RIV_IN because it isn't in baseline

plt_cols = ['dSTORAGE_sum','LAK_IN', 'ET_OUT','GHB_NET', 'SFR_IN', 'SFR_OUT']
plt_labels=['Cumulative\nStorage Change', 'Floodplain\nRecharge', 'GW ET',
            'Net\nGW Flow','Stream\nLosses', 'Stream\nBaseflow']


# %%
fig,ax= plt.subplots(len(plt_cols),1, sharex=True,  figsize=(6.5, len(plt_cols)*1),dpi=300)
# plt_wb(wb.resample('MS').sum(), wb0.resample('MS').sum(), plt_cols, plt_labels, ax)
plt_wb(wb.resample('MS').mean(numeric_only=True), wb0.resample('MS').mean(numeric_only=True),
       plt_cols, plt_labels, ax, scale=1E-3)

# for n, var in enumerate(plt_cols):
#     sig_fill(ttest_all[ttest_all.term==var], ax[n])
    # sig_fill(ttest1_all[ttest1_all.term==var], ax[n])

fig.legend([label_restoration,label_baseline], ncol=2, loc='outside upper center', bbox_to_anchor=(0.5, 1.05),)
# fig.supylabel('Flux (MCM)')
fig.supylabel('Flux (thousand $m^3$/day)')
fig.tight_layout(h_pad=0.1)

# %%
# also want to review the floodplain recharge which only exists in the 
# concept version
riv_plt = wb.resample('MS').mean(numeric_only=True).RIV_IN.multiply(scale)
fig,ax =plt.subplots(figsize=(6,3))
riv_plt.plot(ax=ax)

# %% [markdown]
# Floodplain recharge (RIV_IN) shows declining rates with time likely as more water is recharged there is a reduced gradient. However, the rates are about 1/5 of the total stream seepage so it may be slightly overestimating recharge.

# %%
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
dif_lgd = [
    # Patch(facecolor='tab:blue', alpha=0.5, label='Reconnected Floodplain'),
    Line2D([0], [0],color='black',label='Difference'),
    # Line2D([0], [0], color='grey', label='Difference $\pm 1\sigma$'),
    Patch(facecolor='grey', alpha=0.5, label='Difference $\pm 1\sigma$'),

]

# %%
fig,ax= plt.subplots(len(plt_cols),1, sharex=True,  figsize=(6.5, len(plt_cols)*1),dpi=300)
scale = 1E-3

wb_diff = wb.resample('MS').mean(numeric_only=True) - wb0.resample('MS').mean(numeric_only=True)
wb_std = (wb[plt_cols]-wb0[plt_cols]).resample('MS').std(numeric_only=True)

plt_wb_diff(wb_diff, plt_cols, plt_labels, ax, color='black')

# plt_wb_diff(wb_diff+wb_std, plt_cols, plt_labels, ax, color='gray')
# plt_wb_diff(wb_diff-wb_std, plt_cols, plt_labels, ax, color='gray')
for n, var in enumerate(plt_cols):
    ax[n].fill_between(wb_diff.index, (wb_diff-wb_std)[var].multiply(scale),
                       (wb_diff+wb_std)[var].multiply(scale), color='gray', alpha=0.3)

# for n, var in enumerate(plt_cols):
#     # sig_fill(ttest_all[ttest_all.term==var], ax[n])
#     sig_fill(ttest_all[ttest_all.term==var], ax[n])

# fig.legend(['Difference','Std Dev'], ncol=2, loc='outside upper center', bbox_to_anchor=(0.5, 1.05),)
fig.legend(handles = dif_lgd, ncol=2, loc='outside upper center', bbox_to_anchor=(0.5, 1.05),)
#     ax[0].legend(['No Reconnection',label_baseline], ncol=2)
fig.supylabel('Difference in Flux (thousand $m^3$/day)')
fig.tight_layout(h_pad=0.1)


# %% [markdown]
# # SFR plotting
# - flow-duration curve similar to previous work on levee removal, hopefully more noticeable difference in low-flows. Not sure how to integrate change in peak flows (post-process SFR output by subtracting seepage from river package by cell?)
# - stream profile to show the change in groundwater elevation and seepage between concepts.

# %%
import mf_utility
from importlib import reload
reload(mf_utility)
from mf_utility import clean_sfr_df

# %%
# baseline sfr output
sfrdf0=  clean_sfr_df(model_ws, dt_ref, name='MF') #pd_sfr, 
# concepts
sfrdf=  clean_sfr_df(concept_ws, dt_ref, name='MF') #pd_sfr, 


# %%
# identify the sfrdf data where the floodplain is
riv_sfrdf0 = sfrdf0.reset_index().merge(riv_df[['i','j']], left_on=['row','column'], right_on=['i','j'])
# interest in summing the river leakage to compare with leakage from the floodplain
sfr_leak_sum0 = riv_sfrdf0.groupby('dt')['Qaquifer'].sum()
sfr_stage_mean0 = riv_sfrdf0.groupby('dt')['stage'].mean()
# sfrdf0

# %%
# riv_sfr_sum.plot()

# %%
# wb.resample('MS').mean().RIV_IN.plot()

# %%
# sfrdf0

# %%
ax.set_yscale('log')
# ax.set_xlim(.7,.8)
# ax.set_xlim(0, .7)
# ax.set_xlim(0.55,0.65)
ax.set_ylabel('Discharge ($m^3/d$)')
ax.set_xlabel('Exceedance Probability')
# plt_flow_duration(sfrdf0, sfrdf, [2015], ax)
fig.legend([label_restoration,label_baseline], ncol=2, loc='outside upper center', bbox_to_anchor=(0.5, 1.0),)


# %% [markdown]
# # Change in conditions at floodplain
# Sample groundwater elevations where there are RIV cells
# Plot RIV_IN and the stage assigned to the floodplain

# %%
hdobj0 = flopy.utils.HeadFile(model_ws+'/MF.hds')
hdobj = flopy.utils.HeadFile(concept_ws+'/MF.hds')

# %%
# get the recharge rate for the floodplain
riv_plt = wb.resample('MS').mean(numeric_only=True).RIV_IN


# %%
# floodplain cell locations for sampling groundwater levels
riv_df = pd.read_csv(join(concept_ws, 'river_ijk.csv'),index_col=0)
riv_idx = list(zip(riv_df.k, riv_df.i, riv_df.j))
# input floodplain stage
spd_stage = pd.read_csv(join(concept_ws,'river_stage.csv'))

# %%
riv_head0 = get_lak_head(hdobj0, riv_idx)
riv_head = get_lak_head(hdobj, riv_idx)

# %%

# %%
cols = ['SW-Inflw', 'Stage(H)','GW-Outflw']
labels=['SW Inflow \n(million $m^3/d$)', 'Floodplain\nInundation\nStage (m)', 'GW Recharge \n(million $m^3/d$)']

# scale = 1E6
fig,ax = plt.subplots(len(cols)+1,1, figsize=(6.5,6.5), dpi=300, sharex=True)

ax[1].plot(wb.RIV_IN.multiply(1E-3).values) # concept
ax[1].plot(sfr_leak_sum0.multiply(1E-3).values) # baseline
ax[1].set_ylabel('Floodplain Recharge\n(thousand $m^3/day$)')

# ax[2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
ax[-2].plot(spd_stage.flow_scale.values)
ax[-2].plot(sfr_stage_mean0.values)
ax[2].set_ylabel('Floodplain\nStage (m)')

ax[-1].plot(riv_head['value'].values,alpha=0.7)
ax[-1].plot(riv_head0['value'].values)
ax[-1].set_ylabel('GW Elevation (m)')
ax[1].legend([label_restoration,label_baseline], ncol=2, loc='upper right')
fig.tight_layout(h_pad=-0.1)
plt.xlabel('Datetime')

# %% [markdown]
# - The test concept of Blodgett Dam shows the concept floodplain stage is way higher than the pre-existing river stage which isn't possible. The floodplain stage should not be greater than 27 on average, this will lower max stage by 7 m.
# - Adjusting the max stage to 27 m produced more reasonable results where stage doesn't exceed river stage and we don't see as extreme recharge and rise in GWL. Likely we should remove the log scale calculation as it allows the floodplain to remain elevated for longer.
