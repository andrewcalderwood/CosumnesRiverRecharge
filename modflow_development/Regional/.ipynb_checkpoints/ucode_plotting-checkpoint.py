# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
# standard python utilities
import os
from os.path import basename, dirname, join, exists
import sys
from importlib import reload
import glob
import pandas as pd
import numpy as np
import calendar
import time

# standard python plotting utilities
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# standard geospatial python utilities
import geopandas as gpd
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator

# import flopy
# import flopy.utils.binaryfile as bf
from importlib import reload


# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

bas_dir = join(gwfm_dir, 'BAS6')
proj_dir = join(gwfm_dir,'Mapping')
plt_dir = join(proj_dir,'figures/')


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
# from map_cln import gdf_bnds, plt_cln
import map_obs_plt as mop
from map_obs_plt import plt_bc_hk, plot_head_simple, plot_dtw_simple
from map_cln import gdf_bnds, plt_cln
from report_cln import base_round

from flopy_utilities import zone_clean, reach_data_gdf
from mf_utility import get_dates, clean_hob
# from mf_utility import clean_wb # not yet because use GW_IN instead

# %%
run_dir = 'C://WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'
base_model_ws = join(loadpth, model_nam)
model_nam = 'historical_geology_cal'

model_ws = loadpth+model_nam


# %%
load_only = ['DIS','BAS6','UPW','OC','SFR','LAK',
            # 'RCH',
            #  'WEL'
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
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')
# round of the steady state period
dt_ref['dt'] = dt_ref.dt.dt.round('D')
dt_ref = dt_ref[~dt_ref.steady]

# %%
# save base arrays for comparison
hk_base = gel.hk.array
ss_base = gel.ss.array
sy_base = gel.sy.array

# %%
gel_dict = {'hk':hk_base,'ss':ss_base,'sy':sy_base}

# %% [markdown]
# # Check parameters

# %%
r_params

# %%
# review folders
m_dirs = pd.Series(os.listdir(model_ws))
m_dirs = m_dirs[m_dirs.str.contains(r'r\d{3}')]

base_params = pd.read_csv(join(model_ws, 'ZonePropertiesInitial.csv'))
param_check = pd.DataFrame(columns=['Zone','Lithology', 'Variable','Value'])
for n in np.arange(0, len(m_dirs)):
    r_params = pd.read_csv(join(model_ws, m_dirs.iloc[n], 'ZonePropertiesInitial.csv'))
    np_diff = ((r_params[['K_m_s','Ss','Sy']]/base_params[['K_m_s','Ss','Sy']])>1)
    row, col = np.where(np_diff)
    p_diff = r_params[['K_m_s','Ss','Sy']].iloc[row,col].iloc[0]  
    p0_diff = base_params[['K_m_s','Ss','Sy']].iloc[row,col].iloc[0]  
    print(n, p_diff.index[0], p_diff.values[0], '--', p0_diff.values[0])
    # save values in dataframe for reference
    row_add = np.append(r_params.iloc[row][['Zone','Lithology']].values, [p_diff.index[0], p_diff.values[0]])
    param_check.loc[n] = row_add
# verifies there is at least one parameter adjusted by 1.01
# also shows that the perturbamnt was done in parameter space, not log space
# I think the steps are just too small to see an effect, still need to verify change in geology

# %%
# clean variable names for referencing
param_check['gel_var'] = param_check.Variable.str.lower()
param_check.gel_var = param_check.gel_var.str.replace('k_m_s','hk')

# %%
# iterate over geologic parameters to verify if changed
for n in np.arange(0, len(m_dirs)):
    upw = flopy.modflow.ModflowUpw.load(join(model_ws, m_dirs.iloc[n], 'MF.upw'), model=m)
    gel_var = param_check.loc[n].gel_var
    gel_arr = getattr(upw, gel_var).array/gel_dict[gel_var]
    print(n, param_check.loc[n].Lithology, gel_var, (gel_arr!=1).sum(), gel_arr[gel_arr!=1].mean())

# %% [markdown]
# The geologic arrays are being changed by the parameters, and typically have the correct 1.01 change.

# %%
# plot surface average
plt.imshow((gel_arr!=1).sum(axis=0))
plt.colorbar(shrink=0.6)
plt.show()
# plot cross section along columns
plt.imshow((gel_arr!=1).sum(axis=1), aspect=5)
plt.colorbar(shrink=0.6)
plt.show()



# %% [markdown]
# - spatially the change is distributed across all layers and rows larger impact on the sands/gravels
# - the change is also distributed across layers, although change looks more significant near surface

# %%
# plot surface average
plt.imshow(gel_arr.mean(axis=0))
plt.colorbar(shrink=0.6)
plt.show()

# %% [markdown]
# - plotting the mean change by layer shows difference up to 1.008 (~1.01) after accounting for the upscaling that impacts changed cells with those nearby

# %%
k=0
fig,ax = plt.subplots(2,1, sharex=True)
im = ax[0].imshow(getattr(upw, gel_var).array[k])
plt.colorbar(im)
im = ax[1].imshow(gel_dict[gel_var][k])
plt.colorbar(im)


# %%
all_obs = pd.read_csv(model_ws+'/input_data/all_obs_grid_prepared.csv',index_col=0, parse_dates=['date'])
all_obs.index = all_obs.index.rename('date')
all_obs = all_obs.reset_index()

# %%
hobout = clean_hob(model_ws, dt_ref, split_c = '.')


# %%
hobout_all_n = pd.DataFrame()
for n in np.arange(0,len(m_dirs)):
    hobout_n = clean_hob(join(model_ws,m_dirs.iloc[n]), dt_ref, split_c = '.').assign(r=n)
    hobout_all_n = pd.concat((hobout_all_n, hobout_n))
# add correct columns for dt
hobout_all_n= hobout_all_n.drop(columns=['spd','kstpkper','dt','wy']).merge(all_obs)

# %%

# %%
hobout_all_n.groupby(['r']).mean(numeric_only=True)

# %%
# all_obs

# %% [markdown]
# Looking at mean error for all sites by parallel run, it really shows the difference is on the 0.001-0.01 m scale.

# %%
# sns.relplot(hobout_all_n[hobout_all_n.Sensor=='N13407'],
#             x='dt',y='sim_val', hue='r')
sns.relplot(hobout_all_n,
            x='date',y='sim_val', hue='r', col_wrap=4, col='Sensor',
                      facet_kws={'sharex':True, 'sharey':False}
    )

# %% [markdown]
# We only see a shift in results for select wells, verify that geology change is across arrays, not just localized.
# - 2926, 10384, 10418, 10383, 10161, 11078, 14626 are examples that show a slight change
#     - it may be that these are in a certain type of geology?

# %%
chk = [2926, 10384, 10418, 10383, 10161, 11078, 14626]

# %%
hob_gpd = all_obs[all_obs.node.isin(chk)].drop_duplicates('node')
hob_gpd = gpd.GeoDataFrame(hob_gpd, geometry = gpd.points_from_xy(hob_gpd.longitude, hob_gpd.latitude), crs='epsg:4326')

# %%
# hob_gpd.plot()
