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
# Load in segmented modflow input to check recharge package was updated properly

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
# import pyproj # for converting proj4string
# import shapely
import geopandas as gpd
# import rasterio

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

from mf_utility import get_layer_from_elev
from map_cln import gdf_bnds



# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
run_dir = 'C:/WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
# run_dir = 'D://WRDAPP/GWFlowModel'


loadpth = run_dir +'/Cosumnes/Regional/'

model_nam = 'historical_simple_geology_reconnection'
model_nam = 'input_write_2014_2020'
model_nam = 'input_write_2014_2022'

base_model_ws = loadpth+model_nam


# %%
loadpth = run_dir +'/Cosumnes/Economic/'
model_nam = 'input_write_2014_2020'
model_nam = 'input_write_2014_2022'

model_ws = join(loadpth, model_nam, 'crop_modflow')


scen_nam = 'input_write_2014_2022_R3'
scen_ws = join(loadpth, scen_nam, 'crop_modflow')


# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')


# %%
# dem_data = np.copy(m.dis.top.array)
# botm = np.copy(m.dis.botm.array)

# %%
# For the tab files the left column is time (in model units) and the right column is flow (model units)
# Time is days, flow is cubic meters per day
# USGS presents flow in cfs (cubic feet per second)
inflow = pd.read_csv(join(gwfm_dir, 'SFR_data', 'MB_daily_flow_cfs.csv'), index_col = 'datetime', parse_dates = True)
# covnert flow from cubic feet per second to cubic meters per day
inflow['flow_cmd'] = inflow.flow_cfs * (86400/(3.28**3))



# %%
# deer creek doesn't flow in dry-season
# have minimum flow of 100/200 cfs to start deer creek flowing
# deer creek is approximated as about 10% of Cosumnes flow during wet season
dc_min_flow = 100*(86400/(3.28**3))
inflow_dc = inflow.copy()
# anything above the flow threshold is 10% of Cosumnes
inflow_dc.flow_cmd *= 0.1
# anything below the flow threshold is 0
inflow_dc.loc[inflow_dc.flow_cmd<dc_min_flow,'flow_cmd'] = 0

# %%
import h5py
uzf_dir = join(gwfm_dir, 'UZF_data')
nrow_p, ncol_p = 100,230


# %%
wel_dir = join(gwfm_dir, 'WEL_data')
uzf_dir = join(gwfm_dir, 'UZF_data')

# %% [markdown]
# Write static modflow files into the main directory including LPF, GHB, CHD. LPF (21 MB) with LAK (10 MB) and SFR (40 KB) will not need to be written as there is no dependence on stress periods. GHB, CHD, EVT, LAK, and SFR (20 MB) will need to be overwritten or saved multiple times as they have a change due to stress periods with ITMP. Pre-processing and writing output for each of these will save runtime later, but take up about 1.5 GB of storage.
#
#
# The RCH package is 492 MB and well package is 1.59 GB but these file sizes will be subdivided for each period so won't take up much more storage than before.

# %%
all_run_dates = pd.read_csv(join(model_ws, 'all_run_dates.csv'),parse_dates=['date'])
# all_run_dates

# %% [markdown]
# Load input to check key changes

# %%
load_only=['DIS','RCH']

# %%
nrow,ncol = grid_p.row.max(), grid_p.column.max()



# %%
nruns = all_run_dates.shape[0]-1

# %%
##############################################################################################
## write out the irrigation independent inputs (GHB, CHD, UPW, OC, NWT, DIS)
print('Reading  input files')
rech_all = np.full((0, nrow,ncol),np.nan)
rech_sum = np.full((nruns, nrow,ncol),np.nan)

# for m_per in [0]:
for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date

    m_model_ws = join(model_ws,str(m_strt.date()))
    print(basename(m_model_ws), end=',\t')
    # switch to modflow nwt to enable option bloack for use in owhm
    m_month = flopy.modflow.Modflow.load(join(m_model_ws,'MF.nam'),
                                        load_only=load_only,
                                        )

    rech = m_month.rch.rech.array[:,0]
    rech_all = np.vstack((rech_all, rech))
    rech_sum[m_per] = rech.sum(axis=0)


# %%
rech_all_s = np.full((0, nrow,ncol),np.nan)
rech_sum_s = np.full((nruns, nrow,ncol),np.nan)

for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date

    m_model_ws = join(scen_ws,str(m_strt.date()))
    # switch to modflow nwt to enable option bloack for use in owhm
    m_month_s = flopy.modflow.Modflow.load(join(m_model_ws,'MF.nam'),
                                        load_only=load_only,
                                        )

    rech = m_month_s.rch.rech.array[:,0]
    rech_all_s = np.vstack((rech_all_s, rech))
    rech_sum_s[m_per] = rech.sum(axis=0)

# %%
# plt_temp = rech_sum.sum(axis=(1,2))[1:]

# %%

# %%
# ny=3
# fig,ax = plt.subplots(ny,1, figsize=(4,4), sharex=True, dpi=300)
# for n in np.arange(0,ny):
#     im = ax[n].imshow(rech_sum[n])
#     plt.colorbar(im, shrink=0.8)

# %%
rech_sum.sum(axis=(1,2))
# rech_sum_s.sum(axis=(1,2))

# %%
# scale by number of cells to get distributed recharge
area = 100*230
# scale=
rech_plt_s = rech_sum_s.sum(axis=(1,2))[1:]/area
plt.plot(rech_plt_s, label='Scenario')

rech_plt = rech_sum.sum(axis=(1,2))[1:]/area
plt.plot(rech_plt, label='Baseline')
# plt.plot(plt_temp/area, label='Baseline')
# plt.plot(rech_sum.sum(axis=(1,2))[1:]/area, label='Scenario')
plt.legend()
dx=4
plt.xticks(np.arange(0, len(plt_temp))[::dx], all_run_dates.date.dt.year[1:-1][::dx]);
plt.ylabel('Annual total recharge rate (m)')

# %%
