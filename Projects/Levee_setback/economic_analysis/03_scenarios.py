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
# Script to write out scenarios to run under 03_copy_model_modflow and then 03_model_connect.py

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
# import rasterio


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

# loadpth = run_dir +'/Cosumnes/levee_setback/streamflow/'
# # model_nam = 'setback_streamflow'
# model_nam = 'historical_streamflow'

loadpth = run_dir +'/Cosumnes/Regional/'

model_nam = 'historical_simple_geology_reconnection'
model_nam = 'input_write_2014_2020'

base_model_ws = loadpth+model_nam


# %%
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt')

nrow,ncol,nlay,delr,delc = (m.dis.nrow, m.dis.ncol, m.dis.nlay, m.dis.delr, m.dis.delc)

if 'LPF' in m.get_package_list():
    gel_nam = 'LPF'
else:
    gel_nam = 'UPW'
gel = m.__getattr__(gel_nam)

# %%
# loadpth = run_dir +'/Cosumnes/Economic/'

# m.model_ws = join(loadpth, model_nam, 'crop_modflow')
# # drop HOB since we don't want to update it
# m.remove_package('HOB')
# # re-write name file before copying
# # causes issue with file path for tab file
# m.write_name_file()

# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')


# %%
dem_data = np.copy(m.dis.top.array)
botm = np.copy(m.dis.botm.array)

# %%
# For the tab files the left column is time (in model units) and the right column is flow (model units)
# Time is days, flow is cubic meters per day
# USGS presents flow in cfs (cubic feet per second)
inflow = pd.read_csv(join(gwfm_dir, 'SFR_data', 'MB_daily_flow_cfs.csv'), index_col = 'datetime', parse_dates = True)
# covnert flow from cubic feet per second to cubic meters per day
inflow['flow_cmd'] = inflow.flow_cfs * (86400/(3.28**3))

# save tabfiles dict to variable for reference
tabfiles_dict = m.sfr.tabfiles_dict

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
all_strt_date = pd.to_datetime(m.dis.start_datetime)
all_dates = all_strt_date + (m.dis.perlen.array.cumsum()-1).astype('timedelta64[D]')
all_end_date = all_dates[-1]
print(all_strt_date, all_end_date)
months = pd.date_range(all_strt_date, all_end_date, freq='MS')
years = pd.date_range(all_strt_date, all_end_date, freq='YS').year.values

# %% [markdown]
#
# ## update scenarios
# Create different recharge scenarios

# %%
div_flow = 300
mar = inflow.copy()
# the max permitted diversion rate is 15.6 cfs
div_rate_cfs = 15.6
print('Total estimated annual diversion %.1f' %(div_rate_cfs*86400/43560*90), 'AF at max diversion')

mar['rch_cfs'] = 0
mar.loc[mar.flow_cfs>div_flow, 'rch_cfs'] = div_rate_cfs
# diversion should only be Dec 15 to Mar 15 or more simply Jan - Mar
mar.loc[~mar.index.month.isin(np.arange(1,4)),'rch_cfs'] = 0
# calculate recharge recharge rate in m3/day
mar['rch'] = mar.rch_cfs*86400*0.3048**3
# need to assign recharge to cells

# %%
# load kautz vineyards to assign to cells
teichert = gpd.read_file(join(gwfm_dir,'Mapping','Kautz_shapefiles', 'Kautz Property.shp'))
rooney = gpd.read_file(join(gwfm_dir,'Mapping','Kautz_shapefiles', 'Kautz Property2.shp'))
vineyards = pd.concat((teichert, rooney)).to_crs(grid_p.crs)
vineyards_grid = gpd.overlay(grid_p, vineyards).drop(columns=['Id','Name','Area'])
vineyards_grid = vineyards_grid[vineyards_grid.geometry.area>(200*200*0.5)]


# %%
# scale the total recharge flux by the area to get the average rate
mar['rch_rate'] = mar.rch/vineyards_grid.geometry.area.sum()
# convert mar to grid level
# need date as index, row, column, rch_rate (m/day)
mar_grid = mar[['rch_rate']].assign(id=0).reset_index().merge(vineyards_grid.assign(id=0)).drop(columns=['id'])
mar_grid_out = mar_grid[['datetime','row','column','rch_rate']].rename(columns={'datetime':'date'})
mar_grid_out.to_csv(join(proj_dir, 'scenarios', 'R1_MAR_max_diversion_for_available_flow.csv'),index=False)

# %%
