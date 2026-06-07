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
# Goal: compare results of original SWB (daily irrigation) to results of optimized SWB with spaced out irrigation

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

loadpth = run_dir +'/Cosumnes/Regional/'

model_nam = 'input_write_2014_2022'

base_model_ws = loadpth+model_nam

# %%
# runs much slower with rch package but better to be sure on the input
load_only = ['DIS','BAS6','SFR', 'RCH']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only,
                              )

nrow,ncol,nlay,delr,delc = (m.dis.nrow, m.dis.ncol, m.dis.nlay, m.dis.delr, m.dis.delc)



# %%
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')


# %%
end_date

# %%
dt_ref['totim'] = (dt_ref.dt-dt_ref.dt.iloc[0]).dt.days+1


# %% [markdown]
# To compare recharge/percolation it might be easier to just grab the MODFLOW input to evaluate spatially where issues might be.

# %%
uzf_dir = join(gwfm_dir, 'UZF_data')

# %%
# function can't handle data less than a year
# from swb_utility import load_swb_data, yr_lim_from_dates
# data is already in array format
# perc = load_swb_data(strt_date, pd.to_datetime('2022-9-30'), 'field_percolation', uzf_dir)
# the output from UZF/Basic soil budget - fields is the fields but then is converted to array format



# %%
# percolation is scaled as a fraction of cells before conversion to arrays so units are correct m/day
# perc_sum = perc.sum(axis=(1,2))

# %%
# gde_dir = join(uzf_dir,'shp_GDE_TFT')

# # locations of GDes and native vegetation
# # the pre-calculated rooting depth assumes deeper roots for GDEs than regular natives
# GDE_all = gpd.read_file(join(gde_dir,'GDE_and_native_cell.shp'))

# # convert rooting depth to array format for modflow input, hydrographs in wells show drawdown to about 10 m
# # below ground so should use 10 m for all gde
# ext_dp = np.full((nrow,ncol),2)
# ext_dp[(GDE_all.row-1).astype(int), (GDE_all.column-1).astype(int)] = GDE_all.rtg_dp


# %%
# # where GDEs are active the rain should be directly applied instead of SWB percolation
# # as we don't want to double count ET
# adj_perc = perc.copy()
# # initially rain was applied only where there was a deep rooting depth
# et_rain_bool = (evt_active)&(ext_dp>2)
# adj_perc[:, et_rain_bool] = rain_arr[:, et_rain_bool]
# # model results show under prediction of levels in the floodplain
# # so it may make sense to use rain everywhere there is ET
# # only do it for days with floodplain inundation
# adj_perc[evt_active_all] = rain_arr[evt_active_all]
# adj_perc_fp = adj_perc.copy()

# %%
# also loading applied water if interested in comparison
# AW = load_swb_data(strt_date, end_date, 'field_applied_water', uzf_dir)


# %%
# grab recharge array from flopy object
rech_arr = m.rch.rech.array[:,0]

# %%
rech_arr = rech_arr[:,0]

# %%

# %% [markdown]
# # Load data from baseline economic run
#
# Goal: 
# - check that percolation is not too different in irrigation vs wet season
# - see if lack of soil moisture passing is driving this issue. a quick check would be running the reference swb ag winter with field capacity instead of wilting point to see the fractional change
# - should also look at spatial extent of recharge to see if an area is missing, can just use modflow for this

# %%
loadpth = run_dir +'/Cosumnes/Economic/'
model_nam = 'input_write_2014_2022'
scen = 'R200'
# it is probably better to create a slightly different file name then to copy these over for a set scenario
econ_model_ws = join(loadpth, model_nam+'_'+scen)

all_run_dates = pd.read_csv(join(econ_model_ws, 'crop_modflow', 'all_run_dates.csv'))

# I should update this file path
swb_ws = join(econ_model_ws, 'rep_crop_soilbudget')


# %%
# to be sure we evaluate the actual input, combine recharge across all input models then evaluate spatial variability
rech_arr_per_all = np.zeros(rech_arr.shape)
npos = 0
for n, d in enumerate(all_run_dates.date[:-1]):
# runs much slower with rch package but better to be sure on the input
    load_only = ['DIS', 'RCH']
    m_per = flopy.modflow.Modflow.load('MF.nam', model_ws= join(econ_model_ws, 'crop_modflow',d), 
                                    exe_name='mf-owhm', version='mfnwt',
                                  load_only = load_only,
                                  )
    rech_arr_per = m_per.rch.rech.array[:,0]
    rech_arr_per_all[npos:npos+len(rech_arr_per)] = rech_arr_per
    npos += len(rech_arr_per)

# %%
npos, npos+len(rech_arr_per)

# %%
rech_per_sum = rech_arr_per_all.sum(axis=0)
rech_sum = rech_arr.sum(axis=0)
vmax = np.max((rech_sum.max(),rech_per_sum.max()))

fig,ax = plt.subplots(2,1)
ax[0].imshow(rech_per_sum, vmax=vmax)
# plt.colorbar(shrink=0.6)

im = ax[1].imshow(rech_sum, vmax=vmax)

# Adjust layout to make room for colorbar
fig.subplots_adjust(right=0.85)
# Define [left, bottom, width, height] for the colorbar
cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])

fig.colorbar(im, shrink=0.6, cax= cbar_ax)

# %% [markdown]
# Recharge map results shows that the spatial pattern is similar but the overall recharge rates are considerably lower. It also looks like the "flooding" recharge in the lower floodplain isn't passed through all the way. May need to add this back? Or check native vs ag lands distribution.

# %%

soil_path = join(uzf_dir,'clean_soil_data')
# reference for field to grid cell
field_ids = 'parcel'
grid_soil_ag = pd.read_csv(join(soil_path, field_ids+'_field_to_cell.csv'), index_col=0)

field_ids = 'native'
# reference for field to grid cell
grid_soil = pd.read_csv(join(soil_path, field_ids+'_field_to_cell.csv'), index_col=0)

# %%
# columns are the ntaive land use fields
pc_df = pd.read_csv( join(proj_dir,'native_parcel_zonalstats', 
                   'swb_percolation.csv'),index_col=0)
# long format to join with row,column
pc_df_long = pc_df.melt(ignore_index=False, 
                        var_name='UniqueID',
                       value_name='pc')
pc_df_long.UniqueID=pc_df_long.UniqueID.astype(int)

# %%
# join field to row,column
native_pc = pc_df_long.reset_index().merge(grid_soil)
native_pc = native_pc.rename(columns={'index':'date'})
# calculate recharge rate by scaling for fraction of cell covered 
native_pc['rch_rate'] = native_pc.pc*native_pc.field_area/(200*200)
# update index to datetime
native_pc = native_pc.set_index('date')
native_pc.index=pd.to_datetime(native_pc.index)

# %%
pc_per_all = pd.DataFrame()
pc_non_irr_all = pd.DataFrame()
for n, d in enumerate(all_run_dates.date[1:-1]):
    year = pd.to_datetime(d).year
    # create locally referenced dt_ref file to avoid reloading dis
    dt_ref_per = dt_ref[(dt_ref.dt<=all_run_dates.date[n+1])&(dt_ref.dt>=d)].copy()
    dt_ref_per.kstpkper = list(zip(np.zeros(len(dt_ref_per),dtype=int), np.arange(len(dt_ref_per))))

    econ_ws_yr =join(econ_model_ws, d)
    # wb_per, out_cols_per, in_cols_per  = clean_wb(econ_ws_yr, dt_ref_per)
    pc_df_all = pd.read_csv(join(swb_ws, 'output', 'pc_all'+str(year)+'.csv'),index_col=0)
    pc_per_all = pd.concat((pc_per_all, pc_df_all))
    # data from non-irrigation period (rain applied to fields)
    pc_df = pd.read_csv(join(econ_model_ws,
                  'ag_non_irr_swb_percolation_'+str(year)+'.csv'),index_col=0)
    pc_non_irr_all = pd.concat((pc_non_irr_all, pc_df))

# econ_ws_yr

# %%
# need to look at combined portion of these and compare to original estimates
pc_per_all
pc_non_irr_all
pc_df_long
