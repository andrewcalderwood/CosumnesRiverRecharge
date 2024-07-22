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
# Load the existing historical modflow model and adapt only the packages necessary for the EcoFIP anlaysis. 
# - If restoration includes existing farmland then WEL/RCH need to have those cells removed
# - Convert SFR to RIV for initial testing at project locations
# - Research further the use of MODFLOW-API to connect to HEC-RAS. This won't be fully done, because EcoFIP doesn't run HEC-RAS in a continuous series but assumes a water surface given a flow time series so it would be that we could try to have an adjustment factor in the model by estimating seepage so that downstream segments have lower flow so lower water surface elevations.

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
# import seaborn as sns
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates

# standard geospatial python utilities
import geopandas as gpd
import rasterio

# mapping utilities
# import contextily as ctx
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as fm
# from matplotlib.ticker import MaxNLocator


# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = join(usr_dir,'Box/research_cosumnes/GWFlowModel')

bas_dir = join(gwfm_dir, 'BAS6')
# proj_dir = join(gwfm_dir,'Projects','EcoFIP')
# plt_dir = join(proj_dir,'figures/')

# %%
lwa_dir = join(usr_dir, 'LWA Dropbox','01_Project-Teams')
proj_dir = join(lwa_dir, '669.03 - DWR Cosumnes Floodplain Recharge')
main_concept_dir = join(proj_dir, 'Concepts')
gis_dir = join(main_concept_dir,'GIS')
concept_dir = join(main_concept_dir,'Blodgett_Dam')


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
run_dir = 'C://WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'


model_ws = loadpth+model_nam


# %%
copy_files=False

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

# %%
# load datetime reference for the model
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')

# %% [markdown]
# Use Blodgett Dam shapefile as a test case. Also this needs to be generalized if possible to be re-ran as a complete script for each different concept.

# %%
concept_name = 'Blodgett_Dam'
concept_ws = join(loadpth,'EcoFIP',concept_name)
m.model_ws = concept_ws

# %%
import shutil

# %%
if copy_files:
    # copy over modflow files that aren' updated
    # glob.glob(join(model_ws, 'MF.*'))
    files = pd.DataFrame(os.listdir(model_ws),columns=['fn'])
    # only keep modflow related files
    files = files[files.fn.str.contains('MF')]
    # drop files to be updated
    update = 'rch|wel|sfr'
    files = files[~files.fn.str.contains(update)]
    # copy files over to new workspace
    for f in files.fn:
        shutil.copy(join(model_ws,f), join(concept_ws, f))


# %%
concept_dir = join(gwfm_dir,'Projects',concept_name)
# load shapefile that defines geospatial extent for test
gdf_extent = gpd.read_file(join(concept_dir,'geospatial', 'blodgett_dam_restoration_extent', 'blodgett_dam_restoration_extent.shp'))
# gdf_extent = gdf_extent[gdf_extent.id ==0]
# load model grid for the test
# gridp = gpd.read_file(join(
zon_stats = gpd.read_file(join(gwfm_dir, 'DIS_data', 'grid_zonal_stats','elevation_m_statistics.shp'))

# %%
# get local elevation information for the concept
proj_zon_stats = gpd.overlay(zon_stats, gdf_extent)
proj_zon_stats[['i','j']] = proj_zon_stats[['row','column']] - 1



# %%
# proj_zon_stats.plot('mean', legend=True)
# proj_zon_stats.area

# %% [markdown]
# Just looking at the Blodgett Dam project footprint makes me think it could make sense to downsize but that is primarily if we wan to see change in groundwater elevations on a sub-cel scale.
# - some of these concepts may take place where there is already no pumping

# %% [markdown]
# # WEL/RCH Updates
# Remove or scale WEL/RCH where new floodplain interaction is:
# - ideally we scale the the pumping based on the area of a cell converted to floodplain
# - recharge is easier since it is already distributed cell-by-cell.
#
#

# %% [markdown]
# ## WEL

# %%
# reload old well package for testing
# wel = flopy.modflow.ModflowWel.load(join(model_ws,'MF.wel'),model=m)

# %%
wel = m.wel


# %%
fields = gpd.read_file(join(gis_dir, 'modflow_related','ag_fields_with_well_loc.shp'))


# %%
# identify fields in the AOI that need to have pumping scaled
field_well = gpd.overlay(fields, gdf_extent)
# calculate the scaling of pumping based on the field removed
field_well['area_scale'] = (1-field_well.area)/field_well.full_area
# get the zero-based row-column
field_well[['i','j']] = field_well[['row','column']] - 1
# field_well

# %%
def update_well(wel_arr, field_well):
    """
    Given an array of well package input, scale the pumping if
    the spatial extent of a field is reduced
    wel_arr: modflow well array with k,i,j, flux
    field_well: dataframe with i,j and area_scale to reduce pumping
    """
    wel_cols = ['k','i','j','flux']
    # turn array into dataframe for easier referencing
    wel_df = pd.DataFrame(wel_arr)
    # identify the wells with pumping to update based on scaled area
    wel_update = wel_df.merge(field_well, how='inner')
    wel_update.flux *= wel_update.area_scale
    # identify well data to keep as is
    wel_keep = wel_df.merge(field_well, how='outer', indicator=True)
    wel_keep = wel_keep[wel_keep._merge=='left_only']
    # concatenate the datasets to write back into the wel dictionary
    wel_new = pd.concat((wel_keep, wel_update))[wel_cols]
    return(wel_new)


# %%
# this script takes close to a minute to run
time0 = time.time()
# it is likely in-efficient with the two merge function
wel_dict = dict()

# for each stress period check if any pumping is in the concept area
for t in np.arange(0, m.dis.nper):
    wel_arr = wel.stress_period_data[t]
    wel_new = update_well(wel_arr, field_well)
    wel_dict[t] = wel_new.values

time1 = time.time()
print(time1-time0)

# %%
# Create well flopy object
wel = flopy.modflow.ModflowWel(m, stress_period_data=wel_dict,ipakcb=55)

# %%
wel.write_file()

# %% [markdown]
# ## RCH

# %%
rch = m.rch
# open recharge array
rech = rch.rech.array
rech = rech[:,0,:]

# %%
# select row,columns to use to subset recharge
proj_grid = proj_zon_stats[['node','row','column','geometry']].copy()
proj_grid[['i','j']] = proj_grid[['row','column']] - 1
# the recharge should be scaled based on how much of a cell will be
# taken up by the floodplain
proj_grid['rch_scale'] = 1-(proj_grid.area/(200*200))
# proj_grid.rch_scale

# %%
rech_new = np.copy(rech)
# repeat recharge scaling for each stress period
rch_scale = np.repeat(np.reshape(proj_grid.rch_scale, (1,-1)),
                      rech_new.shape[0], axis=0)
# scale the cells where the concept is
rech_new[:, proj_grid.i, proj_grid.j] *= rch_scale

# %%

rch_dict = dict()

# for each stress period check if any pumping is in the concept area
for t in np.arange(0, m.dis.nper):
    rch_dict[t] = rech_new[t]


# %%
# Create well flopy object
rch = flopy.modflow.ModflowRch(m, rech=rch_dict, nrchop=3, ipakcb=55)

# %%
rch.write_file()


# %% [markdown]
# # SFR Updates
# Two options:
# 1. Totally replace the streamflow routing with RIV cells with stage based on output from HEC-RAS
# 2. Replace cells at proposed restoration sites with RIV cells
#     - downstream segment needs to divert flow from upstream to skip over the RIV cells
#     - if the RIV is removed in the dry season then the segment data could be redefined to allow flow passing through the segments
#     - **likely need to define strhc1 by segment to be able to assign a value of 0 when RIV cells overlap.** Easiest fix is this, to simply just remove SFR seepage but maintain flow passing through.

# %%
# just for testing
# m.sfr = flopy.modflow.ModflowSfr2.load(join(model_ws,'MF.sfr'),model=m)

# %%
sfr = m.sfr

# %%
reach_data = pd.DataFrame(sfr.reach_data)

# %%
# identify stream reaches impacted by the updated concept
update_reach = reach_data.merge(proj_grid[['row','column','i','j']])
drop_strhc1 = reach_data.reachID.isin(update_reach.reachID)
# make the stream hydraulic conductivity 0 where it will be overlapped
# by floodplain cells driven by the river package
reach_data.loc[drop_strhc1, 'strhc1'] = 0


# %%
sfr.reach_data.strhc1 = reach_data.strhc1

# %%
sfr.write_file()

# %% [markdown]
# # Simulate floodplain
# Use the RIV package to simulate floodplain inundation with the output from cbec.
#
# When defining the hydraulic conductivity, it may make sense to use the estimated values from the AEM data to improve consistency. Then we could also perhaps go back to a simpler geology of just 3 layers.

# %%
from mf_utility import get_layer_from_elev

# %% [markdown]
# It may make sense to move the hydraulic conductivity work to pre-processing or have the standard model save as needed for each tprogs realization

# %%
riv_df = proj_zon_stats.copy()

# for the test we can use the average max value
riv_df['rbot'] = (riv_df['mean']+riv_df['min'])/2


# get model layer (0-based) for river package
botm_slice = m.dis.botm.array[:, riv_df.i, riv_df.j]
riv_df['k'] = get_layer_from_elev(riv_df['rbot'].values, botm_slice, m.dis.nlay)


# %%
riv_df[['row','column','i','j','k','rbot']].to_csv(join(m.model_ws, 'river_ijk.csv'))

# %%
from mf_utility import param_load

# %%
gel_dir = join(gwfm_dir,'UPW_data')
# # load data from Steven
params = param_load(concept_ws, gel_dir, 'ZonePropertiesInitial_Maples.csv')  
params = params.set_index('Zone')
# temporarily increase Mud from 0.0017 to 0.01 to ease convergence
params.loc[4,'K_m_d'] = 0.01

# %%
tprogs_fxn_dir = doc_dir+'/GitHub/CosumnesRiverRecharge/tprogs_utilities'
if tprogs_fxn_dir not in sys.path:
    sys.path.append(tprogs_fxn_dir)
# sys.path
import tprogs_cleaning as tc

reload(tc)

# %%
import h5py
tprogs_name = 'tprogs_final'
tprogs_fn = join(gel_dir, tprogs_name+'.hdf5')
r=5
tprogs_info = [80, -80, 320]

m_top = m.dis.top.array
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')

# tprogs_line = np.loadtxt(tprogs_files[r])
# # # # using m.dis.top to account for cases when top is above dem
# masked_tprogs= tc.tprogs_cut_elev(tprogs_line, m.dis.top.array, tprogs_info)

with h5py.File(tprogs_fn, mode='r') as f:
    tprogs_arr = f['tprogs']['r'+str(r).zfill(3)][:].astype(float)
masked_tprogs= tc.tprogs_arr_cut_elev(tprogs_arr, m_top, tprogs_info)
K, Sy, Ss = tc.int_to_param(masked_tprogs, params, porosity=False)


# %%
from scipy.stats import gmean
riv_dem = np.copy(dem_data)
riv_dem[riv_df.i, riv_df.j] = riv_df['rbot']
# use 10 mean below land surface or historic wse whichever is shallower
riv_vka_bot = riv_dem-10

# sample unsaturated zone conductivity for uhc
unsat_K_all  = tc.get_tprogs_for_elev(K, riv_dem, riv_vka_bot, tprogs_info)
# calculate geometric mean for the unsat zone routing
unsat_K = gmean(unsat_K_all, axis=0)

# sample unsat K for each floodplain cell
riv_df['uhc'] = unsat_K[riv_df.i, riv_df.j]


# %% [markdown]
# Calculate the conductance (C = K * l * w/ b)  
# - assume a streambed thickness of 2 m
# - cond will vary overtime based on stage and wetted area 
#
# Input for river package: \[lay, row, col, stage, cond, rbot\]
#

# %%
# load sfr tab file for streamflow relation if needed
sfr_tab = pd.read_csv(join(m.model_ws, 'MF.tab'), delimiter='\t', header=None)
sfr_tab.columns=['time','flow']
sfr_tab['log10_flow']=np.log10(np.where(sfr_tab.flow==0,1,sfr_tab.flow))
sfr_tab['flow_scale'] = sfr_tab.log10_flow/sfr_tab.log10_flow.max()

# %%
max_stage = riv_df['max'].quantile(0.9)
# based on review of results of riv vs sfr the floodplain
# stage can't be greater than 28 m and floodplain would likely lower this
max_stage = 27
min_stage = riv_df['min'].min()
min_stage, max_stage
# approximate stage for each stress period based on inundation
spd_stage = min_stage+sfr_tab.flow_scale*(max_stage-min_stage)


# %%
spd_stage.to_csv(join(m.model_ws,'river_stage.csv'))

# %%
# test version
riv_df['stage'] = riv_df['max'].mean()
# approximate wetted area (floodplain length x width
# should be provided by cbec as function of WSE
riv_df['wetted_area'] = riv_df.area


# %%
m.dis.nper, len(sfr_tab)

# %%
t0 = time.time()
riv_dict = dict()
for t in np.arange(0, m.dis.nper):
    riv_df_in = riv_df.copy()

    riv_df_in['stage'] = spd_stage[t]
    riv_df_in['wetted_area'] = riv_df_in.area*sfr_tab.flow_scale.iloc[t]
    # drop river cells if none is saturated
    riv_df_in = riv_df_in.loc[riv_df_in.stage > riv_df_in.rbot]
    # calculate the conductance (C = K * l * w/ b)
    riv_df_in['cond'] = riv_df_in.uhc * riv_df_in.wetted_area / 2
    
    # sample relevant columns for input
    riv_dict[t] = riv_df_in[['k','i','j','stage','cond','rbot']].values


t1 = time.time()
print(t1-t0)

# %%
# riv_dict

# %%
riv = flopy.modflow.ModflowRiv(model = m, stress_period_data = riv_dict, ipakcb = 55)

# %%
# river check noted that there were instances with stage below rbot
# solution was to remove stress periods where stage is below rbot
riv.check()

# %%
riv.write_file()

# %% [markdown]
# Update name file to account for new river package
