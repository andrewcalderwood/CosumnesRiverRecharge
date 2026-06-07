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
# Load the existing historical modflow model:
# - the historical modflow model can be created using scripts under andrewcalderwood/CosumnesRiverRecharge/modflow_development/regional
# - alternatively a baseline version of this model will be uploaded in a zip folder on Dropbox to create a static version to use for the scenario analysis. Two version will be used: 1. 2014-2020 for testing purposes and 2. 2000-2022 which will be considered the complete analysis for determining the impacts of floodplain restoration
#
# Adapt only the packages necessary for the Economic anlaysis. 
# - If restoration includes existing farmland then WEL/RCH need to have those cells removed
# - for MAR project inject extra recharge at the site

# %%
# standard python utilities
import os
from os.path import basename, dirname, join, exists, expanduser
import sys
from importlib import reload
import glob
import shutil

import pandas as pd
import numpy as np
import time

# standard python plotting utilities
# import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.dates as mdates

# standard geospatial python utilities
import geopandas as gpd
#import rasterio

# mapping utilities
# import contextily as ctx
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as fm
# from matplotlib.ticker import MaxNLocator


# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')


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
from mf_utility import get_layer_from_elev

# %%
run_dir = 'C://WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'
model_nam = 'input_write_2014_2020'


model_ws = loadpth+model_nam


# %%
load_only = ['DIS','BAS6',
             # 'UPW','OC','SFR','LAK',
            'RCH', 
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
m_top = m.dis.top.array
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
# load datetime reference for the model
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')

# %%
# load model grid for the test
# gridp = gpd.read_file(join(
zon_stats = gpd.read_file(join(gwfm_dir, 'DIS_data', 'grid_zonal_stats','elevation_m_statistics.shp'))
grid_p = zon_stats[['node','row','column','geometry']].copy()
grid_p[['i','j']] = grid_p[['row','column']]-1

# %% [markdown]
# Use Blodgett Dam shapefile as a test case. Also this needs to be generalized if possible to be re-ran as a complete script for each different concept.

# %%
concept_name = 'Blodgett_Dam'
concept_name = 'Hanford_Gravel'
concept_name = 'EG'
# concept_name = 'S60C1_LevSetback_BigWithChans'

concept_ws = join(loadpth,'EcoFIP',concept_name)
m.model_ws = concept_ws
print(concept_name)

# %%
copy_files=False
# copy_files=True

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


# %% [markdown]
# ## Load elevation data

# %%
concept_data = join(main_concept_dir,'data')
os.listdir(concept_data)
# load processed data from cbec
# the units appears to be in feet 
proj_zon_stats = pd.read_csv(join(concept_data,concept_name, concept_name+'_All_Elev_Grid_Cells.csv'))
# the old version for Hanford was only for the concept area
# proj_zon_stats = pd.read_csv(join(concept_data,concept_name, concept_name+'_Elev_Grid_Cells.csv'))
grp_cols = ['node','row','column']
# convert numeric float columns from ft to meters
scale_cols = proj_zon_stats.columns.difference(grp_cols)
proj_zon_stats[scale_cols] *= 0.3048

# %%
# convert dataframe to geodataframe in case needed
proj_zon_stats = zon_stats[['row','column','geometry']].merge(proj_zon_stats)
# standard addition of 0-based columns for use in python
proj_zon_stats[['i','j']] = proj_zon_stats[['row','column']] - 1

# can only use data where there are elevation statistics (limited by DEM/HEC-RAS extent)
proj_zon_stats = proj_zon_stats[~proj_zon_stats[['min','mean']].isna().all(axis=1)]

# %%
# proj_zon_stats.plot('mean', legend=True)
# proj_zon_stats.area

# %% [markdown]
# ## Load WSE data
# The original test version was only for the concept area of Hanford Gravel, the new version will provide the data for the entire domain

# %%
wse_dir = join(concept_data, concept_name,'wide_format')
print(os.listdir(wse_dir))

# %%
# takes about 1.5 minutes to load a csv
area = pd.read_csv(join(wse_dir, 'Cell_Area_TS.csv'))
wse = pd.read_csv(join(wse_dir, 'Cell_WSE_TS.csv'))
# these include data for every single cell in the domain that is wet for at least one timestep

# %%
print('Area dims, full:',area.shape[1], ', any wet:',area.dropna(axis=1, how='all').shape[1])
print('WSE dims, full:',wse.shape[1], ', any wet:', wse.dropna(axis=1, how='all').shape[1])

# %%
# remove cells that are always dry (cbec now only provides wet cells but leave just in case)
area_cln = area.dropna(axis=1, how='all')
wse_cln = wse.dropna(axis=1, how='all')


# %%
# translating the node to row and column

# %%
# convert to long format for easier referencing
area_long = area_cln.melt(id_vars='date')
wse_long = wse_cln.melt(id_vars ='date')

# %%
# drop any individual cell/timestep combos that are dry to reduce data set size
# rename columns to prepare for join
area_long_wet = area_long.dropna(subset='value').rename(columns={'value':'area'})
wse_long_wet = wse_long.dropna(subset='value').rename(columns={'value':'wse'})

# %%
# merge WSE and area datasets to inform RIV package
wse_grid = area_long_wet.merge(wse_long_wet, on= ['date', 'variable'])

# create a column with a clean node number
wse_grid['node'] = wse_grid['variable'].str.extract(r'(\d+)').astype(int)

## new version (EG, LevSetback)
# # wse is in ft and area is in ft^2
wse_grid.area *= (0.3048**2)
wse_grid.wse *= (0.3048)

# prepare the index data to join to the WSE data
row_column= proj_zon_stats[['node','row','column']]
# join the row,column data to the long dataset
wse_grid= wse_grid.merge(row_column, on='node')

# %%
# # old version (Hanford_Gravel)
# wse_grid = pd.read_csv(join(concept_data,concept_name, concept_name+'_WSE_Grid_Cells.csv'))
# # wse is in ft and area is in ft^2
# wse_grid.area *= (0.3048**2)
# wse_grid.wse *= (0.3048)

# %% [markdown]
# Just looking at the Blodgett Dam project footprint makes me think it could make sense to downsize but that is primarily if we wan to see change in groundwater elevations on a sub-cel scale.
# - some of these concepts may take place where there is already no pumping

# %% [markdown]
# # Concept extent
# Geodataframe to define concept extent and use to adjust other BC

# %%
gdf_extent= gpd.read_file(join(gis_dir, 'tier3_alts','cosumnes_T3_siteBoundaries_Final.shp')).to_crs(proj_zon_stats.crs)
# drop all na columns
gdf_extent = gdf_extent.dropna(how='all', axis=1)
if concept_name !='EG':
    #Extract the number to filter the data based on the concept name
    id_str=pd.Series(concept_name).str.extract(r"(S\d+)").values
    id_num= int(id_str[0][0].replace('S', ''))
    gdf_extent = gdf_extent[gdf_extent['id'] == id_num]


# %%
if concept_name != 'EG':
    # select row,columns to use to subset recharge
    proj_grid = proj_zon_stats[['node','row','column','geometry']].copy()
    proj_grid[['i','j']] = proj_grid[['row','column']] - 1
    
    # only adjust recharge where the concepts are
    proj_grid = gpd.overlay(proj_grid, gdf_extent)
    
    # the recharge should be scaled based on how much of a cell will be
    # taken up by the floodplain
    proj_grid['rch_scale'] = 1-(proj_grid.area/(200*200))


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
wel = flopy.modflow.ModflowWel.load(join(model_ws,'MF.wel'),model=m)

# %%
wel = m.wel


# %%
fields = gpd.read_file(join(gis_dir, 'modflow_related','ag_fields_with_well_loc.shp'))


# %%
if concept_name !='EG':
    # identify fields in the AOI that need to have pumping scaled
    field_well = gpd.overlay(fields, gdf_extent)
    # use whole cells from zone stats if gdf_extent is unavailable
    # field_well = gpd.overlay(fields, proj_zon_stats.iloc[0:1][['mean','geometry']])
    # calculate the scaling of pumping based on the field removed
    field_well['area_scale'] = (1-field_well.area)/field_well.full_area
    # get the zero-based row-column
    field_well[['i','j']] = field_well[['row','column']] - 1


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
if concept_name !='EG':
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
    # Create well flopy object
    wel = flopy.modflow.ModflowWel(m, stress_period_data=wel_dict,ipakcb=55)


# %%
# flexible to use existing or updated
wel.write_file()

# %% [markdown]
# ## RCH

# %%
rch = flopy.modflow.ModflowRch.load(join(model_ws,'MF.rch'),model=m)

# %%
if concept_name != 'EG':
    # open recharge array
    rech = rch.rech.array
    rech = rech[:,0,:]
        
    rech_new = np.copy(rech)
    # repeat recharge scaling for each stress period
    rch_scale = np.repeat(np.reshape(proj_grid.rch_scale, (1,-1)),
                          rech_new.shape[0], axis=0)
    # scale the cells where the concept is
    rech_new[:, proj_grid.i, proj_grid.j] *= rch_scale

# %%
if concept_name != 'EG':
    rch_dict = dict()
    
    # for each stress period check if any pumping is in the concept area
    for t in np.arange(0, m.dis.nper):
        rch_dict[t] = rech_new[t]
    
    # Create well flopy object
    rch = flopy.modflow.ModflowRch(m, rech=rch_dict, nrchop=3, ipakcb=55)



# %%
# flexible to use existing or updated
rch.write_file()

# %%
concept_ws

# %% [markdown]
# # Make sure OC saves budget and heads

# %%
# For later model runs when all the data is needed to be saved
month_intervals = (pd.date_range(strt_date,end_date, freq="MS")-strt_date).days
spd = {}
spd = { (j,0): ['save head', 'save budget'] for j in np.arange(0,m.dis.nper,1)}
# spd = { (j,0): ['save head'] for j in np.arange(0,nper,1)}

for j in month_intervals:
    spd[j, 0] = ['save head', 'save budget','print budget']
    # spd[j, 0] = ['save head', 'print budget']
    
oc = flopy.modflow.ModflowOc(model = m, stress_period_data = spd, compact = True)

# %%
oc.write_file()

# %% [markdown]
# # SFR Updates
# Three options:
# 1. Totally replace the streamflow routing with RIV cells with stage based on output from HEC-RAS
# 2. Replace cells at proposed restoration sites with RIV cells
#     - downstream segment needs to divert flow from upstream to skip over the RIV cells
#     - if the RIV is removed in the dry season then the segment data could be redefined to allow flow passing through the segments
# 3. **likely need to define strhc1 by segment to be able to assign a value of 0 when RIV cells overlap.** Easiest fix is this, to simply just remove SFR seepage but maintain flow passing through.
#
# Consider updating strhc1 to use the leakage rates identified with the EcoFIP Tier 2 (soil map + AEM Ksat)

# %%
avg_K_grid = gpd.read_file(join(gis_dir, 'modflow_related','avg_K_grid.shp'))


# %%
# just for testing
# m.sfr = flopy.modflow.ModflowSfr2.load(join(model_ws,'MF.sfr'),model=m)

# %%
sfr = m.sfr

# %%
reach_data = pd.DataFrame(sfr.reach_data)
# merge reaches with new conductivity data
reach_data = reach_data.merge(avg_K_grid[['i','j','K_m_d']], how='left')

# %%
# do a nearest spatial join to look at potential infill from nearby
# create a spatial reach data object
reach_grid = grid_p[['i','j','geometry']].merge(reach_data)
# identify the cells that didn't find conductivity data
reach_grid_fill = reach_grid[reach_grid.K_m_d.isna()].copy().drop(columns=['K_m_d'])
# look to see how far away the nearest cell with data is
# there is still some missing at 200 m (48 reaches) and 400 m (31 reaches)
reach_grid_filled = reach_grid_fill.sjoin_nearest(avg_K_grid[['K_m_d','geometry']], how='inner', max_distance=400, distance_col='dist_m')

# average estimated conductivity properties for each reach
reach_grid_filled = reach_grid_filled.groupby(['iseg','ireach','i','j'])[['K_m_d']].mean(numeric_only=True).reset_index()

# %%
# identify which data is covered by the spatial join
reach_keep = reach_grid.merge(reach_grid_filled[['i','j', 'iseg','ireach']], how='left', indicator=True)
reach_keep = reach_grid[(reach_keep['_merge']=='left_only')]

# append the merged data with the nearest join data
reach_grid_updated = pd.concat((reach_keep, reach_grid.drop(columns=['K_m_d']).merge(reach_grid_filled)))
reach_grid_updated = reach_grid_updated.sort_values(['iseg','ireach']).drop(columns=['geometry'])


# %%
# the majority of the model grid stream cells are covered with AEM data
# and it generally aligns with the existing data
# reach_data.plot(x='reachID',y=['strhc1','K_m_d'])

# reach_grid_updated.plot(x='reachID',y=['strhc1','K_m_d'])

# %%
# replace strhc1, TPROGs ksat with the value from AEM
reach_grid_updated.loc[reach_grid_updated.K_m_d.isna(),'K_m_d'] = reach_grid_updated.loc[reach_grid_updated.K_m_d.isna(),'strhc1']


# %%
# create final reach data for input
reach_data_f = reach_data.drop(columns='K_m_d').merge(reach_grid_updated[['iseg','ireach','K_m_d']])

# %%
# load deep geology data and conductivity to prevent over-write of foothill geology
# need to replace the // with /// to remove URL error
deep_geology = np.loadtxt(join(model_ws,'input_data', 'deep_geology.tsv').replace('//','///'))
m_shape = m.dis.botm.shape
deep_geology = np.reshape(deep_geology, m_shape).astype(bool)

# %%
# vka_deep = m.upw.vka.array
# # filter out to deep geology
# vka_deep[~deep_geology.astype(bool)] = np.nan

# %%
# the existing conductivity in strhc1 is the unsat K from tprogs so for the deep geology this includes
# the land surface formation and the deeper units
reach_data_f['deep_geology'] = deep_geology[reach_data_f.k, reach_data_f.i, reach_data_f.j]
# where there is deep_geology we could overwrite the final vka
# reach_data_f.loc[reach_data.deep_geology,'K_m_d'] = reach_data_f.loc[reach_data.deep_geology,'strhc1']

# %% [markdown]
# One question is whether we want to reset streambed conductances in the upper reaches of Deer Creek and Cosumnes River where we overly what we considered to be lower conductivity foothill formations.

# %%
# identify stream reaches impacted by the updated concept
update_reach = proj_zon_stats[['row','column','i','j','mean','min','geometry']].merge(reach_data_f)
drop_strhc1 = reach_data_f.reachID.isin(update_reach.reachID)
# make the stream hydraulic conductivity 0 where it will be overlapped
# by floodplain cells driven by the river package
reach_data_f.loc[drop_strhc1, 'strhc1'] = 0


# %%
# review sfr cells that are overlapped
reach_chk = update_reach[['row','column','iseg','ireach','reachID','strtop','mean','min']]
reach_chk = reach_chk.sort_values(['iseg','ireach'])
print('Difference in MF stream top and concept cell minimum')
print((reach_chk.strtop-reach_chk['min']).apply([np.mean, np.std]))

# %%
reach_chk.to_csv(join(m.model_ws, 'sfr_riv_df.csv'),index=False)

# %%
# plot to visualize the difference in elevation across space
fig,ax=plt.subplots(figsize=(4,2))
reach_chk.plot(x='reachID',y=['strtop','min'],ax=ax)
plt.ylabel('Eleveation (m)')
plt.legend(['SFR','cbec'])

# %% [markdown]
# It appears that the fine resolution DEM cbec uses identifies significantly deeper stream channel bottoms for each cell or the cell minimum is not representative of the useable channel bottom (e.g., it's a scour hole).

# %%
# # create short summary of which SFR stream cells we are removing the conductance
# import matplotlib.pyplot as plt 
# fig,ax=plt.subplots(figsize=(3,2))
# update_reach.plot(color='darkblue',ax=ax)
# proj_zon_stats.plot(color='none',edgecolor='black',ax=ax)

# %%
sfr.reach_data.strhc1 = reach_data_f.strhc1

# %%
sfr.write_file()

# %% [markdown]
# # Simulate floodplain
# Use the RIV package to simulate floodplain inundation with the output from cbec.
# - originally this was only going to be for the local project area (e.g., new floodplain added by levee adjustment)
# - cbec suggested doing this for the entire domain as they expect to see downstream impacts of concepts so we may want to represent the whole model floodplain to see downstream decreases in recharge
#
# When defining the hydraulic conductivity we will use combined  estimated values from the AEM data and soil data to be consistent with the Tier 2 results. Depending how we look at initial results we could go back to a simpler geology of just 3 layers.

# %% [markdown]
# Our river bottom and thickness is most important for defining our hydraulic gradient starting point and it is the elevatoin where we'll see flows drain out.
# - the absolute minimum in a cell may not be well connected but the first quartile likely represents an appropriate starting location. rbot has to be set to the minimum or flows in the stream will constantly have WSE below rbot

# %%
riv_rch_grid = gpd.read_file(join(gis_dir, 'modflow_related','avg_riv_rch_rate_grid.shp'))


# %%
botm = m.dis.botm.array
ibound = m.bas6.ibound.array

# %%
riv_df = proj_zon_stats.copy()

# for the test we can use the average max value
# riv_df['rbot'] = (riv_df['mean']+riv_df['min'])/2
riv_df['rbot'] = riv_df['min']

# get model layer (0-based) for river package
botm_slice = botm[:, riv_df.i, riv_df.j]
riv_df['k'] = get_layer_from_elev(riv_df['rbot'].values, botm_slice, m.dis.nlay)
# this is producing some RIV cells in active (layers above land surface)

# compare rbot to layer elev
riv_df['m_botm'] = botm[riv_df.k, riv_df.i, riv_df.j]
riv_df['ibound'] = ibound[riv_df.k, riv_df.i, riv_df.j]

# add unsaturated recharge rate to river data
riv_df = riv_df.merge(riv_rch_grid[['i','j','K_m_d']], how='left')
# there are a few cells that have NA values


# %%

# the issue with elevation happens with cells that are just above the model layer bottom
# which means the issue is comparison of the 10 m DEM vs fine resolution
print('There are ', riv_df[riv_df.ibound==0].shape[0],'RIV cells in inactive cells')
# where the layer is to high set lower and check if fixed
riv_df.loc[riv_df.ibound==0, 'k'] += 1
riv_df['ibound'] = ibound[riv_df.k, riv_df.i, riv_df.j]

print('There are ', riv_df[riv_df.ibound==0].shape[0],'RIV cells still in inactive cells')


# %%
# check the typical difference from the minimum to 1st quartile
# 2 meters is not bad since that is similar uncertainty to the
# streambed thickness
(riv_df.Q1st-riv_df['min']).mean()


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
# tprogs realization
r=5
tprogs_info = [80, -80, 320]


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


# %%
# sample unsat K for each floodplain cell
riv_df['gmean_K'] = unsat_K[riv_df.i, riv_df.j] # geology based version
# updated conductivity based on
riv_df['uhc'] = riv_df.K_m_d


# %%
riv_df.plot(x='node', y=['K_m_d','gmean_K'], kind='line', linestyle='', marker='.', label=['EcoFIP','Geology'])
plt.yscale('log')
plt.ylabel('K (m/d)')
# the use of the recharge data instead of geology seems to have slightly increased values in most places with some large decreases
# this makes sense as the soil data will generally be more conductive than the geologic mud

# %% [markdown]
# Calculate the conductance (C = K * l * w/ b)  
# - assume a streambed thickness of 2 m (update with soil map data)
# - cond will vary overtime based on stage and wetted area 
#
# Input for river package: \[lay, row, col, stage, cond, rbot\]
#

# %%
# drop times when there is no WSE reported
#wse_grid = wse_grid.dropna(subset='wse')
# convert to datetime
wse_grid.date = pd.to_datetime(wse_grid.date)
#wse_long_wet.loc[:, 'date'] = pd.to_datetime(wse_long_wet['date'])
print('Number of WSE cells from all time', wse_grid.node.unique().shape[0])
#print('Number of WSE cells from all time', wse_long_wet.variable.unique().shape[0])
# subset WSE data to model period
wse_model = wse_grid[(wse_grid.date>=strt_date)&(wse_grid.date<end_date)].copy()
#wse_model = wse_long_wet[(wse_long_wet.date>=strt_date)&(wse_long_wet.date<end_date)].copy()
wse_model=wse_model.set_index('date')
print('Number of WSE cells from modeled period', wse_model.variable.unique().shape[0])


# %%
keep_riv_cols= ['row','column','i','j','k','rbot','uhc']


# %%
riv_df_all = riv_df[keep_riv_cols].copy()
# identify dates for stage data
riv_df_all = wse_model.reset_index().merge(riv_df_all).set_index('date')

# there is a slight issue in the test dataset
# where WSE is less than the minimum cell level
riv_df_all = riv_df_all[riv_df_all.wse>riv_df_all.rbot]

# the dataset also has dates where the WSE > rbot but the area is null
riv_df_all = riv_df_all.dropna(subset=['wse','area'])


# %%
# no longer should save due to large file size (250 MB)
# won't work for 20 year run (will need to use hdf5 or something
# or else should save external
# still have to save because with post-processing it's useful to have EG and have flexibility
riv_df_all.to_csv(join(m.model_ws,'river_spd_df.csv'))

# %%
# # subset to concept area for saving summarized WSE output
# riv_df_concept = riv_df_all.reset_index().merge(proj_grid[['row','column']])
# riv_df_concept.to_csv(join(m.model_ws,'river_concept_spd_df.csv'))

# %%
# verifies that there are periods when there are more than 53 cells
riv_cell_count = riv_df_all.groupby('date')[['wse']].count()
riv_cell_count.wse.plot()
plt.ylabel('Number of wetted cells')
print('There are',riv_cell_count[riv_cell_count.wse>riv_cell_count.wse.min()].shape[0],'days above the min number of cells.')

# %%
dates = pd.date_range(strt_date, end_date)
# dates

# %%
# runs pretty slow
t0 = time.time()
riv_dict = dict()
for t, d in enumerate(dates):
    # riv_df_in = riv_df[keep_riv_cols].copy()
    # # identify dates for stage data
    # riv_df_in = wse_model.loc[d].merge(riv_df_in)
    if d in riv_df_all.index:
        riv_df_in = riv_df_all.loc[d].copy()

        # # drop river cells if none is saturated
        # calculate the conductance (C = K * l * w/ b)
        riv_df_in['cond'] = riv_df_in.uhc * riv_df_in.area / 2
        # sample relevant columns for input
        riv_dict[t] = riv_df_in[['k','i','j','wse','cond','rbot']].values

t1 = time.time()
print(t1-t0)

# %%
# the WSE doesn't seem to change significantly between any period of the dataframe
# but after plotting time series it shows there are more cells in winter
n=350
print(riv_dict[n].shape[0], riv_df_all.loc[dates[n]].shape[0])
riv_df_all.loc[dates[n]].wse.quantile([0,.5,1])

# %%
riv = flopy.modflow.ModflowRiv(model = m, stress_period_data = riv_dict, ipakcb = 55)

# %%
# river check noted that there were instances with stage below rbot
# solution was to remove stress periods where stage is below rbot
riv.check()

# using updated input dataset from cbec found
# instances of Not a number and stage below rbot

# %%
riv.write_file()

# %%
m.model_ws

# %% [markdown]
# Update name file to account for new river package
#

# %% [markdown]
# ## Add river obs
# Use river obs to help post-process between recharge in floodplain and in the river.  
# Could also use stream stage, bottom, and gw elevation to back-calculate on cell-by-cell basis.  
# - the rvob input crashed modflow so it likely isn't going to be worth including. it will be easier to just use zone budget for periods of interest and perhaps to use the external zone tool rather than the python one or set up the scripts like I did before to do the calculation outside of jupyter

# %%
# differentiate between concept cells in the river vs floodplain
zon_merge = proj_zon_stats[['row','column','i','j','mean','min','geometry']]
flag_cells = zon_merge.merge(reach_data[['i','j','iseg','ireach','reachID']],how='left')


# %%
# flag_cells.reachID.isna().sum()
# floodplain cells are group 1
flag_cells['rvob_grp'] = 1
# river cells are group 2
flag_cells.loc[~flag_cells.reachID.isna(),'rvob_grp'] = 2

# %%
nper = m.dis.nper

# %%
# number of cell groups (number of flow observations), 
# for cosumnes can be 2: stream or floodplain
nqfb = len(flag_cells.rvob_grp.unique())

# number of cells total for all cell groups
# nqcfb = riv.stress_period_data[0].shape[0] # under-estimate depending the stress period
nqcfb = flag_cells.shape[0]
# number of time obs for all cell groups
# does this need to be multiplied by the number of cells or cell groups?
# it needs to be multipled by number of groups because this defines all following inputs
nqtfb = nper*nqfb
# number of times of flow observations for a cell group, list of length nqfb
nqobfb = [nper, nper]
# number of cells in a group, list of length nqgb
nqclfb = []
for n in flag_cells.rvob_grp.unique():
    nqclfb += [int((flag_cells.rvob_grp==n).sum())]


# obs name, list of length nqtfb
fp_nam = ['fldpln' + str(s).zfill(4) for s in np.arange(0,nper)]
strm_nam = ['strm' + str(s).zfill(4) for s in np.arange(0,nper)]
obsnam = fp_nam + strm_nam
# number of reference stress period, list
irefsp = np.arange(0,nper).tolist()*nqfb
# time offset from beginning of stress period, list
toffset = ([0]*nper)*nqfb
# flow observation for each cell group, list
# negative for flow out of aquifer
# use 0 since we don't need this for reference
flwobs = ([0]*nper)*nqfb

# updated above

# layer, row, col, factor list of length(nqfb, nqclfb) for each
# for RIV these change each spd depending on flooding so they should be dynamic
riv_spd = riv.stress_period_data
layer, row, column = [[],[],[]]
for s in np.arange(0, nper):
    layer += [list(riv_spd[s].k)]
    row += [list(riv_spd[s].i)]
    column += [list(riv_spd[s].j)]

# factor is the portion of the simulated gain or loss in the cell that is included in the total gain or loss for this cell group (fn of eq. 5).
# this should be one for all, insert for each rvobs group
factor = []
for n in nqclfb:
    factor += [list(np.ones(n))]



# %%
# flowtype (string) – String that corresponds to the head-dependent flow boundary condition type (CHD, GHB, DRN, RIV)
flowtype = 'RIV'
extension = 'rvob'
chob_unit = [40,42]
chob = flopy.modflow.ModflowFlwob(model = m, nqfb=nqfb,nqcfb=nqcfb,nqtfb=nqtfb,nqobfb=nqobfb,nqclfb=nqclfb,
                          obsnam=obsnam,irefsp=irefsp,toffset=toffset,flwobs=flwobs,
                          layer=layer,row=row,column=column,factor=factor,
                          flowtype = flowtype, extension = extension,
                                 filenames = ['MF.rvob','MF.rvob.out'])

# %%
chob.write_file()

# %%
m.model_ws
