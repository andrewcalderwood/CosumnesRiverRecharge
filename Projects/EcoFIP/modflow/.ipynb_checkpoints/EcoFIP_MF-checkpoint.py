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
copy_files=False
# copy_files=True

# %%
concept_name = 'Blodgett_Dam'
concept_name = 'Hanford_Gravel'

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


# load model grid for the test
# gridp = gpd.read_file(join(
zon_stats = gpd.read_file(join(gwfm_dir, 'DIS_data', 'grid_zonal_stats','elevation_m_statistics.shp'))

# %%
# load shapefile that defines geospatial extent for test
# gdf_extent = gpd.read_file(join(concept_dir,'geospatial', 'blodgett_dam_restoration_extent', 'blodgett_dam_restoration_extent.shp'))
# gdf_extent = gdf_extent[gdf_extent.id ==0]
# # get local elevation information for the concept
# proj_zon_stats = gpd.overlay(zon_stats, gdf_extent)



# %%
concept_data = join(main_concept_dir,'data')
os.listdir(concept_data)
# load processed data from cbec
# the units appears to be in feet 
proj_zon_stats = pd.read_csv(join(concept_data,concept_name+'_Elev_Grid_Cells.csv'))
grp_cols = ['node','row','column']
# convert numeric float columns from ft to meters
scale_cols = proj_zon_stats.columns.difference(grp_cols)
proj_zon_stats[scale_cols] *= 0.3048

# %%
wse_grid = pd.read_csv(join(concept_data,concept_name+'_WSE_Grid_Cells.csv'))
# wse is in ft and area is in ft^2
wse_grid.area *= (0.3048**2)
wse_grid.wse *= (0.3048)

# %%
# convert dataframe to geodataframe in case needed
proj_zon_stats = zon_stats[['row','column','geometry']].merge(proj_zon_stats)

# %%
# standard addition of 0-based columns for use in python
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
# field_well = gpd.overlay(fields, gdf_extent)
# use whole cells from zone stats if gdf_extent is unavailable
field_well = gpd.overlay(fields, proj_zon_stats[['mean','geometry']])
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
# Two options:
# 1. Totally replace the streamflow routing with RIV cells with stage based on output from HEC-RAS
# 2. Replace cells at proposed restoration sites with RIV cells
#     - downstream segment needs to divert flow from upstream to skip over the RIV cells
#     - if the RIV is removed in the dry season then the segment data could be redefined to allow flow passing through the segments
#     - **likely need to define strhc1 by segment to be able to assign a value of 0 when RIV cells overlap.** Easiest fix is this, to simply just remove SFR seepage but maintain flow passing through.
# 3. Consider updating strhc1 to use the leakage rates identified with the EcoFIP Tier 2 (soil map + AEM Ksat)

# %%
# just for testing
# m.sfr = flopy.modflow.ModflowSfr2.load(join(model_ws,'MF.sfr'),model=m)

# %%
sfr = m.sfr

# %%
reach_data = pd.DataFrame(sfr.reach_data)

# %%
# proj_zon_stats

# %%
# identify stream reaches impacted by the updated concept
update_reach = proj_zon_stats[['row','column','i','j','mean','min','geometry']].merge(reach_data)
drop_strhc1 = reach_data.reachID.isin(update_reach.reachID)
# make the stream hydraulic conductivity 0 where it will be overlapped
# by floodplain cells driven by the river package
reach_data.loc[drop_strhc1, 'strhc1'] = 0


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
# create short summary of which SFR stream cells we are removing the conductance
import matplotlib.pyplot as plt 
fig,ax=plt.subplots(figsize=(3,2))
update_reach.plot(color='darkblue',ax=ax)
proj_zon_stats.plot(color='none',edgecolor='black',ax=ax)

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
#
# Our river bottom and thickness is most important for defining our hydraulic gradient starting point and it is the elevatoin where we'll see flows drain out.
# - the absolute minimum in a cell may not be well connected but the first quartile likely represents an appropriate starting location. rbot has to be set to the minimum or flows in the stream will constantly have WSE below rbot

# %%
riv_df = proj_zon_stats.copy()

# for the test we can use the average max value
# riv_df['rbot'] = (riv_df['mean']+riv_df['min'])/2

riv_df['rbot'] = riv_df['min']

# get model layer (0-based) for river package
botm_slice = m.dis.botm.array[:, riv_df.i, riv_df.j]
riv_df['k'] = get_layer_from_elev(riv_df['rbot'].values, botm_slice, m.dis.nlay)


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




# %%
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

# %% [markdown]
# Temporary code to substitute in for floodplain stage before data from cbec

# %%

max_stage = riv_df['max'].quantile(0.9)
# based on review of results of riv vs sfr the floodplain
# stage can't be greater than 28 m and floodplain would likely lower this
max_stage = 27
min_stage = riv_df['min'].min()
min_stage, max_stage
# approximate stage for each stress period based on inundation
spd_stage = min_stage+sfr_tab.flow_scale*(max_stage-min_stage)

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
# t0 = time.time()
# riv_dict = dict()
# for t in np.arange(0, m.dis.nper):
#     riv_df_in = riv_df.copy()
#     # for test dataset
#     riv_df_in['stage'] = spd_stage[t]
#     riv_df_in['wetted_area'] = riv_df_in.area*sfr_tab.flow_scale.iloc[t]    
#     # drop river cells if none is saturated
#     riv_df_in = riv_df_in.loc[riv_df_in.stage > riv_df_in.rbot]
#     # calculate the conductance (C = K * l * w/ b)
#     riv_df_in['cond'] = riv_df_in.uhc * riv_df_in.wetted_area / 2
#     # sample relevant columns for input
#     riv_dict[t] = riv_df_in[['k','i','j','stage','cond','rbot']].values

# t1 = time.time()
# print(t1-t0)

# %%
# convert to datetime
wse_grid.date = pd.to_datetime(wse_grid.date)
# subset WSE data to model period
wse_model = wse_grid[(wse_grid.date>=strt_date)&(wse_grid.date<end_date)]
wse_model=wse_model.set_index('date')



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

riv_df_all.to_csv(join(m.model_ws,'river_spd_df.csv'))

# %%
dates = pd.date_range(strt_date, end_date)
# dates

# %%
t=0
d=dates[t]
# riv_df_in.area

# %%
# riv_dict[0][:10]

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
# riv_dict

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

# %% [markdown]
# Update name file to account for new river package
#

# %% [markdown]
# ## Add river obs
# Use river obs to help post-process between recharge in floodplain and in the river.

# %%
# riv_df.shape
# keep_bool

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


nqclfb

# %%
nper = m.dis.nper

# %%
# number of cell groups (number of flow observations), 
# for cosumnes can be 2: stream or floodplain
nqfb = 2
# number of cells total for all cell groups
nqcfb = riv.stress_period_data[0].shape[0]
# number of time obs for all cell groups
nqtfb = nper
# number of times of flow observations for a cell group, list of length nqgb
nqobfb = [nper, nper]
# number of cells in a group, list of length nqgb
nqclfb = []
for n in flag_cells.rvob_grp.unique():
    nqclfb += [int((flag_cells.rvob_grp==n).sum())]
    
# obs name, list of length nqtfb
obsnam = ['klam_gain' + str(s).zfill(3) for s in np.arange(0,nper)]
# obsnam = ['klam_lower', 'klam_upper']
# number of reference stress period, list
irefsp = np.arange(0,nper).tolist()
# time offset from beginning of stress period, list
toffset = [0]*nper
# flow observation for each cell group, list
# negative for flow out of aquifer
flwobs = [-klam_gain_sum]*nper
# flwobs = [klam_gain_lower,klam_gain_upper]
# layer, row, col, factor list of length(nqfb, nqclfb) for each
layer = [chd.stress_period_data[0].k]*nper
row = [chd.stress_period_data[0].i]*nper
column = [chd.stress_period_data[0].j]*nper
fact = np.ones(nqclfb)/nqclfb
factor = [fact]*nper

# flowtype (string) – String that corresponds to the head-dependent flow boundary condition type (CHD, GHB, DRN, RIV)
flowtype = 'RIV'
extension = 'rvob'
chob_unit = [40,42]
chob = flopy.modflow.ModflowFlwob(model = m, nqfb=nqfb,nqcfb=nqcfb,nqtfb=nqtfb,nqobfb=nqobfb,nqclfb=nqclfb,
                          obsnam=obsnam,irefsp=irefsp,toffset=toffset,flwobs=flwobs,
                          layer=layer,row=row,column=column,factor=factor,
                          flowtype = flowtype, extension = extension,
                                 filenames = ['MF.rvob','MF.rvob.out'])
