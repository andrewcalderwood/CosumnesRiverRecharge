# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Script to run the code for the crop choice model, irrigation optimizer and modflow in a cycle for the years of interest.
#
# 1. Run the crop choice model
# 2. Run the soil water budget optimization
# 3. Update WEL/RCH packages then run MF
# 4. Re-run the soil water budget without optimization applying the previous irrigation rates to get the actual profit and to account for crop change
# 5. Start the next year
#
# And run different management and/or recharge scenarios. This is done via the batch file and specifying which model workspace is used.

# %%
import sys
import os
from os.path import join,exists, dirname, basename, expanduser
import time 
import re

import h5py
import numpy as np
import pandas as pd
import geopandas as gpd

from importlib import reload

# %%
# for batch file case we want to ignore the consant h5py deprecation warning
import warnings

warnings.filterwarnings("ignore")
# warnings.filterwarnings("ignore", category=DeprecationWarning) 


# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')


# %%
# updated version specifies concept_name and copy_files here so it can be easily
# seen as these are the main update to make in a script
in_data = sys.argv

if 'ipykernel' in in_data[0]:
    # m_nam = 'input_write_2014_2020_R3'
    # m_nam = 'input_write_2014_2020'
    m_nam = 'input_write_2014_2025_R300'

    # input_name = 'static_model_inputs.xlsx'
    input_name = 'static_model_inputs_no_p_o.xlsx'
    year_load_var_in = "False"
    year_load_var_in = 2019
    # identify which method for using/creating the optimized SWB
    # new_local, existing_local, existing_shared
    # and another option existing_local_and_crop to not write predict crop
    create_rep_swb = "new_local"

else:
    m_nam = in_data[1]
    # added option to specify different static_model_inputs to allow easier adjustment of operating costs, revenue
    input_name = in_data[2]
    year_load_var_in = in_data[3]
    if year_load_var_in != "False":
        year_load_var_in = int(year_load_var_in)
    create_rep_swb =  in_data[4]



print('sys.argv[1] (m_nam) is...')
print(m_nam)

print('sys.argv[2] (input_name) is...')
print(input_name)

print('sys.argv[3] (year_load_var_in) is...')
print(year_load_var_in)

print('sys.argv[4] (create_rep_swb) is...')
print(create_rep_swb)

t_start = time.time()


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)


# %%
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# %%
git_dir = join(doc_dir,'GitHub','CosumnesRiverRecharge')
add_path(join(git_dir,'python_utilities'))

# %%
from report_cln import base_round
from mf_utility import get_layer_from_elev
import flopy_utilities
reload(flopy_utilities)
from flopy_utilities import update_bas_options

# %%
import parcelchoicemodelupdate.f_predict_landuse
reload(parcelchoicemodelupdate.f_predict_landuse)
from parcelchoicemodelupdate.f_predict_landuse import predict_crops, predict_crops_rand

# %%
import functions.Basic_soil_budget_monthly as swb
# reload(swb)

# import f_rep_swb_profit_opt
# reload(f_rep_swb_profit_opt)
from f_rep_swb_profit_opt import load_run_swb

# %%
# import functions.f_gw_dtw_extract
# reload(functions.f_gw_dtw_extract)
from functions.f_gw_dtw_extract import sample_dtw, avg_heads
# from functions.f_gw_dtw_extract import get_dtw
from functions.f_gw_dtw_extract import calc_simple_dtw

# %%
# from functions import data_functions
# reload(functions.data_functions)
from functions.data_functions import read_crop_arr_h5

from functions.data_functions import init_h5


# %%
# import functions.output_processing
# reload(functions.output_processing)
from functions.output_processing import get_local_data, out_arr_to_long_df
from functions.output_processing import get_wb_by_parcel

# %%
import reference_swb_ag_winter 
reload( reference_swb_ag_winter )
from reference_swb_ag_winter import run_swb_ag_winter


# %%
# option 2 is to call a function
import f_summarize_output
reload(f_summarize_output)
from f_summarize_output import summarize_output_year

# %% [markdown]
# If we are going to run the crop/swb for a certain period then there should already have been a MODFLOW model run for that same period with estimates of streamflow and precipitation to drive the other boundary conditions.

# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
# loda model scenario reference sheet
scenario_summary = pd.read_excel(join(data_dir, 'scenario_summary.xlsx'))
# select the current model run scenario based on model and scenario name from batch file
scenario_info = scenario_summary[scenario_summary.m_nam==m_nam] # & (scenario_summary.m_nam==m_nam)
# create series to make easier data sampling
scenario_info = scenario_info.iloc[0]
# create into variables for easier referencing
sw_con = (scenario_info.sw_con).astype(float)
gw_con = (scenario_info.gw_con).astype(float)
scenario_name = scenario_info.scenario_name
base_m_nam = scenario_info.base_m_nam

# print out additional summary info
print('SW Constraint %.2f' %sw_con, 'inches')
print('GW Constraint %.2f' %gw_con, 'inches')

print('base_nam is...')
print(base_m_nam)
print('scenario_name is...')
print(scenario_name)
print('\n\n')
# scenario_info


# %%
# load parcel data for reference as needed
parcels = gpd.read_file(join(proj_dir,'Parcels shapefile/parcels.shp'))
parcels['area_m2'] = parcels.geometry.area

# load dataframe that assigns each field to cells for recharge
field_id = 'parcel'
field_df = pd.read_csv(join(uzf_dir, 'clean_soil_data', field_id+'_field_to_cell.csv'),index_col=0)

# load dataframe with location (row,column) of wells
well_loc_df = gpd.read_file(join(gwfm_dir, 'WEL_data','parcels_to_wells','parcels_to_wells.shp'))
well_loc_merge = well_loc_df[['UniqueID','row','column','depth_m']].copy()
well_loc_merge = well_loc_merge.merge(parcels[['UniqueID', 'area_m2']])
well_loc_merge.UniqueID = well_loc_merge.UniqueID.astype(int)

# %%
# should we may the loadpth economic instead of Regional?
# should make drive extension more flexible/in batch file
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'F:/WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'D:/WRDAPP/GWFlowModel/Cosumnes/Economic'

# update to different modflow models here, next step is using the 20 year model
# m_nam = 'historical_simple_geology_reconnection'
# m_nam = 'input_write_2014_2020'

# m_nam = 'input_write_2000_2022'
# define modflow model WS to reference for modflow input
m_model_ws = join(dirname(loadpth), 'Regional', base_m_nam)
m_model_ws

# %%
load_only=['DIS', 'BAS6']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= m_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only)

# %%
# save ibound array for when updating bas6 for each run
ibound = m.bas6.ibound.array
# bottom array is needed for referencing well layer
botm = m.dis.botm.array

# %%
m_dim = (m.dis.nlay, m.dis.nrow,  m.dis.ncol)
np.savetxt(join(m_model_ws, 'model_metadata.txt'), m_dim)


# %% [markdown]
# ## General model input
# Data that is not impacted by crop choice or irrigation decisions
# - Well data for pumping
# - Native land use recharge
# - recharge scenarios (requires use of a different model_ws)

# %%
# may want to add a year start/end or put in subfolder of modflow model
# so we can run different versions
model_ws = join(loadpth, m_nam)

# %%
## add domestic pumping
wel_dir = join(gwfm_dir, 'WEL_data')
# load prepared daily domestic use data
dom_use = pd.read_csv(join(wel_dir, 'domestic_water_use.csv'), index_col=0, parse_dates=True)
# load data of locations of domestic wells
dom_loc = pd.read_csv(join(wel_dir, 'ag_res_parcel_domestic_wells.csv'), index_col=0)
# make row,column 0 based
dom_loc[['row','column']] = (dom_loc[['row','column']]-1).astype(int)
# aggregate to the cell level, summing area will keep water usage scaling correct
dom_loc = dom_loc.groupby(['node','row','column', 'CITY']).sum(numeric_only=True).reset_index()
# get domestic well layers
dom_wel_bot = (m.dis.top[dom_loc.row, dom_loc.column]- dom_loc.fill_depth_m).values
dom_loc['layer'] = get_layer_from_elev(dom_wel_bot, botm[:,dom_loc.row, dom_loc.column], m.dis.nlay)

# use either the total area or expected fraction of irrigated area
dom_loc['pump_scale'] = dom_loc.used_area_acres # avg of 25% less than area_acres (sometimes more though)
# alternate to avoid over-estimation just stick with 2AF/year
dom_loc['pump_scale'] = 1 # avg of 25% less than area_acres (sometimes more though)


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
soil_path = join(uzf_dir,'clean_soil_data')
field_ids = 'native'
# reference for field to grid cell
grid_soil = pd.read_csv(join(soil_path, field_ids+'_field_to_cell.csv'), index_col=0)


# %%
field_ids = 'parcel'
grid_soil_ag = pd.read_csv(join(soil_path, field_ids+'_field_to_cell.csv'), index_col=0)

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
# drop dates and cells with zero recharge as the default recharge is 0
native_pc = native_pc[native_pc.rch_rate!=0]

# %%
# quick budget checking
# native_pc['year'] = native_pc.index.year
# yr_sum = native_pc.groupby(['year','row','column'])['rch_rate'].sum()
# yr_sum = yr_sum.reset_index()
# yr_sum.groupby('year').mean()*3.28*12 # m/yr, 0.67 inches/yr for max year
# we are underestimating recharge

# %%
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')
# well_loc_merge.depth_m
## sample well layer
# get top elevations to calculate well screen elevation
# well_loc_merge['top_elev'] = m.dis.top.array[well_loc_merge.row-1, well_loc_merge.column-1]
# well_loc_merge['well_elev'] = well_loc_merge.top_elev - well_loc_merge.depth_m
# switched to DEM since model top can sometimes be different than reality
well_loc_merge['dem'] = dem_data[well_loc_merge.row-1, well_loc_merge.column-1]
well_loc_merge['well_elev'] = well_loc_merge.dem - well_loc_merge.depth_m*0.9
# make layer 1 based
well_loc_merge['layer'] = get_layer_from_elev(well_loc_merge.well_elev.values, 
                    botm[:, well_loc_merge.row-1, well_loc_merge.column-1], m.dis.nlay)+1
# make 0-based columns
well_loc_merge[['k','i','j']] = well_loc_merge[['layer','row','column']]-1
# save file with well layers so it can be easily referenced by other scripts (e.g., summarize_output)
parcel_wells = well_loc_merge[['UniqueID','dem','layer', 'row','column']].copy()
parcel_wells.to_csv(join(model_ws, 'crop_modflow', 'parcel_wells_with_layer.csv'), index=False)


# %% [markdown]
# Specify dates for the entire model period and seasonal dates to start/stop the model.

# %%
all_strt_date = pd.to_datetime(m.dis.start_datetime)
all_dates = all_strt_date + (m.dis.perlen.array.cumsum()-1).astype('timedelta64[D]')
all_end_date = all_dates[-1]
print(all_strt_date, all_end_date)
months = pd.date_range(all_strt_date, all_end_date, freq='MS')
years = pd.date_range(all_strt_date, all_end_date, freq='YS').year.values

# %% [markdown]
# ## Scenarios
# Need to start thinking of the best way to save different model scenarios with the high level model_ws change. Even with different static_model_input.xlsx as long as scenario_name is kept consistent then it can reference back to this.

# %%

if scenario_name=='none':
    # # baseline scenario is no extra recharge
    mar_grid = pd.DataFrame()
else:
    # # this will need to be coded to be flexible and to change the model_ws
    mar_grid = pd.read_csv(join(proj_dir, 'scenarios', scenario_name+'.csv'), index_col='date', parse_dates=['date'])
    print('Using mar_grid for', scenario_name)



# %%

# load summary excel sheet on irrigation optimization
# this will specify the date ranges to run and pause
# fn = join(data_dir,'static_model_inputs.xlsx') # original
fn = join(data_dir ,input_name) # new version allows alternates


season = pd.read_excel(fn, sheet_name='Seasons', comment='#')
irr_est = pd.read_excel(fn, sheet_name='Irr_establish', comment='#')

os.makedirs(join(model_ws, 'output_clean'), exist_ok=True)
os.makedirs(join(model_ws,'crop_soilbudget','field_dtw'),exist_ok=True)

# save log by date so we can see old versions
os.makedirs(join(model_ws, 'log'), exist_ok=True)
if 'ipykernel' in in_data[0]:
    print('print to jupyter')
else:
    print('printing remaining output to log file')
    sys.stdout = open(join(model_ws, 'log', 'model_connect_log_'+str(pd.to_datetime('today').date())+'.txt'), 'w')

print(model_ws)
print(scenario_name)
print('\n\n')

# print out additional summary info
print('SW Constraint %.2f' %sw_con, 'inches')
print('GW Constraint %.2f' %gw_con, 'inches')
print('\n')

swb_ws = join(model_ws, 'rep_crop_soilbudget')
os.makedirs(join(swb_ws, 'field_SWB'), exist_ok=True)
os.makedirs(join(swb_ws, 'output'), exist_ok=True)
# secondary directory for re-run with adjusted irrigation
os.makedirs(join(model_ws, 'crop_soilbudget_est', 'field_SWB'), exist_ok=True)

# simple code to set dates for april 1
all_run_dates = pd.read_csv(join(model_ws, 'crop_modflow', 'all_run_dates.csv'), parse_dates=['date'])

# %% [markdown]
# Review season dates, the plan was to change the start of the irrigation date for misc. grain and hay since irrigators don't typically start until summer even though it grows in winter. Simplify to crop choice Jan 1 and Irrig. run Apr 1
# - also at some point we may have discussed just doing them at the same time? in that case we just need to iterate over years

# %%
print(all_run_dates)
print('\n\n')
print('Starting yearly iterations')
print('\n\n')
sys.stdout.flush()

# %% [markdown]
# Cameron asked about cross-validation and hind-sight in the model. Maybe look at best way to re-run SWB or modflow after seeing the model run.

# %% [markdown]
# # Begin yearly iteration here 

# %% [markdown]
# For the first period we want to run things simply to get started. Thus we won't have output from the irrigation model.
# - check which soil budget data is available from the pre-existing dataset or if we should set up something simple to run on the hydrology of the current year assuming full irrigation on rain-ETc deficit.
#     - this should be a question for Yusuke, in the off-season of October to March we are going to assume no irrigation so there will no groundwater pumping. Just dispersed recharge.

# %%
# # this loop was set to run for the years of interest
# start at 1 instead of 0 to skip first pre-period
for m_per in np.arange(1, all_run_dates.shape[0]-1):
# for m_per in [1]:
# for m_per in [4]:

    m_strt = all_run_dates.iloc[m_per].date
    m_end = all_run_dates.iloc[m_per+1].date-pd.DateOffset(days=1)
    use = all_run_dates.iloc[m_per].use
    dates = pd.date_range(m_strt, m_end)

    # %%
    print(m_strt.date(),'to', m_end.date())

    # %%
    year = m_strt.year
    # crop='Corn'
    # define whether the year for loading variables should change with year or be static (use average)
    # year_load_var = None
    if year_load_var_in == "False":
        year_load_var = year
    else:
        year_load_var = year_load_var_in
    print('Using',str(year_load_var_in),'for crop variables')

# %% [markdown]
# # Define groundwater elevation sampling
#
# **The stata coefficients were estimated using DTW in fall that is in meters NOT in feet**  
#
# **Fall DTW is for crop choice and should be DTW in meters**

    # %%
    model_ws_last = join(model_ws, 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date()))
    hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')

    # %%
    # model output dates
    m_dates = all_run_dates.loc[m_per-1].date+np.array(hdobj.get_times()).astype('timedelta64[D]')
    m_dates = pd.DataFrame(m_dates, columns=['dates']).set_index('dates')
    m_dates['kstpkper'] = hdobj.get_kstpkper()
    # subset to 1 month of output
    # determine dates for fall sampling
    fall_dates = m_dates[m_dates.index.month==10]



# %% [markdown]
# Irrigation decisions are made with spring levels but the crop choice model was developed with fall levels

    # %%
    # get head value from last 30 days to avoid using extreme single day value
    fall_heads = avg_heads(fall_dates.kstpkper.values, hdobj, m_dim)
    
    # the dtw conversion runs a little slow
    # get the DTW for the wels in the simulation from the last period
    well_dtw = sample_dtw(fall_heads, parcel_wells)
    # need to make integer for join with crop choice
    well_dtw.UniqueID = well_dtw.UniqueID.astype(int)

# %% [markdown]
# # Crop choice model
# This uses the updated DTW from each previous year.
#
# Some of the file load in of the crop choice could be moved external of the loop (e.g., WY types and logit_coefs)  
#
# Need to also update to use the estimated revenue rather than historical values.

    # %%
    # load Sac Valley WYT
    wyt = pd.read_csv(join(proj_dir, 'model_inputs', 'sacramento_WY_types.txt'))
    # whether it is a critical or dry year or not
    wyt['wy_dc'] = 0
    # for wet or above normal years the boolean will be 0 
    wyt.loc[wyt['Yr-type'].isin(['C','D']),'wy_dc'] = 1

    # %%
    crop_choice_dir = 'parcelchoicemodelupdate'
    # Read logit coefficients
    logit_coefs = pd.read_csv(join(crop_choice_dir, 'logit_coefs.csv'))
    logit_coefs = pd.read_csv(join(crop_choice_dir, 'logit_coefs_last_crop.csv'))
    crop_dict = logit_coefs.set_index('Crop_Eq_new')['Crop_Eq'].to_dict()

    # the parcel data needs the dtwfa (avg fall dtw in feet for the parcel) and wy_dc (pulled from Sac wy type dataset and switched to dry boolean)
    # missing WY type prediction? 
    # Read parcel data
    data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
    # still needs to be updated to auto update DTW and WY type
    # data['wy_dc'] = np.where(data['year'] == 2020, 1, 0) # should be pulled from Sac WY type data instead
    data['wy_dc'] = wyt.loc[wyt.WY==year, 'wy_dc'].values[0]
    # update DTW to use simulated value instead of contoured
    well_dtw_merge = well_dtw[['UniqueID','dtw_ft']]
    # the DTW used in crop choice should be in meters not ft
    well_dtw_merge['dtwfa'] = well_dtw_merge.dtw_ft * 0.3048

    data = data.drop(columns=['dtwfa', 'dtwsp'])
    data = data.merge(well_dtw_merge, left_on='parcel_id',right_on='UniqueID')
    
    # Import prior year revenue data by crop
    rev_prior_yr_df = pd.read_csv(join(crop_choice_dir, "data_model/rev_prior_yr.csv"))

    # TO DO: Yusuke wants this to be to the nearest year
    rev_prior_yr_all = pd.read_csv(join(crop_choice_dir, "data_model/rev_prior_yr_all.csv"))
    rev_prior_yr_all = rev_prior_yr_all.rename(columns={'crop_cat':'Crop_Eq'})
    # identify nearest year to use revenue
    # should do the same for crop parameters and use year_load_var
    rev_prior_yr_all['year_offset'] = (rev_prior_yr_all.year - year_load_var).abs()
    # select the year that is closest in and use the lower year in a tie.
    rev_prior_yr_df = rev_prior_yr_all.sort_values(['year_offset', 'year']).drop_duplicates('Crop_Eq', keep='first')
    rev_prior_yr_df = rev_prior_yr_df.reset_index(drop=True).drop(columns=['year_offset', 'year'])


    # %%
    # expect updated column name for pod
    # add column to POD that was available in previous dataset
    pod = data[['parcel_id','pod']].copy().rename(columns={'pod':'pod_bool'})
    pod['pod'] = 'No Point of Diversion on Parcel'
    pod.loc[pod.pod_bool==1, 'pod'] = 'Point of Diversion on Parcel'


    # %%
    # need to add columns to data to inform last years crop from loop
    if m_per == 1:
        # the STATA data has results for 3726 fields, not the full 4,000 fields
        # stata_first = pd.read_csv(join(proj_dir, 'model_inputs', 'stata_input_2018.csv'))
        # the full 4,000 fields comes from a combination of years using the first year where a crop is specified
        print('Using stata input initial')
        stata_first = pd.read_csv(join(proj_dir, 'model_inputs', 'stata_input_initial.csv'))
        # where there is not a crop_cat_last, assume that crop_cat is valid
        stata_first.loc[stata_first.crop_cat_last.isna(), 'crop_cat_last'] = stata_first.loc[stata_first.crop_cat_last.isna(), 'crop_cat']
        crop_cat_last = stata_first[['parcel_id','crop_cat_last', 'chosen']].replace(crop_dict)
        crop_cat_last = crop_cat_last.pivot_table(index='parcel_id',columns='crop_cat_last', values='chosen', fill_value=0)
        # if a field is missing from the table then assume it is unclassified fallow?

    else:
        crop_cat_last =  pd.read_csv(join(swb_ws, 'field_SWB', 'parcel_crop_choice_'+str(year-1)+'.csv'), index_col=0)
        crop_cat_last = crop_cat_last.assign(chosen=1).replace(crop_dict)
        crop_cat_last = crop_cat_last.pivot_table(index='parcel_id',columns='Crop_Choice', values='chosen', fill_value=0)


    # %%
    # crop choice model uses "_" instead of " "
    # the irrigation model set up uses " "
    # data_out = predict_crops(data.copy(), rev_prior_yr_df, logit_coefs)
    data_out = predict_crops_rand(data.copy().merge(crop_cat_last.reset_index()), 
                                  rev_prior_yr_df, logit_coefs)
    data_out['Crop_Choice'] = data_out.Crop_Choice.str.replace('_',' ')
    # update naming of Corn
    data_out.Crop_Choice = data_out.Crop_Choice.str.replace('Corn  ','Corn, ')
    # we should include the DTW data point used for comparison
    # add depth to water information to the output
    data_out = data_out.merge(well_dtw[['UniqueID','dtw_ft']].rename(columns={'UniqueID':'parcel_id'}))
    # TODO static crop/AW: if skipping crop choice, would not write and instead read in this 
    # save output with only parcel and crop choice
    data_out.to_csv(join(swb_ws, 'field_SWB', 'parcel_crop_choice_'+str(year)+'.csv'))


    # %%
    # for the first year of crop prediction copy this backward three years (assume that it was static)
    # to allow for establishment irrigation code to run
    if m_per ==1:
        for n in np.arange(1,4).astype(int):
            data_out.to_csv(join(swb_ws, 'field_SWB', 'parcel_crop_choice_'+str(year-n)+'.csv'))

# %% [markdown]
# # Load MF output
# The depth to water function should sample from the previous model run which may be a year or less long.
#
# The original get_dtw function was set up assuming a continuous model run which won't be the case.
#

    # %%
    # determine dates for spring sampling
    spring_dates = m_dates[m_dates.index.month==3]
    # get head value from last 30 days to avoid using extreme single day value
    spring_heads = avg_heads(spring_dates.kstpkper.values, hdobj, m_dim)
    
    # the dtw conversion runs a little slow
    # get the DTW for the wels in the simulation from the last period
    well_dtw = sample_dtw(spring_heads, parcel_wells)
    # need to make integer for join with crop choice
    well_dtw.UniqueID = well_dtw.UniqueID.astype(int)
    well_dtw.to_csv(join(swb_ws,'field_SWB', 'modflow_spring_dtw_ft_WY'+str(year)+'.csv'))

# %% [markdown]
# # Irrigation submodel

# %% [markdown]
# The crop choice model didn't predict any corn.
# - this submodel may need to be in a separate script to run multiprocessing (parallel) as this re-loads the entire active script each time.
# - this would involve calling a batch file that then call the python file to do multiprocess of the python script
# - the alternative would be to simply pre-define the SWB for the maximum possible range of DTW for each year then they would only need to be called instead of re-optimized.

    # %%
    crop_list = ['Corn','Alfalfa','Pasture','Misc Grain and Hay', 'Grape']
    # this is set up to essentially skip over Unclassified Fallow since we don't want to calculate the depth to water for it
    # this is also skipping Other which should be accounted for by pulling the output from the original model
    # for the crop, or should this be something else?


    # %%
    #
    crop_in = data_out.merge(pod)
    crop_in = crop_in.rename(columns={'Crop_Choice':'name'})
    crop_in.to_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv')) # not sure why I save it again just with POD info
    # crop_in[crop_in.name==crop]
    pred_crops = crop_in.name.unique()
    print(pred_crops)

    # add in the crops from previous years
    for year_previous in np.arange(all_run_dates.date.dt.year.min()-2, year):
        crop_in_previous = pd.read_csv(join(swb_ws, 'field_SWB', 'parcel_crop_choice_'+str(year_previous)+'.csv'), index_col=0)
        # to avoid conflict, add previous years as new columns with convention name_year
        crop_in_previous = crop_in_previous[['parcel_id','Crop_Choice']].rename(columns={'Crop_Choice':'name_'+str(year_previous)})
        crop_in = crop_in.merge(crop_in_previous)

    # %%
    print(crop_in.groupby('name')[['parcel_id']].count())
    # paper by Snyder, statewide values for pasture or range, ended up using pasture as 
    # crop coefficients are either 0.9 or 0.95

# %% [markdown]
# ## Simplified representation of DTW range from min to max.
# Instead of using the full record of the DTW, the mdoel should just sample the average DTW for the month of interest

    # %%
    ## simple representative DTW for linear steps from dtw_min to dtw_max
    ## with a 5 ft decline from June to December based on observed data
    dtw_simple_df = calc_simple_dtw(well_dtw, year, dtw_step = 10)
    dtw_simple_df.to_csv(join(swb_ws,'field_SWB', 'dtw_ft_WY'+str(year)+'.csv'))


# %% [markdown]
# ## Iterate over all crops to save the representative profiles
# The way the irrigation models are set up, they still run for the entire season but they are simulated in spring which means that they can't decide winter water pumping for the first year. In each year following the simulation runs april to april so they can influence pumping in the following year.
# - just need to note that the first period of 6 months is spin-up in a way.

# %% [markdown]
# After updating function need to check everything still works.
#
# - fixed issue hdf5 overwrite
# - found that Corn, sorghum name match was off

    # %%
    print(pred_crops)
    for crop in crop_list:
    # for crop in ['Alfalfa']:
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year = year_load_var, input_name=input_name)
        # need to account for when crops aren't predicted and skip them
        # if pred_dict[crop] in pred_crops: 
        print(crop, ':',pred_dict[crop])

    # %%
    # initialize HDF5 files for the year
    # if skipping SWB optimization then don't need to intialize these
    # initialize SWB folder
    if create_rep_swb in ['existing_local', 'existing_local_and_crop', 'existing_shared']:
        print('Using existing hdf5 files for SWB optimization')
    else:
        print('Initiating new hdf5 files for SWB optimization')
        for var in ['profit', 'cost', 'yield', 'percolation','GW_applied_water', 'SW_applied_water']:
            name = join(swb_ws, 'field_SWB', var + '_WY'+str(year)+'.hdf5')
            init_h5(name)
    
        for var in ['swb_output']:
            name = join(swb_ws, 'field_SWB', var + '_WY'+str(year)+'.hdf5')
            init_h5(name, groups=['wc','ETa', 'rp'])

# %%
    # if skipping SWB optimize then skip this
    if create_rep_swb in ['existing_local', 'existing_local_and_crop', 'existing_shared']:
        print('Using existing optimized SWB  results')
    else:
        print('Calculating new optimized SWB results')
        # for crop in ['Alfalfa']:
        for crop in crop_list:
            # will need to add year to swb.load_var(crop, year) if we want to use year specific profit and cost
            # variables, this would be useful for comparing against baseline while future should use average
            # within load_run_swb, these variables are called again specified by year
            var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year = year_load_var, input_name=input_name)
            # need to account for when crops aren't predicted and skip them
            if pred_dict[crop] in pred_crops: 
                # to equalize the situation we might use a simple DTW profile
                load_run_swb(crop, year, crop_in, swb_ws,
                             dtw_simple_df, 
                             soil_rep=True,
                             sw_con=sw_con, gw_con=gw_con,
                             input_name = input_name
                             ) 
        sys.stdout.flush()

    # %%
    fn = join(swb_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
    print('Crops with SWB results')
    with h5py.File(fn) as dset:
        finished_crops = list(dset['array'].keys())
        # arr = dset['array']['Corn'][:]
        print(finished_crops)
    # only grape was completed?

# %% [markdown]
# 1. Load the representative results and sample for each field by crop type to back calculate the irrigation requirements. use the estimated irrigation as an input to the modflow model for pumping and percolation for recharge.
#     -   use the DTW id to reference to the irrigation in the full array, need to group by SW, GW or mixed.
#     -   if we wanted we could re-run the SWB one time with the specified irrigation rates to get the exact recharge rates with field specific values (this should be done, and is partly implemented from when re-running for establishment irrigation)
# 2. Run the modflow model to get the resultant DTW profile
# 3. re-calculate the profit using the irrigation and actual DTW profile on a soil by soil basis (non-optimization) after running the next modflow chunk. Actually the re-run for the true profit could be done if profits aren't needed mid-simulation
#

    # %%
    # TODO static crop/AW: if skipping crop optimize then should skip this
    # load the processed dataframe with all datas
    # this is where we could specify the pre-calculated optimized SWB output
    # instead of swb_ws specify the path on box
    swb_version = input_name.replace('static_model_inputs','').replace('.xlsx','')
    if year_load_var_in != "False":
        swb_version += '_'+year_load_var
    swb_rep_ws = join(proj_dir, 'model_inputs', 'swb_rep', 'version'+swb_version)

    pc_df_all, irr_gw_df_all, irr_sw_df_all = get_wb_by_parcel(swb_ws, year, 
                     crop_in, finished_crops, dtw_simple_df, well_dtw)
    # this output with the parcel data needs to be saved as well - now done after updating below
    # pc_df_all.to_csv(join(swb_ws, 'output', 'pc_all'+str(year)+'.csv'))
    # irr_gw_df_all.to_csv(join(swb_ws, 'output', 'irr_gw_all'+str(year)+'.csv'))
    # irr_sw_df_all.to_csv(join(swb_ws, 'output', 'irr_sw_all'+str(year)+'.csv'))


# %% [markdown]
# # Correct for new crops and apply estimated irrigation to actual field conditions for true percolation
#
# - For irrigation this involves scaling the currently proposed irrigation rates by the new total divided by the current total.  
# - For percolation this involves re-running the SWB model for all fields to get new percolation rates.  
# - For yield, we need to then specify no yield for these years and 1/2 yield for Grape year 3. This could maybe be done in post-processing? Or by reading in hdf5 file and manually correcting based on year and field ID in the loop (I forgot to update this originally)
#
# for simplicity take predicted crops from year 1 and copy backward three years to ease input set up (assume first prediction has no changes from past year)

    # %%
    import functions.summarize_functions
    reload(functions.summarize_functions)
    # functions related to post-processing updated irrigation
    from functions.summarize_functions import format_irr_all, adj_irr_rates

    # %%
    # TODO static crop/AW: if skipping previous step then need to skip this as well
    # update irrigation rates to account for establishment where crops change
    # does it make sense that these are addressed separately? In the chance there is mixed irrigation then
    # the total scaling would need to be adjusted
    irr_gw_df_all = adj_irr_rates(irr_gw_df_all, year, crop_list, 
                                  irr_est, crop_in, pred_dict, report_change=True)
    # quick way to use existing function to identify change and years back
    crop_change_info = irr_gw_df_all[['UniqueID','changed','year_offset']].drop_duplicates()
    irr_gw_df_all = irr_gw_df_all.drop(columns=['changed','year_offset'])
    irr_sw_df_all = adj_irr_rates(irr_sw_df_all, year, crop_list, 
                                  irr_est, crop_in, pred_dict)
    # save files for referencing when re-running the SWB to calculate the actual profit after modflow runs
    irr_gw_df_all.to_csv(join(swb_ws, 'output', 'irr_gw_all'+str(year)+'.csv'))
    irr_sw_df_all.to_csv(join(swb_ws, 'output', 'irr_sw_all'+str(year)+'.csv'))
    # save crop change info for reference in updating profit/yield
    crop_change_info.to_csv(join(swb_ws, 'output', 'crop_change_info'+str(year)+'.csv'), index=False)

    # %%
    # need to get dtw_df (even if just estimated), unlike f_summarize_output modflow hasn't been run yet
    # need data from November (-1 year) to December
    # dtw_df doesn't really matter as no opt is run, and the profit will need to be re-run with actually dtw later
    dtw_wide = np.repeat(np.reshape(well_dtw.dtw_ft.values, (1,-1)), len(dtw_simple_df), axis=0)
    # assume same simple 5 ft dtw decline
    dtw_decline = dtw_simple_df.iloc[:,0]-10
    dtw_decline = np.repeat(np.reshape(dtw_decline.values, (-1,1)), len(well_dtw), axis=1)
    dtw_wide_df = pd.DataFrame(dtw_wide-dtw_decline)
    dtw_wide_df.index = dtw_simple_df.index
    dtw_wide_df.columns = well_dtw.UniqueID

    # %%
    # initiate hdf5 files to save outputs
    for var in ['profit', 'yield', 'cost','percolation','GW_applied_water', 'SW_applied_water']:
        name = join(model_ws,'crop_soilbudget_est', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        init_h5(name)
    # should review and decide if it helps to save this separately or if it is better to overwrite
    # original files to avoid confusion

    # %%
    # check profit/yield if annual or daily, they are annual
    # added hdf5 to track SWB output
    for var in ['swb_output']:
        name = join(model_ws,'crop_soilbudget_est', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        init_h5(name, groups=['wc','ETa', 'rp'])

    # %%
    print('Re-running SWB with updated irrigation values for fields with establishment irrigation rates')
    print('Calculates actual percolation for fields')
    # for crop in ['Alfalfa']:
    for crop in crop_list:
        # will need to add year to swb.load_var(crop, year) if we want to use year specific profit and cost
        # variables, this would be useful for comparing against baseline while future should use average
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year = year_load_var, input_name=input_name)
        # need to account for when crops aren't predicted and skip them
        if pred_dict[crop] in pred_crops:     
            # irr_all is the combination for each crop
            irr_all = format_irr_all(crop, year, crop_in, pred_dict, irr_gw_df_all, irr_sw_df_all)
            # from code f_summarize_output that re-runs SWB at end of each year with true DTW for each field
            # should rename the path so we can have both saved
            load_run_swb(crop, year, crop_in, join(model_ws,'crop_soilbudget_est'), dtw_wide_df, input_name = input_name, soil_rep = False,
                            run_opt=False, irr_all=irr_all, field_id = 'parcels')

            sys.stdout.flush()

    # creates the updated irrigation rates so that percolation can be re-calculated
    # next need to extract out the water budget data to create the updated percolation data

    # %%
    # may want to use code from plot_output that has read_hdf5  array and assigns UniqueID
    # irrigation rates are already defined for establishment so only need to read in percolation
    df_all = pd.DataFrame()
    for var in ['percolation']: # ,'GW_applied_water', 'SW_applied_water'
        print(var, end=',')
        name = join(model_ws, 'crop_soilbudget_est', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        for crop in crop_list:
            var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year = year_load_var, input_name=input_name)
            yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
            yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
            # get the total extent of the irrigation season (calculation period)
            crop_dates = pd.date_range(yield_start.min(), yield_end.max(), freq='D')
            # need to account for when crops aren't predicted and skip them
            if pred_dict[crop] in pred_crops: 
                # extract output and convert to dataframe with ID columns
                arr = read_crop_arr_h5(crop, name)
                df = pd.DataFrame(arr, columns=crop_dates)
                # add parcel information back
                df = pd.concat((df,crop_in[crop_in.name==pred_dict[crop]].reset_index(drop=True)),axis=1)
                # melt to long format for easier appending
                df = df.melt(var_name='date', id_vars=crop_in.columns)
                df = df.assign(crop=crop, year=year, var=var)
                df_all = pd.concat((df_all, df))

    # %%
    # columns for  pc_df_all 'UniqueID', 'dtw_id', 'date', 'rate' is dtw_id needed?
    pc_df_all = df_all.rename(columns={'parcel_id':'UniqueID', 'value':'rate'})[['UniqueID','date','rate']]
    pc_df_all.to_csv(join(swb_ws, 'output', 'pc_all'+str(year)+'.csv'))


    # %%
    # load in data on init_wc
    # the best way to save water budget output is likely to put them all in the same hdf5 file
    # should this be the crop_soilbudget (isn't updated until end) or crop_soilbudget_est (correct dims)? 
    fn = join(model_ws, 'crop_soilbudget_est', 'field_SWB', "swb_output_WY"+str(year)+".hdf5")
    init_wc_all = pd.DataFrame()
    for crop in crop_list:
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year = year_load_var, input_name=input_name)
        yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
        yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
        # get the total extent of the irrigation season (calculation period)
        crop_dates = pd.date_range(yield_start.min(), yield_end.max(), freq='D')
        wc_all = read_crop_arr_h5(crop, fn, group='wc')
        df = pd.DataFrame(wc_all, columns=crop_dates)
        # sample the water content from the final time step
        # init_wc = wc_all[-1,:]
        init_wc = df[[yield_end.max()]]
        # add info on UniqueID for ordering parcels
        init_wc = pd.concat((init_wc,crop_in[crop_in.name==pred_dict[crop]].reset_index(drop=True)[['parcel_id']]),axis=1)
        init_wc = init_wc.assign(crop=crop, year=year)
        init_wc_all = pd.concat((init_wc_all, init_wc.melt(id_vars=['crop','year','parcel_id'], var_name='end_date',value_name='init_wc')))
    # sort values to prepare for using as initial conditions for winter SWB
    init_wc_all = init_wc_all.sort_values('parcel_id').reset_index(drop=True)

    # for the ag SWB we also need to pull in the WC from fallow fields from the last date to start the next as there won't be irrig 
    # reuslts for these fields

# %%
# for testing in function
# init_wc_all.to_csv('temp_wc.csv', index=False)

# %% [markdown]
# ## RCH input
#
# Should explore to see if this is missing some of the extra recharge I include in the original model to account for floodplain inundation
    # %%
    print('Beginning MODFLOW input update')

    # %%
    # join water budget data to field area for scaling to cell for recharge
    rch_df = pc_df_all.merge(field_df)
    # calculate the fraction of a cell covered by the field
    rch_df['field_scale'] = rch_df.field_area/(200*200)
    # calculate the recharge rate for each field
    rch_df['rch_rate'] =  rch_df.rate*rch_df.field_scale
    # set date as index for sampling
    rch_df = rch_df.set_index('date')
    # we only need to assign rch when non-zero
    # rch_df = rch_df[rch_df.rch_rate!=0]
    # actually need to keep the zero values in place to avoid being overwritten by "winter swb" rates

# %%
# # TODO: import irrigation and recharge estimates for row, cell (or fields) marked as other for the year
# and apply input from original model

# %% [markdown]
# Take the combined dataframe to inform the recharge array with recharge from the irrigation season. The next step is to insert recharge from non-irrigated times and native lands. The ag recharge was originally defined for the WY which meant recharge during the period of interest was missing, updated to run the same period as each year 4/1 to 3/31 next year.


# %%
# sample ending water content from irrigation season

    # %%
    ## add off-season (winter recharge) farmlands recharge 
    ## need to use the crop-type selected for the year and assume field content at start
    ## nothing needs to be optimized so the soil water budgets can be run at once
    ## the only difference is that the ETc might need to be resampled

    # function that saves a csv with the output by unique field similar to the native land
    ag_pc_df = run_swb_ag_winter(year, m_nam = m_nam, loadpth = loadpth, init_wc = init_wc_all)

    # long format to join with row,column
    ag_pc_df_long = ag_pc_df.melt(ignore_index=False, 
                            var_name='UniqueID',
                           value_name='pc')
    ag_pc_df_long.UniqueID=ag_pc_df_long.UniqueID.astype(int)

    # join field to row,column
    ag_pc = ag_pc_df_long.reset_index().merge(grid_soil_ag)
    ag_pc = ag_pc.rename(columns={'index':'date'})
    # calculate recharge rate by scaling for fraction of cell covered 
    ag_pc['rch_rate'] = ag_pc.pc*ag_pc.field_area/(200*200)
    # update index to datetime
    ag_pc = ag_pc.set_index('date')
    ag_pc.index=pd.to_datetime(ag_pc.index)
    # drop dates and cells with zero recharge as the default recharge is 0
    ag_pc = ag_pc[ag_pc.rch_rate!=0]


# %% [markdown]
# Identify where the irrigated period SWB already provides percolation data and only provide in fill for the outer period. For fields where the land use is not defined or fallow then this will fill in the entire period. For fields with the irrigation SWB it should only fill in outside of the irrigation period.

    # %%
    # this code seems to not be fully working because there are overlap dates for the same field
# I accidentally included an inverse boolean term on the right_only statement
    # need to only keep dates with percolation where there is not already data from the optimized swb
    outer_join = rch_df[['UniqueID']].reset_index().merge(ag_pc[['UniqueID']].reset_index(), how = 'outer', indicator = True)
    # identify data that is on the right side only (ag winter recharge) 
    ag_join = outer_join[(outer_join._merge == 'right_only')].drop('_merge', axis = 1).drop_duplicates()
    
    # filter out the ag winter recharge to join to the dataframe
    ag_pc_include = ag_pc.reset_index().merge(ag_join, how='inner')
    # make index the date
    ag_pc_include = ag_pc_include.set_index('date')

    # %%
    ## add native lands recharge and winter ag land precip recharge
    rch_all = pd.concat((rch_df, native_pc, ag_pc_include))

    ## add recharge due to MAR or other scenario
    # have dataframe with recharge rates for cells and row,column
    mar_grid_join = mar_grid[mar_grid.index.isin(dates)]
    rch_all = pd.concat((rch_all, mar_grid_join))

    # aggregate to the row, column level
    rch_all = rch_all.groupby(['date','row','column'])[['rch_rate']].sum()
    rch_all = rch_all.reset_index(['row','column'])
    # is groupby faster or just adding individual datasets to the array??


    # %%
    # allocate dictionary for recharge
    rch_dict = dict()
    for t, d in enumerate(dates):
        # get data for the stress period
        rch_arr = np.zeros((m.dis.nrow, m.dis.ncol))
        if d in rch_all.index:
            # spd_df = rch_df.loc[d]
            spd_df = rch_all.loc[d]
            rch_arr[spd_df.row-1, spd_df.column-1] = spd_df.rch_rate
        # assign values to recharge dict
        rch_dict[t] = rch_arr
    # print to help ID errors
    print('Finished updating RCH dictionary')
# %% [markdown]
# ## WEL input

# %%
    
    # join water budget data to field area for application to cell for well
    wel_df = irr_gw_df_all.merge(well_loc_merge)
    # calculate the volumetric flux of the well 
    wel_df['flux'] = -wel_df.rate*wel_df.area_m2
    # # set date as index for sampling
    wel_df = wel_df.set_index('date')
    # wel_df

    # %%
    # allocate dictionary for recharge
    # I think the issue here with writing input was indentation was incorrect
    wel_dict = dict()
    for t, d in enumerate(dates):
        # get data for the stress period
        if d in wel_df.index:
            # spd_df = wel_df.loc[d]
            # need to use boolean to avoid ending up with single value (series instead of dataframe)
            spd_df = wel_df.loc[wel_df.index==d]
            wel_arr = spd_df[['k','i','j','flux']].values
        else:
            wel_arr = [[0,0,0,0]]
        # append domestic well pumping to pumping dataset
        dom_loc['flux'] = - dom_use.loc[d,'flux_m3d']*dom_loc.pump_scale
        wells_dom = dom_loc[['layer','row','column','flux']].values
        
        # got error that says wel_arr is 1D when wells_dom is 2D
        # when manually running the code this is because there is a date where there is only one well pumping
        # and it becomes a 1D array
        # assign input to well dictionary
        wel_dict[t] = np.append(wel_arr, wells_dom, axis=0)
    # print to help ID errors
    print('Finished updating WEL dictionary')

# %%
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# plt.imshow(rch_dict[100])
# plt.colorbar(norm = mpl.colors.LogNorm)

# %% [markdown]
# # Run update modflow

    # %%
    print('Loading in period MODFLOW model')
    # load the existing modflow model for the next year
    mf_ws = join(model_ws, 'crop_modflow/'+str(m_strt.date()))
    m_month = flopy.modflow.Modflow.load(join(mf_ws,'MF.nam'), model_ws = mf_ws, load_only=['DIS'])
    # m.model_ws = mf_ws
    # The model should start in hydraulic connection
    if m_per > 0:
        mf_ws_last = join(model_ws, 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date()))
        hdobj = flopy.utils.HeadFile(mf_ws_last + '/MF.hds')
        sp_last = hdobj.get_kstpkper()[-1]
        strt = hdobj.get_data(sp_last)
    
    # re-create basic file with last period heads as new starting heads
    bas_month = flopy.modflow.ModflowBas(model = m_month, ibound=ibound, strt = strt)
    # write bas input with updated starting heads
    bas_month.write_file()

    # update basic input to save flow budget summary
    update_bas_options(mf_ws, m_month.name, added_options = "BUDGETDB flow_budget.txt")

# %%

    # update flopy packages with new data
    wel_month = flopy.modflow.ModflowWel(model=m_month, stress_period_data = wel_dict, ipakcb=55)
    rch_month = flopy.modflow.ModflowRch(model=m_month, nrchop = 3, rech =  rch_dict, ipakcb=55)
    
    wel_month.write_file()
    rch_month.write_file()
    print('Finish writing new modflow input')
    


    # %%
    # need to specify the path of the modflow exe or it defaults to mf2005 and crashes
    m_month.exe_name = 'mf-owhm.exe'
    success, buff = m_month.run_model()

    # %%
    sys.stdout.flush()


# %% [markdown]
# # Re-run soil water budget to calculate actual profit
# The soil water budget needs to be re-run at the end of each year to calculate the profit based on the actual water availability for surface water and groundwater. This could be done at the very end of the simulation if the farmer doesn't need previous year profits between years.
# easiest way is likely to call summarize_output as a function or by having it read-in a csv with the year to make a pseudo-function
#

    # %%
    # calculate the actual profit and yield for each parcel then calculate the average
    # to inform next years crop choice
    summarize_output_year(loadpth, m_nam, m_per, parcels, input_name)

    # would need to update yield within this final function to account for establishment irrigation

# %%
t_final = time.time()
print('Total time was %.2f hours' %((t_final-t_start)/3600))

# %% [markdown]
# The MODFLOW run-time for one year was 49 min. This was for the version with 19 layers which we can cut down to 10 or so by removing TPROGs and this could be sped up by reducing the amont of output saved since we may not need to review the CBC output.

