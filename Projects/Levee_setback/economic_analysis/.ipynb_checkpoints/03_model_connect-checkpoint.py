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
# Script to run the code for the crop choice model, irrigation optimizer and modflow in a cycle for the years of interest.
#
# 1. Run the crop choice model
# 2. Run the soil water budget optimization
# 3. Update WEL/RCH packages then run MF
# 4. Start the next year

# %%
import h5py
import numpy as np
import pandas as pd
import geopandas as gpd
import os
from os.path import join,exists, dirname, basename, expanduser

# %%
from importlib import reload

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')

# %%
import sys
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

# %%
from parcelchoicemodelupdate.f_predict_landuse import predict_crops

# %%
import functions.Basic_soil_budget_monthly as swb
reload(swb)
# import f_swb_profit_opt
# reload(f_swb_profit_opt)
# from f_swb_profit_opt import load_run_swb

import f_rep_swb_profit_opt
reload(f_rep_swb_profit_opt)
from f_rep_swb_profit_opt import load_run_swb

# %% [markdown]
# If we are going to run the crop/swb for a certain period then there should already have been a MODFLOW model run for that same period with estimates of streamflow and precipitation to drive the other boundary conditions.

# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

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
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Regional/'

# base_model_ws = loadpth + 'crop_soilbudget'
m_model_ws = loadpth + 'historical_simple_geology_reconnection'

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
# well_loc_merge.depth_m
## sample well layer
# get top elevations to calculate well screen elevation
well_loc_merge['top_elev'] = m.dis.top.array[well_loc_merge.row-1, well_loc_merge.column-1]
well_loc_merge['well_elev'] = well_loc_merge.top_elev - well_loc_merge.depth_m
# make layer 1 based
well_loc_merge['layer'] = get_layer_from_elev(well_loc_merge.well_elev.values, 
                    botm[:, well_loc_merge.row-1, well_loc_merge.column-1], m.dis.nlay)+1
# make 0-based columns
well_loc_merge[['k','i','j']] = well_loc_merge[['layer','row','column']]-1

# %% [markdown]
# Specify dates for the entire model period and seasonal dates to start/stop the model.

# %%
all_strt_date = pd.to_datetime(m.dis.start_datetime)
all_dates = all_strt_date + (m.dis.perlen.array.cumsum()-1).astype('timedelta64[D]')
all_end_date = all_dates[-1]
print(all_strt_date, all_end_date)
months = pd.date_range(all_strt_date, all_end_date, freq='MS')
years = pd.date_range(all_strt_date, all_end_date, freq='YS').year.values

# %%

# load summary excel sheet on irrigation optimization
# this will specify the date ranges to run and pause
fn = join(data_dir,'static_model_inputs.xlsx')
season = pd.read_excel(fn, sheet_name='Seasons', comment='#')

# %%
# choose crops on first day of year
month_crop = pd.Series(1)
day_crop = pd.Series(1)
month_irr = pd.Series(4)
day_irr = pd.Series(1)

## specify dates where modflow will start 
all_run_dates = pd.DataFrame()
# yn = 0
# y = years[yn]
for y in years:
    run_dates = swb.ymd2dt(y, season.month_run, season.day_run, season.start_adj)
    run_dates = run_dates.drop_duplicates().sort_values()
    run_dates = pd.DataFrame(run_dates).assign(use='irrigation')
    crop_date = swb.ymd2dt(y, month_crop, day_crop, season.start_adj)
    crop_date = pd.DataFrame(crop_date).assign(use='crop').dropna()
    irr_date = swb.ymd2dt(y, month_irr, day_irr, season.start_adj)
    irr_date = pd.DataFrame(irr_date).assign(use='irrigation').dropna()
    # all_run_dates = pd.concat((all_run_dates, crop_date, run_dates))
    all_run_dates = pd.concat((all_run_dates, crop_date, irr_date))

# all_run_dates
# simple code to set dates for april 1
all_run_dates = pd.date_range(all_strt_date, all_end_date,freq='AS-Apr')
all_run_dates = pd.DataFrame(all_run_dates).assign(use='irrigation')
# add total start and end dates
all_run_dates = pd.concat((pd.DataFrame([all_strt_date]).assign(use='start'), all_run_dates))
all_run_dates = pd.concat((pd.DataFrame([all_end_date]).assign(use='end'), all_run_dates))
all_run_dates=all_run_dates.sort_values(0).reset_index(drop=True).rename(columns={0:'date'})

# %% [markdown]
# Review season dates, the plan was to change the start of the irrigation date for misc. grain and hay since irrigators don't typically start until summer even though it grows in winter. Simplify to crop choice Jan 1 and Irrig. run Apr 1
# - also at some point we may have discussed just doing them at the same time? in that case we just need to iterate over years

# %%
all_run_dates

# %%
# # this loop was set to run for the years of interest
# for m_per in np.arange(1,5): # runs first year to next crop choice
# # for m_per in np.arange(0, all_run_dates.shape[0]-1):
for m_per in [1]:
    m_strt = all_run_dates.iloc[m_per].date
    m_end = all_run_dates.iloc[m_per+1].date-pd.DateOffset(days=1)
    use = all_run_dates.iloc[m_per].use
    dates = pd.date_range(m_strt, m_end)

# %%
print(m_strt.date(),'to', m_end.date())

# %%
# the model will run the irrigation optimizer on specified dates (multiple crops can be done at once or in sequence)
# the modflow model will be run for the periods between the specified irrigation optimizer dates
# loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Regional/'
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Regional/'

base_model_ws = loadpth + 'crop_soilbudget'
os.makedirs(base_model_ws, exist_ok=True)

# %%

# %%
# if use == 'start':
    # if use is start then just run the model for the period and stop

# %% [markdown]
# # Begin yearly iteration here 

# %% [markdown]
# For the first period we want to run things simply to get started. Thus we won't have output from the irrigation model.
# - check which soil budget data is available from the pre-existing dataset or if we should set up something simple to run on the hydrology of the current year assuming full irrigation on rain-ETc deficit.
#     - this should be a question for Yusuke, in the off-season of October to March we are going to assume no irrigation so there will no groundwater pumping. Just dispersed recharge.

# %%
year = 2015
year = m_strt.year
# crop='Corn'

# %% [markdown]
# # Define groundwater elevation sampling

# %%
model_ws_last = loadpth + 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date())
hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')
# heads of the very last step to start the next period
sp_last = hdobj.get_kstpkper()[-1]
end_heads = hdobj.get_data(sp_last)

# %%
# model output dates
m_dates = all_run_dates.loc[m_per-1].date+np.array(hdobj.get_times()).astype('timedelta64[D]')
m_dates = pd.DataFrame(m_dates, columns=['dates']).set_index('dates')
m_dates['kstpkper'] = hdobj.get_kstpkper()
# subset to 1 month of output
# determine dates for fall sampling
fall_dates = m_dates[m_dates.index.month==10]
# determine dates for spring sampling
spring_dates = m_dates[m_dates.index.month==3]


# %% [markdown]
# Irrigation decisions are made with spring levels but the crop choice model was developed with fall levels

# %%
import functions.f_gw_dtw_extract
reload(functions.f_gw_dtw_extract)
from functions.f_gw_dtw_extract import sample_dtw, avg_heads

# %%
# get head value from last 30 days to avoid using extreme single day value
fall_heads = avg_heads(fall_dates.kstpkper.values, hdobj, m)

# the dtw conversion runs a little slow
# get the DTW for the wels in the simulation from the last period
well_dtw = sample_dtw(fall_heads, botm)
# need to make integer for join with crop choice
well_dtw.UniqueID = well_dtw.UniqueID.astype(int)

# %% [markdown]
# # Crop choice model
# Evetnaully this should use the updated DTW from each previous year.

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

# the parcel data needs the dtwfa (avg fall dtw in feet for the parcel) and wy_dc (pulled from Sac wy type dataset and switched to dry boolean)
# missing WY type prediction? 
# Read parcel data
data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
# still needs to be updated to auto update DTW and WY type
# data['wy_dc'] = np.where(data['year'] == 2020, 1, 0) # should be pulled from Sac WY type data instead
data['wy_dc'] = wyt.loc[wyt.WY==year, 'wy_dc'].values[0]
# update DTW to use simulated value instead of contoured
well_dtw_merge = well_dtw[['UniqueID','dtw_ft']].rename(columns={'dtw_ft':'dtwfa'})
data = data.drop(columns=['dtwfa', 'dtwsp'])
data = data.merge(well_dtw_merge, left_on='parcel_id',right_on='UniqueID')

# Import prior year revenue data by crop
rev_prior_yr_df = pd.read_csv(join(crop_choice_dir, "data_model/rev_prior_yr.csv"))

# %% [markdown]
# Of note is that the simulated DTW for fall is 4x time greater than the observed so we may see a sharp drop in crops because of this.
# - it leads to an increase in unclassified fallow

# %%
# expect updated column name for pod
# add column to POD that was available in previous dataset
pod = data[['parcel_id','pod']].copy().rename(columns={'pod':'pod_bool'})
pod['pod'] = 'No Point of Diversion on Parcel'
pod.loc[pod.pod_bool==1, 'pod'] = 'Point of Diversion on Parcel'

# %%
# crop choice model uses "_" instead of " "
# the irrigation model set up uses " "
data_out = predict_crops(data.copy(), rev_prior_yr_df, logit_coefs)
data_out['Crop_Choice'] = data_out.Crop_Choice.str.replace('_',' ')
# update naming of Corn
data_out.Crop_Choice = data_out.Crop_Choice.str.replace('Corn  ','Corn, ')

data_out.to_csv(join(crop_choice_dir, 'parcel_crop_choice_'+str(year)+'.csv'), index=False)

# %% [markdown]
# # Load MF output
# The depth to water function should sample from the previous model run which may be a year or less long.
#
# The original get_dtw function was set up assuming a continuous model run which won't be the case.

# %%
# ## simple soil water budget without optimization to calculate
# ## recharge rates when there is no irrigation
# ## this still misses fields that are considered native I think

# # soil budget expects DTW array with same shape as soil
# dtw_avg = pd.DataFrame(pd.date_range(str(year-1)+'-11-1', str(year)+'-12-31'), columns=['date'])
# # create a simple array of zero depth to water to inform the soil water budget
# # without optimization
# dtw_zero = pd.DataFrame(np.zeros((len(dtw_avg), len(data_out))))
# # hopefully the parcel ids match as datatype (int vs str)
# dtw_zero.columns = data_out.parcel_id.values
# dtw_zero.index = dtw_avg.date

# for crop in crop_list:
#     var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
#     # need to account for when crops aren't predicted and skip them
#     if pred_dict[crop] in pred_crops: 
#         # to equalize the situation we might use a simple DTW profile

#         load_run_swb(crop, year, crop_in, join(loadpth, 'crop_soilbudget_no_irr'),
#                      dtw_zero, 
#                      soil_rep=False,
#                     run_opt=False) 

# %%
import functions.f_gw_dtw_extract as f_gw_dtw_extract
reload(f_gw_dtw_extract)
from functions.f_gw_dtw_extract import get_dtw

# %%
# if use == 'start':
    # if use is start then just sample initial groundwater level


# %%
# get head value from last 30 days to avoid using extreme single day value
spring_heads = avg_heads(spring_dates.kstpkper.values, hdobj, m)

# the dtw conversion runs a little slow
# get the DTW for the wels in the simulation from the last period
well_dtw = sample_dtw(spring_heads, botm)
# need to make integer for join with crop choice
well_dtw.UniqueID = well_dtw.UniqueID.astype(int)

# %% [markdown]
# The code that read DTW from the complete model run will no longer be needed with the updated input.

# %%
# # loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Regional/'
# dtw_model_ws = loadpth+'historical_simple_geology_reconnection'
# # dtw_model_ws = loadpth+'strhc1_scale'
# # year = 2020
# dtw_df_in = get_dtw(year, dtw_model_ws)
# dtw_df_in.to_csv(join(base_model_ws, 'field_SWB', 'dtw', 'dtw_ft_parcels_'+str(year)+'.csv'))

# # gw heads sampled from jan 1 to dec 31, heads from previous year can be easily
# # sampled from previous csv

# %%
# # dtw_df
# fp = join(base_model_ws, 'field_SWB', 'dtw','dtw_ft_parcels_'+str(year-1)+'.csv')
# if os.path.isfile(fp):
#     dtw_df_previous = pd.read_csv(fp,
#                                  index_col=0, parse_dates=['dt'])
#     dtw_df_previous.columns = dtw_df_previous.columns.astype(int)
#     dtw_df = pd.concat((dtw_df_previous.loc[str(year-1)+'-11-1':], dtw_df_in), axis=0)

# %% [markdown]
# # Irrigation submodel

# %% [markdown]
# The crop choice model didn't predict any corn.
# - this submodel may need to be in a separate script to run multiprocessing (parallel) as this re-loads the entire active script each time.

# %%
crop_list = ['Corn','Alfalfa','Pasture','Misc Grain and Hay', 'Grape']

# crop = 'Corn'
crop = crop_list[1]
# year = 2020
# load_run_swb(crop, year)

# %%
#
crop_in = data_out.merge(pod)
crop_in = crop_in.rename(columns={'Crop_Choice':'name'})
crop_in.to_csv(join(base_model_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'))
# crop_in[crop_in.name==crop]
pred_crops = crop_in.name.unique()
print(pred_crops)

# %%
crop_in.groupby('name')[['parcel_id']].count()

# %% [markdown]
# Originally we discussed the farmers assuming a standard increase (winter) and decline (summer) in groundwater based on the fall elevation. If we do this then we can run the optimization for 5-10 groundwater level patterns for each crop which will be more efficient than running for each unique sequence of groundwater levels.
# - we can validate this by running several example fields then translating back to each field with lookup table (use alfala which has only 7 selected)

# %%
# # for crop in ['Alfalfa']:
# for crop in ['Grape']:
#     var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
#     # need to account for when crops aren't predicted and skip them
#     if pred_dict[crop] in pred_crops: 
#         # load_run_swb(crop, year, crop_in, base_model_ws, dtw_df)
#         # reduce number for testing
#         load_run_swb(crop, year, crop_in, base_model_ws, dtw_df)

# %%
# dtw_df_crop_out = dtw_df[crop_in[crop_in.name==pred_dict[crop]].parcel_id.values]
# dtw_df_crop_out.to_csv(join(base_model_ws, 'field_SWB', 'dtw_ft_WY'+str(year)+'.csv'))

# %%
# dtw_crop_mean = dtw_df[data_out[data_out.Crop_Choice==pred_dict[crop]].parcel_id].loc['2020-4-1':].mean().values
# fig,ax=plt.subplots(figsize=(4,1))
# ax.plot(dtw_crop_mean)
# dtw_df.shape

# %% [markdown]
# ## Simplified representation of DTW range from min to max.
# Instead of using the full record of the DTW, the mdoel should just sample the average DTW for the month of interest

# %%
import functions.f_gw_dtw_extract
reload( functions.f_gw_dtw_extract)
from functions.f_gw_dtw_extract import calc_simple_dtw

# %%
## simple representative DTW for linear steps from dtw_min to dtw_max
## with a 5 ft decline from June to December based on observed data
dtw_simple_df = calc_simple_dtw(well_dtw, year, dtw_step = 10)

# %%
dtw_simple_df.to_csv(join(loadpth, 'rep_crop_soilbudget','field_SWB', 'dtw_ft_WY'+str(year)+'.csv'))

# %%
# pd.DataFrame(dtw_simple[:,0], dtw_avg.index).loc['2020-02-12':'2020-10-04'].plot()
# plt.plot(dtw_simple.mean(axis=0))
# pd.concat([dtw_simple_df]*2, axis=1)
# dtw_simple_df.iloc[:, ::2]

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
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
    # need to account for when crops aren't predicted and skip them
    # if pred_dict[crop] in pred_crops: 
    print(crop, ':',pred_dict[crop])

# %%
# initialize HDF5 files for the year
from functions.data_functions import init_h5
base_model_ws = join(loadpth, 'rep_crop_soilbudget')
# initialize SWB folder
os.makedirs(join(base_model_ws, 'field_SWB'), exist_ok=True)
for var in ['profit', 'yield', 'percolation','GW_applied_water', 'SW_applied_water']:
    name = join(base_model_ws, 'field_SWB', var + '_WY'+str(year)+'.hdf5')
    init_h5(name)

# check profit/yield if annual or daily

# %%


# for crop in ['Corn']:
for crop in crop_list:
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
    # need to account for when crops aren't predicted and skip them
    if pred_dict[crop] in pred_crops: 
        # to equalize the situation we might use a simple DTW profile
        load_run_swb(crop, year, crop_in, join(loadpth, 'rep_crop_soilbudget'),
                     dtw_simple_df, 
                     soil_rep=True) 

# %%
fn = join(model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
with h5py.File(fn) as dset:
    finished_crops = list(dset['array'].keys())
    print(finished_crops)
# only grape was completed?

# %% [markdown]
# 1. Load the representative results and sample for each field by crop type to back calculate the irrigation requirements. use the estimated irrigation as an input to the modflow model for pumping and percolation for recharge.
#     -   use the DTW id to reference to the irrigation in the full array, need to group by SW, GW or mixed.
#     -   if we wanted we could re-run the SWB one time with the specified irrigation rates to get the exact recharge rates with field specific values
# 2. RUn the modflow model to get the resultant DTW profile
# 3. re-calculate the profit using the irrigation and actual DTW profile on a soil by soil basis (non-optimization) after running the next modflow chunk. Actually the re-run for the true profit could be done if profits aren't needed mid-simulation

# %%
from functions import data_functions
reload(data_functions)
from functions.data_functions import read_crop_arr_h5

# %% [markdown]
# the creation of input files is set up to have to iterate over each crop individually. This means that the recharge and pumping for each crop needs to be appended to a main dataframe that is then used to create the modflow input.

# %%
model_ws = join(loadpth, 'rep_crop_soilbudget')

# %%
import functions.output_processing
reload(functions.output_processing)
from functions.output_processing import get_local_data, out_arr_to_long_df
from functions.output_processing import get_wb_by_parcel


# %%
# def get_wb_by_parcel(crop, crop_in, dtw_simple_df, well_dtw):

# %%
# initialize dataframes to save all the output from irrigation and recharge
pc_df_all = pd.DataFrame()
irr_sw_df_all = pd.DataFrame()
irr_gw_df_all = pd.DataFrame()

# %%
# for crop in crop_list:
# for crop in ['Grape']:
for crop in ['Alfalfa']:
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)

    # need separte hdf5 for each year because total is 300MB, group by crop in array
    fn = join(model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
    pc_all = read_crop_arr_h5(crop, fn)
    
    # # applied water (GW and SW are separate)
    fn = join(model_ws, 'field_SWB', "GW_applied_water_WY"+str(year)+".hdf5")
    irr_gw = read_crop_arr_h5(crop, fn)
    
    fn = join(model_ws, 'field_SWB', "SW_applied_water_WY"+str(year)+".hdf5")
    irr_sw = read_crop_arr_h5(crop, fn)

# %% [markdown]
# For alfalfa there are parcels with POD but the produced results are only for GW only?
# - the issue is that the fields all had POD so the groupby created only soil group, need to include flag or way to pass this through
#

# %%
# it would be useful to reference start/end dates for inputting the recharge/pumping

#  get the dates to update modflow
yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
# get the total extent of the irrigation season (calculation period)
strt_date = yield_start.min()
end_date = yield_end.max()
# strt_date, end_date

# %%
crop_ref = crop_in[crop_in.name==pred_dict[crop]]


# %%
dtw_df_mean = dtw_simple_df.mean().values
# if fields constains both with and without pod then double to account alternate run
if (crop_ref.pod_bool==1).any()&(crop_ref.pod_bool==0).any():
    dtw_df_mean = np.hstack([dtw_df_mean]*2)
    pod_bool = np.repeat([0, 1], int(len(irr_gw)/2))
elif (crop_ref.pod_bool==1).all():
    pod_bool = 1
else:
    pod_bool = 0
    # # summary output from hdf5 into csv

out_summary = pd.DataFrame(np.transpose((irr_gw.sum(axis=1), irr_sw.sum(axis=1), pc_all.sum(axis=1), dtw_df_mean)), 
             columns=['irr_gw_m','irr_sw_m','pc_m', 'dtw_mean_ft'])
out_summary['dtw_id'] = np.arange(0, len(out_summary))
# specify whether the representative profile is for with or without a POD
out_summary['pod_bool'] = pod_bool

# %%
# return the array that references DTW to the irr and pc for each represetnative
# soil water budget
out_lin = get_local_data(dtw_simple_df, well_dtw, crop_ref, irr_gw, irr_sw, pc_all)



# %%
# the dtw_id is the index for the representative outputs of irrigation and percolation
# need to translate these into modflow inputs
# calculate the total irrigation demand for each parcel for each day in the irrigation season
irr_tot = irr_gw[out_lin.dtw_id]+irr_sw[out_lin.dtw_id]
# sample the percolation rate for each parcel
rch_tot = pc_all[out_lin.dtw_id]

# %%

# create a dataframe with the irr or pc for each UniqueID parcel
pc_df_long = out_arr_to_long_df(pc_all, out_lin, strt_date, end_date)
irr_gw_df_long = out_arr_to_long_df(irr_gw, out_lin, strt_date, end_date)
irr_gw_df_long = out_arr_to_long_df(irr_sw, out_lin, strt_date, end_date)

# %% [markdown]
# It may make the most sense to have an intermediate dataframe with the output from each soil budget with dates to create a common reference.

# %%
import functions.output_processing
reload(functions.output_processing)
from functions.output_processing import get_local_data, out_arr_to_long_df
from functions.output_processing import get_wb_by_parcel

# %%
# load the processed dataframe with all datas
pc_df_all, irr_gw_df_all, irr_sw_df_all = get_wb_by_parcel(model_ws, crop, year, 
                 crop_in, finished_crops, dtw_simple_df, well_dtw)

# %%
# add crop type to array to check recharge rates
pc_check = pc_df_all.merge(crop_in, left_on='UniqueID',right_on='parcel_id')
# pc_check = pc_check[pc_check.name.str.contains('Vineyard')]
# pc_check.rate.max()
pc_check.groupby('name').mean(numeric_only=True)
# may want to do a little plotting to review distribution looks okay for 
# irrigation season recharge (should be low rates)

# %%
# join water budget data to field area for scaling to cell for recharge
rch_df = pc_df_all.merge(field_df)
# calculate the fraction of a cell covered by the field
rch_df['field_scale'] = rch_df.field_area/(200*200)
# calculate the recharge rate for each field
rch_df['rch_rate'] =  rch_df.rate*rch_df.field_scale
# set date as index for sampling
rch_df = rch_df.set_index('date')

# %% [markdown]
# Take the combined dataframe to inform the recharge array with recharge from the irrigation season. The next step is to insert recharge from non-irrigated times and native lands.

# %%

# join water budget data to field area for application to cell for well
wel_df = irr_gw_df_all.merge(well_loc_merge)
# calculate the volumetric flux of the well 
wel_df['flux'] = -wel_df.rate*wel_df.area_m2
# # set date as index for sampling
wel_df = wel_df.set_index('date')
wel_df

# %%
# allocate dictionary for recharge
rch_dict = dict()
wel_dict = dict()
# rch_df.set_index('date').loc[dates.min():dates.max()]
for t, d in enumerate(dates):
    # get data for the stress period
    rch_arr = np.zeros((m.dis.nrow, m.dis.ncol))
    if d in rch_df.index:
        spd_df = rch_df.loc[d]
        rch_arr[spd_df.row-1, spd_df.column-1] = spd_df.rch_rate
    # assign values to recharge dict
    rch_dict[t] = rch_arr
    if d in wel_df.index:
        spd_df = wel_df.loc[d]
        wel_dict[t] = spd_df[['k','i','j','flux']].values
# spd_df

# %%
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# plt.imshow(rch_dict[100])
# plt.colorbar(norm = mpl.colors.LogNorm)

# %% [markdown]
# # Run update modflow

# %%
for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date
    model_ws = loadpth + 'crop_modflow/'+str(m_strt.date())
    # The model should start in hydraulic connection
    if m_per > 0:
        model_ws_last = loadpth + 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date())
        hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')
        sp_last = hdobj.get_kstpkper()[-1]
        strt[:,:,:] = hdobj.get_data(sp_last)

    # for the first period we use dem as starting point
    # if solver criteria are not met, the model will continue if model percent error is less than stoperror
    bas_month = flopy.modflow.ModflowBas(model = m_month, ibound=ibound, strt = strt)
    # write bas input with updated starting heads
    bas_month.write_file()

# %%

# overwrite files that change
wel_spd = dict()
# rech_spd = dict()
# the first model period is October to April where we will assume pumping is 0
# for Ag
for n, t in enumerate(spd):
    # need to decide how to add in the pumping from domestic/municipal
    wel_spd[n] = np.array([[0,0,0,0]]) # m.wel.stress_period_data[t]
    # need to add in recharge from native lands
    # rech_spd[n] = m.rch.rech.array[t,0]

wel_month = flopy.modflow.ModflowWel(model=m_month, stress_period_data = wel_spd, ipakcb=55)
rch_month = flopy.modflow.ModflowRch(model=m_month, nrchop = 3, rech =  rech_spd, ipakcb=55)

wel_month.write_file()
rch_month.write_file()

#     success, buff = m_month.run_model()

# %% [markdown]
# # Re-run soil water budget to calculate actual profit
# The soil water budget needs to be re-run at the end of each year to calculate the profit based on the actual water availability for surface water and groundwater. This could be done at the very end of the simulation if the farmer doesn't need previous year profits between years.
