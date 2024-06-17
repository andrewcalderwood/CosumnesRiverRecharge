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

# %%
from parcelchoicemodelupdate.f_predict_landuse import predict_crops

# %%
import Basic_soil_budget_monthly as swb


# %% [markdown]
# If we are going to run the crop/swb for a certain period then there should already have been a MODFLOW model run for that same period with estimates of streamflow and precipitation to drive the other boundary conditions.

# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Regional/'

# base_model_ws = loadpth + 'crop_soilbudget'
m_model_ws = loadpth + 'historical_simple_geology_reconnection'

# %%
load_only=['DIS']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= m_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only)

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
    
all_run_dates = pd.concat((pd.DataFrame([all_strt_date]).assign(use='start'), all_run_dates))
all_run_dates = pd.concat((pd.DataFrame([all_end_date]).assign(use='end'), all_run_dates))
all_run_dates=all_run_dates.sort_values(0).reset_index(drop=True).rename(columns={0:'date'})

# %% [markdown]
# Review season dates, the plan was to change the start of the irrigation date for misc. grain and hay since irrigators don't typically start until summer even though it grows in winter. Simplify to crop choice Jan 1 and Irrig. run Apr 1
# - also at some point we may have discussed just doing them at the same time? in that case we just need to iterate over years

# %%
# all_run_dates
# simple code to set dates for april 1
all_run_dates = pd.date_range(all_strt_date, all_end_date,freq='AS-Apr')

# %%
# this loop was set to run for the years of interest
for m_per in np.arange(1,5): # runs first year to next crop choice
# for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date
    m_end = all_run_dates.iloc[m_per+1].date

    dates = pd.date_range(m_strt, m_end)

# %%
year = 2020
crop='Corn'

# %%
# the model will run the irrigation optimizer on specified dates (multiple crops can be done at once or in sequence)
# the modflow model will be run for the periods between the specified irrigation optimizer dates
# loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Regional/'
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Regional/'

base_model_ws = loadpth + 'crop_soilbudget'
os.makedirs(base_model_ws, exist_ok=True)

# %% [markdown]
# # Crop choice model
# Evetnaully this should use the updated DTW from each previous year.

# %%
crop_choice_dir = 'parcelchoicemodelupdate'
# Read logit coefficients
logit_coefs = pd.read_csv(join(crop_choice_dir, 'logit_coefs.csv'))

# the parcel data needs the dtwfa (avg dtw in feet for the parcel) and wy_dc (pulled from Sac wy type dataset and switched to dry boolean)
# missing WY type prediction? 
# Read parcel data
data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
data['wy_dc'] = np.where(data['year'] == 2020, 1, 0) # should be pulled from Sac WY type data instead

# Import prior year revenue data by crop
rev_prior_yr_df = pd.read_csv(join(crop_choice_dir, "data_model/rev_prior_yr.csv"))

# %%
pod = data[['parcel_id','pod']].copy().rename(columns={'pod':'pod_bool'})
pod['pod'] = 'No Point of Diversion on Parcel'
pod.loc[pod.pod_bool==1, 'pod'] = 'Point of Diversion on Parcel'

# %%
# crop choice model uses "_" instead of " "
# the irrigation model set up uses " "
data_out = predict_crops(data.copy(), rev_prior_yr_df, logit_coefs)
data_out['Crop_Choice'] = data_out.Crop_Choice.str.replace('_',' ')
data_out.to_csv(join(crop_choice_dir, 'parcel_crop_choice_'+str(year)+'.csv'), index=False)

# %% [markdown]
# # Load MF output
# The depth to water function should sample from the previous model run which may be a year or less long.
#
# The original get_dtw function was set up assuming a continuous model run which won't be the case.

# %%
import functions.f_gw_dtw_extract as f_gw_dtw_extract
reload(f_gw_dtw_extract)
from functions.f_gw_dtw_extract import get_dtw

# %%
# loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Regional/'
dtw_model_ws = loadpth+'historical_simple_geology_reconnection'
# dtw_model_ws = loadpth+'strhc1_scale'
year = 2020
dtw_df_in = get_dtw(year, dtw_model_ws)
dtw_df_in.to_csv(join(base_model_ws, 'field_SWB', 'dtw_ft_parcels_'+str(year)+'.csv'))

# gw heads sampled from jan 1 to dec 31, heads from previous year can be easily
# sampled from previous csv

# %%
# dtw_df
fp = join(base_model_ws, 'field_SWB', 'dtw','dtw_ft_parcels_'+str(year-1)+'.csv')
if os.path.isfile(fp):
    dtw_df_previous = pd.read_csv(fp,
                                 index_col=0, parse_dates=['dt'])
    dtw_df_previous.columns = dtw_df_previous.columns.astype(int)
    dtw_df = pd.concat((dtw_df_previous.loc[str(year-1)+'-11-1':], dtw_df_in), axis=0)

# %%
import matplotlib.pyplot as plt

# %%
# plt.plot((dtw
# _df_previous.loc[str(year-1)+'-11-1':].values - dtw_df_in.loc[str(year)+'-11-1':].values));

# %%
# example depth to water simulation results to show variability in trends
# this seems to show the model has some extreme locations of drawdown that are unrealistic
# dtw_df.iloc[:,::205].plot(legend=False);
# dtw_df_in.iloc[:,::205].plot(legend=False);
# plt.ylabel('Depth to water (ft)')

# %%
# # the plot makes it seem like there is a sharp discontinuity here
# # maybe its the same values as 2020 on accident, maybe issue with hdobj read in
# # dtw_df.iloc[:,0].plot()
# dtw_df_in.iloc[:,0].plot()
# dtw_df_in.iloc[:,200].plot()
# dtw_df_in.iloc[:,600].plot()
# # dtw_df_in.iloc[:,800].plot()
# dtw_df_in.iloc[:,1200].plot()


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

# %%
# subset for parcels for the current crop


# %%
import Basic_soil_budget_monthly as swb
reload(swb)
# import f_swb_profit_opt
# reload(f_swb_profit_opt)
# from f_swb_profit_opt import load_run_swb

import f_rep_swb_profit_opt
reload(f_rep_swb_profit_opt)
from f_rep_swb_profit_opt import load_run_swb

# %% [markdown]
# Alfalfa ran as expected <1 min each. Misc grain and hay is running very slow, some taking up to 20 min with the norm being 3-5 minutes it seems (negative profits confuse function?)
#
# - the insane run time makes it seem more reasonable now to simply create a linear relationship between irrigation rate and DTW for each field and water year.
#     - lets do DTW with 5 ft steps from 0 to 300 ft for 60 iterations
#     - irrigation dates range from 20-80 per season and within that we need accuracy to the 1 cm and if we have a max irrigation of 150 cm (almonds take ~40 in/yr and 60 inches is about 150 cm) thats 150 steps which would create 150^80 permutations (not accounting for duplicates) or 150!/70! which is equally large and that's before accounting for hydrology so this wouldn't be feasible.
# - but also it takes almost a minute per soil type for vineyards, for some reason the first one is the fastest at 0.15 min, this means it would take 1,000 min or about 16 hrs which is not reasonable even in parallel

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

# %%
## simple representative DTW for linear steps from dtw_min to dtw_max
## with a 5 ft decline from June to December based on observed data

# get min and max dtw to reference
dtw_step = 10
dtw_min, dtw_max = dtw_df.mean(axis=0).quantile([0,1]).apply(lambda x: base_round(x, dtw_step))

# calculate dates for DTW decline, shouldn't really go longer than a year since it 
# will be run at least once per year
dtw_avg = pd.DataFrame(pd.date_range(str(year-1)+'-11-1', str(year)+'-12-31'), columns=['date'])
dtw_avg = dtw_avg.assign(decline = 0).set_index('date')
# dates where a decline date is specified
decline_dates = dtw_avg.index[dtw_avg.index >=str(year)+'-6-1']
decline_total = 5
decline = np.cumsum(np.full(len(decline_dates), decline_total/len(decline_dates)))
dtw_avg.loc[decline_dates, 'decline'] = decline
dtw_simple = np.repeat(np.reshape(np.arange(dtw_min, dtw_max, dtw_step), (1,-1)), len(dtw_avg), axis=0)
dtw_simple = dtw_simple + np.reshape(dtw_avg.decline, (-1,1))
dtw_simple_df = pd.DataFrame(dtw_simple, dtw_avg.index)
# plt.plot(dtw_simple[:,0])

# %%
dtw_simple_df.to_csv(join(loadpth, 'rep_crop_soilbudget','field_SWB', 'dtw_ft_WY'+str(year)+'.csv'))

# %%
# pd.DataFrame(dtw_simple[:,0], dtw_avg.index).loc['2020-02-12':'2020-10-04'].plot()
# plt.plot(dtw_simple.mean(axis=0))
# pd.concat([dtw_simple_df]*2, axis=1)
# dtw_simple_df.iloc[:, ::2]

# %% [markdown]
# ## Iterate over all crops to save the representative profiles

# %%
for crop in crop_list:
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
    # need to account for when crops aren't predicted and skip them
    if pred_dict[crop] in pred_crops: 
        # to equalize the situation we might use a simple DTW profile
        load_run_swb(crop, year, crop_in, join(loadpth, 'rep_crop_soilbudget'),
                     dtw_simple_df, 
                     soil_rep=True) 

# %% [markdown]
# 1. Load the representative results and sample for each field by crop type to back calculate the irrigation requirements. use the estimated irrigation as an input to the modflow model for pumping and percolation for recharge.
#     -   use the DTW id to reference to the irrigation in the full array, need to group by SW, GW or mixed.
# 2. RUn the modflow model to get the resultant DTW profile
# 3. re-calculate the profit using the irrigation and actual DTW profile on a soil by soil basis (non-optimization) after running the next modflow chunk. Actually the re-run for the true profit could be done if profits aren't needed mid-simulation

# %%
from functions import data_functions
reload(data_functions)
from functions.data_functions import read_crop_arr_h5

# %%
model_ws = join(loadpth, 'rep_crop_soilbudget')

# %%
# for crop in crop_list:
for crop in ['Grape']:
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)

    # need separte hdf5 for each year because total is 300MB, group by crop in array
    fn = join(model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
    pc_all = read_crop_arr_h5(crop, fn)
    
    # # applied water (GW and SW are separate)
    fn = join(model_ws, 'field_SWB', "GW_applied_water_WY"+str(year)+".hdf5")
    irr_gw = read_crop_arr_h5(crop, fn)
    
    fn = join(model_ws, 'field_SWB', "SW_applied_water_WY"+str(year)+".hdf5")
    irr_sw = read_crop_arr_h5(crop, fn)

# %%
# it would be useful to reference start/end dates for inputting the recharge/pumping
season

# %%
crop_ref = crop_in[crop_in.name==pred_dict[crop]]


# %%
dtw_df_mean = dtw_simple_df.mean().values
# if list of fields contains a field with pod then double to account alternate run
if (crop_ref.pod_bool==1).any():
    dtw_df_mean = np.hstack([dtw_df_mean]*2)
# # summary output from hdf5 into csv
out_summary = pd.DataFrame(np.transpose((irr_gw.sum(axis=1), irr_sw.sum(axis=1), pc_all.sum(axis=1), dtw_df_mean)), 
             columns=['irr_gw_m','irr_sw_m','pc_m', 'dtw_mean_ft'])
# # actually need to reference the array data directly rather than the sum
# the nearest value merge should identify the index which should then be used to reference the irrigation/percolation
# arrays
out_summary['dtw_id'] = np.arange(0, len(out_summary))
# specify whether the representative profile is for with or without a POD
out_summary['pod_bool'] = np.repeat([0, 1], int(len(irr_gw)/2))


# %%
# take the mean dtw for each field
dtw_id_mean = dtw_df.mean()
# sample the dtw for the field simulated for that crop
dtw_id_mean = pd.DataFrame(dtw_id_mean.loc[crop_ref.parcel_id.values], columns=['dtw_mean_ft'])
# want to join dtw info with pod info
dtw_id_mean = crop_ref.join(dtw_id_mean, on='parcel_id')
# # add ID to identify which parcel is selected in a simpler id than parcel_id
dtw_id_mean['id'] = np.arange(0, len(dtw_id_mean))


# %%
out_summary_dtw = out_summary[['dtw_id','dtw_mean_ft', 'pod_bool']].sort_values('dtw_mean_ft')
# identifies the nearest dtw value, should correct with a linear interpolation
# using the slope in irr/recharge at each point scaled by difference in 
# need to do separately for pod/no pod
out_lin0 = pd.merge_asof( 
    dtw_id_mean[dtw_id_mean.pod_bool==0].sort_values('dtw_mean_ft'),
    out_summary_dtw[out_summary_dtw.pod_bool==0].drop(columns='pod_bool'),
                        on='dtw_mean_ft', direction='nearest')
out_lin1 = pd.merge_asof(
    dtw_id_mean[dtw_id_mean.pod_bool==1].sort_values('dtw_mean_ft'), 
    out_summary_dtw[out_summary_dtw.pod_bool==1].drop(columns='pod_bool'),
                        on='dtw_mean_ft', direction='nearest')
# rejoin complete table of pod and no pod
out_lin = pd.concat((out_lin0, out_lin1))
# sort values for plotting and add back the interpoalted DTW
out_lin = out_lin.rename(columns={'dtw_mean_ft':'dtw_mean_ft_sim'}).merge( out_summary_dtw)

# # should scale the irrigation in some way based on the linear interpolation of DTW

# %%
# the dtw_id is the index for the representative outputs of irrigation and percolation
# need to translate these into modflow inputs
# calculate the total irrigation demand for each parcel for each day in the irrigation season
irr_tot = irr_gw[out_lin.dtw_id]+irr_sw[out_lin.dtw_id]
# sample the percolation rate for each parcel
pc_all[out_lin.dtw_id]


# %%
#  get the dates to update modflow
yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
# get the total extent of the irrigation season (calculation period)
strt_date = yield_start.min()
end_date = yield_end.max()
strt_date, end_date

# %% [markdown]
# # Run update modflow

# %%

# %% [markdown]
# # Re-run soil water budget to calculate actual profit
# The soil water budget needs to be re-run at the end of each year to calculate the profit based on the actual water availability for surface water and groundwater. This could be done at the very end of the simulation if the farmer doesn't need previous year profits between years.
