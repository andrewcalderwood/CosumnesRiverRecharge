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
from parcelchoicemodelupdate.f_predict_landuse import predict_crops

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

# %%
import f_gw_dtw_extract
reload(f_gw_dtw_extract)
from f_gw_dtw_extract import get_dtw

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
fp = join(base_model_ws, 'field_SWB', 'dtw_ft_parcels_'+str(year-1)+'.csv')
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
# the plot makes it seem like there is a sharp discontinuity here
# maybe its the same values as 2020 on accident, maybe issue with hdobj read in
# dtw_df.iloc[:,0].plot()
dtw_df_in.iloc[:,0].plot()
dtw_df_in.iloc[:,200].plot()
dtw_df_in.iloc[:,600].plot()
# dtw_df_in.iloc[:,800].plot()
dtw_df_in.iloc[:,1200].plot()


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
# crop_in.groupby('name')[['parcel_id']].count()

# %%
# subset for parcels for the current crop


# %%
import Basic_soil_budget_monthly as swb

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
for crop in ['Alfalfa']:
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
    # need to account for when crops aren't predicted and skip them
    if pred_dict[crop] in pred_crops: 
        load_run_swb(crop, year, crop_in, base_model_ws, dtw_df)

# %%
dtw_df.loc[:, crop_in[crop_in.name==pred_dict[crop]].parcel_id].mean()

# %%

# %% [markdown]
# Simplified representation of DTW

# %%
## simple representative DTW for linear steps 10 ft to 200 ft
## with a 5 ft decline from June to December based on observed data
dtw_avg = pd.DataFrame(pd.date_range(str(year-1)+'-11-1', str(year)+'-12-31'), columns=['date'])
dtw_avg = dtw_avg.assign(decline = 0).set_index('date')
# dates where a decline date is specified
decline_dates = dtw_avg.index[dtw_avg.index >=str(year)+'-6-1']
decline_total = 5
decline = np.cumsum(np.full(len(decline_dates), decline_total/len(decline_dates)))
dtw_avg.loc[decline_dates, 'decline'] = decline
dtw_simple = np.repeat(np.reshape(np.arange(10, 200, 10), (1,-1)), len(dtw_avg), axis=0)
dtw_simple = dtw_simple + np.reshape(dtw_avg.decline, (-1,1))
# plt.plot(dtw_simple[:,0])

# %%
# pd.DataFrame(dtw_simple[:,0], dtw_avg.index)

# %%
# to equalize the situation we might use a simple DTW profile
load_run_swb(crop, year, crop_in, join(loadpth, 'rep_crop_soilbudget'),
             pd.DataFrame(dtw_simple, dtw_avg.index), soil_rep=True)

# %% [markdown]
# Testing by crop

# %%
# for crop in crop_list[:2]:
for crop in ['Misc Grain and Hay']:
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
    # need to account for when crops aren't predicted and skip them
    if pred_dict[crop] in pred_crops: 
        load_run_swb(crop, year, crop_in, base_model_ws, dtw_df)

# %%
# # it takes a while to load the matlab engine but onces it's loaded then it can call functions
# import matlab.engine
# eng = matlab.engine.start_matlab()
# tf = eng.isprime(37)
# print(tf)

# %%

# %% [markdown]
#

# %% [markdown]
# # Run update modflow

# %%

# %% [markdown]
# # Re-run soil water budget to calculate actual profit
# The soil water budget needs to be re-run at the end of each year to calculate the profit based on the actual water availability for surface water and groundwater. This could be done at the very end of the simulation if the farmer doesn't need previous year profits between years.
