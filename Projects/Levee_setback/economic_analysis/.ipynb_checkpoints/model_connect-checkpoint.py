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
# year = 2020
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
# the plot makes it seem like there is a sharp discontinuity here
# dtw_df.iloc[:,0].plot()

# %% [markdown]
# # Irrigation submodel

# %% [markdown]
# The crop choice model didn't predict any corn.

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

# %%
# subset for parcels for the current crop


# %%
import Basic_soil_budget_monthly as swb

import f_swb_profit_opt
reload(f_swb_profit_opt)

from f_swb_profit_opt import load_run_swb

# %% [markdown]
# Alfalfa ran as expected <1 min each. Misc grain and hay is running very slow, some taking up to 20 min with the norm being 3-5 minutes it seems (negative profits confuse function?)
#
# - the insane run time makes it seem much more reasonable now to simply create a linear relationship between irrigation rate and DTW for each field and water year.
#     - lets do DTW with 5 ft steps from 0 to 300 ft for 60 iterations
#     - irrigation dates range from 20-80 per season and within that we need accuracy to the 1 cm and if we have a max irrigation of 150 cm (almonds take ~40 in/yr and 60 inches is about 150 cm) thats 150 steps which would create 150^80 permutations (not accounting for duplicates) or 150!/70! which is equally large

# %%
for crop in crop_list:
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
    # need to account for when crops aren't predicted and skip them
    if pred_dict[crop] in pred_crops: 
        load_run_swb(crop, year, crop_in, base_model_ws, dtw_df)

# %% [markdown]
# # Run update modflow

# %%
