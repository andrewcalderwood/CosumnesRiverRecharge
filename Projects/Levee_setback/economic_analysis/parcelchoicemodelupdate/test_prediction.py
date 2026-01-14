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
# 4. Start the next year
#
# And run different management and/or recharge scenarios.

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

import seaborn as sns

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
import f_predict_landuse
reload(f_predict_landuse)
from f_predict_landuse import predict_crops, predict_crops_rand

# %%
# load Sac Valley WYT
wyt = pd.read_csv(join(proj_dir, 'model_inputs', 'sacramento_WY_types.txt'))
# whether it is a critical or dry year or not
wyt['wy_dc'] = 0
# # for wet or above normal years the boolean will be 0 
wyt.loc[wyt['Yr-type'].isin(['C','D']),'wy_dc'] = 1

# %%
crop_choice_dir = './'
data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
# expect updated column name for pod
# add column to POD that was available in previous dataset
pod = data[['parcel_id','pod']].copy().rename(columns={'pod':'pod_bool'})
pod['pod'] = 'No Point of Diversion on Parcel'
pod.loc[pod.pod_bool==1, 'pod'] = 'Point of Diversion on Parcel'

# %%

crop_choice_dir = './'
# Read logit coefficients
logit_coefs = pd.read_csv(join(crop_choice_dir, 'logit_coefs.csv'))

# Import prior year revenue data by crop
rev_prior_yr_df = pd.read_csv(join(crop_choice_dir, "data_model/rev_prior_yr.csv"))


# %%
dtw_ft = pd.read_csv("D:/WRDAPP/GWFlowModel/Cosumnes/Economic/input_write_2014_2022_R203/output_clean/dtw_ft_mean_all.csv")


# %% [markdown]
# # compare python max prob vs rand sample methods

# %%
data_out_all = pd.DataFrame()
data_out2_all = pd.DataFrame()
for year in np.arange(2016,2022):
# year = 2017
    # the parcel data needs the dtwfa (avg fall dtw in feet for the parcel) and wy_dc (pulled from Sac wy type dataset and switched to dry boolean)
    # missing WY type prediction? 
    # Read parcel data
    data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
    # still needs to be updated to auto update DTW and WY type
    # data['wy_dc'] = np.where(data['year'] == 2020, 1, 0) # should be pulled from Sac WY type data instead
    data['wy_dc'] = wyt.loc[wyt.WY==year, 'wy_dc'].values[0]
    # update DTW to use simulated value instead of contoured
    well_dtw = dtw_ft[dtw_ft.year==year]
    well_dtw_merge = well_dtw[['UniqueID','dtw_ft']].rename(columns={'dtw_ft':'dtwfa'})
    data = data.drop(columns=['dtwfa', 'dtwsp'])
    data = data.merge(well_dtw_merge, left_on='parcel_id',right_on='UniqueID')
    
    
    data_out = predict_crops(data.copy(), rev_prior_yr_df, logit_coefs)
    data_out2 = predict_crops_rand(data.copy(), rev_prior_yr_df, logit_coefs)
    
    data_out_all = pd.concat((data_out_all, data_out.assign(year=year)))
    data_out2_all = pd.concat((data_out2_all, data_out2.assign(year=year)))

# %%
out_comp = pd.concat((data_out_all.assign(grp = 'max_p'), data_out2_all.assign(grp='rand')))

# %%
comp_count = out_comp.groupby(['Crop_Choice','year','grp'])['parcel_id'].count().reset_index()

# %%
# sns.relplot(data_out_all
sns.catplot(comp_count,x='year',y='parcel_id', col='Crop_Choice', col_wrap=3, hue='grp',
            kind='bar', #color='tab:blue',
            sharey=False
           # facet_kws={'sharey': False, 'sharex': True}
)

# sns.histplo


    # %%
    # crop choice model uses "_" instead of " "
    # the irrigation model set up uses " "

    # data_out['Crop_Choice'] = data_out.Crop_Choice.str.replace('_',' ')
    # # update naming of Corn
    # data_out.Crop_Choice = data_out.Crop_Choice.str.replace('Corn  ','Corn, ')
    # # save output with only parcel and crop choice
    # data_out.to_csv(join(swb_ws, 'field_SWB', 'parcel_crop_choice_'+str(year)+'.csv'))


# %% [markdown]
# # Compare new stata to python

# %%
stata_new = pd.read_csv(join(proj_dir, "model_results/stata_predicted_probabilities.csv"))

# %%
stata_new
