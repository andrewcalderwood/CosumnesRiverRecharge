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
# Script to to test the crop choice model against stata
#

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

import matplotlib.pyplot as plt
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

# %% [markdown]
# # load static data

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
# pod = data[['parcel_id','pod']].copy().rename(columns={'pod':'pod_bool'})
# pod['pod'] = 'No Point of Diversion on Parcel'
# pod.loc[pod.pod_bool==1, 'pod'] = 'Point of Diversion on Parcel'

# %%

crop_choice_dir = './'
# Read logit coefficients
logit_coefs = pd.read_csv(join(crop_choice_dir, 'logit_coefs.csv'))

# Import prior year revenue data by crop
rev_prior_yr_df = pd.read_csv(join(crop_choice_dir, "data_model/rev_prior_yr.csv"))


# %%

logit_coefs.columns.values[:-3]

# %%
data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
# data

# %% [markdown]
# # Predict crops

# %% [markdown]
# ##  with the old test data
#
# If running old script need to run with old test data as those should have stata results (crop_cat) based on before the crop_cat last was used.

# %%
data_out_all = pd.DataFrame()
for year in data.year.unique():
    stata_in = data[data.year==year]
    stata_in['wy_dc'] = wyt.loc[wyt.WY==year, 'wy_dc'].values[0]

    data_out = predict_crops(stata_in.reset_index().copy(), rev_prior_yr_df, logit_coefs, return_prob=False)
    data_out_all = pd.concat((data_out_all, data_out.assign(year=year)))

# %%
# columns relevant for comparison
stata_comp = data[['parcel_id','year','crop_cat']]

# %%
# merge stata input with predicted crop cateogry for quick check
# need to confirm if crop_cat is predicted by stata or observation
data_comp = data_out_all.merge(stata_comp)

# shows agreement but some difference

# %% [markdown]
# ## with new test data
#
# compare against stata results that account for crop_cat lats

# %%
out_dir = join(proj_dir, "model_results", 'stata_comp')

# %%
logit_coefs_adj = pd.read_csv(join(crop_choice_dir, 'logit_coefs_last_crop.csv'))
# logit_coefs_adj


# %%
crop_dict = logit_coefs_adj.set_index('Crop_Eq_new')['Crop_Eq'].to_dict()

# %%
# load results from stata to confirm python is doing as expected
stata_df = pd.read_csv(join(proj_dir, "model_results/stata_predicted_probabilities.csv"))

# replace crop name with original
stata_df = stata_df.replace(crop_dict)


# %%
# # when crop_cat == alternatives that is the same as chosen == 1
# stata_chk = stata_df[stata_df.chosen==1]
# stata_chk[stata_chk.crop_cat != stata_chk.alternatives]

# %%
crop_cat_last_df = stata_df.pivot(index=['parcel_id','year'], values='crop_cat_last', columns='alternatives')
# look to see if value in row is equal to column
cat_last_bool_arr = crop_cat_last_df.values==crop_cat_last_df.columns.values
# convert back to a dataframe for input
crop_bool_df = crop_cat_last_df.copy()
crop_bool_df.loc[:,:] = cat_last_bool_arr.astype(int)
# crop_bool_df

# %%
# for standard coefs use original input as reference
stata_use = stata_df[['parcel_id','year']+logit_coefs.columns[1:-3].tolist()]
stata_parcel = stata_use.groupby(['parcel_id','year']).first().reset_index()
# add boolean pod
stata_parcel = stata_parcel.rename(columns={'pod':'pod_str'})
stata_parcel['pod'] = 1
stata_parcel.loc[stata_parcel.pod_str=='No Point of Diversion on Parcel','pod'] = 0

# bring in the columns that identify the last crop type
stata_parcel = stata_parcel.merge(crop_bool_df, on=['parcel_id','year'])

# need to convert previous crop columns to int/float from object
stata_parcel[logit_coefs_adj.Crop_Eq.values] = stata_parcel[logit_coefs_adj.Crop_Eq.values].astype(float)


# %%
# name of crop_cat_last columns should be easy to match stata_df and to multiply for calculation of probability
logit_coefs_adj
# stata_df output is set up so crop_cat_last is defined for each row rather than with multiple columns
# crop_cat_last is a string but will need to be converted to a 1 or 0 for the multiplication
stata_parcel.columns

# %%
data_out_all = pd.DataFrame()
data_out_sel = pd.DataFrame()

for year in stata_parcel.year.unique():
# for year in [2018]:

    # has issue in np.exp
    stata_in = stata_parcel[stata_parcel.year==year]
    data_out = predict_crops_rand(stata_in.reset_index(drop=True).copy(), rev_prior_yr_df, logit_coefs_adj, return_prob=True)
    # old version runs without issue
    # data_out = predict_crops(stata_in.reset_index().copy(), rev_prior_yr_df, logit_coefs, return_prob=True)
    data_out_all = pd.concat((data_out_all, data_out.assign(year=year)))
    # run prediction again to get just the selected crop
    data_out = predict_crops_rand(stata_in.reset_index(drop=True).copy(), rev_prior_yr_df, logit_coefs_adj, return_prob=False)
    data_out_sel = pd.concat((data_out_sel, data_out.assign(year=year)))


# %%
# columns relevant for comparison
stata_comp = stata_df[['parcel_id','year','crop_cat_last','crop_cat','crop','alternatives', 'chosen', 'pprob_stata']]
# stata_comp = stata_df[['parcel_id','year','crop_cat','alternatives','pprob_stata']]
stata_comp.crop_cat = stata_comp.crop_cat.replace({'Mixed pasture and miscellaneous grasses':'Mixed_pasture'})


# %%
p_cols = data_out_all.columns[data_out_all.columns.str.contains('PP.')].tolist()
# also need to include which alternative is specified?
data_out_long = data_out_all.melt(id_vars=['parcel_id', 'year'], value_vars=p_cols, value_name='pprob_py')
data_out_long['alternatives'] = data_out_long.variable.str.replace('PP.','')

# names for crop are the same between old and new stata except for Mixed_pasture

# %%
data_out_long.alternatives.unique()

# %%
# join data on the selected crop from each for the parcel
stata_sel = stata_comp[stata_comp.chosen==1][['parcel_id','year','crop_cat']].rename(columns={'crop_cat':'crop_choice_stata'})
comp_sel = stata_sel.merge(data_out_sel.rename(columns={'Crop_Choice':'crop_choice_python'}))

# comp_sel[comp_sel.crop_choice_stata != comp_sel.crop_choice_python]
comp_sel.to_csv(join(out_dir, 'selected_crop_stata_python.csv'))

# %%
# alternative to looking at where it matches
# look where the same crop has the maximum probability

max_prob = stata_comp.groupby(['parcel_id','year'])['pprob_stata'].idxmax()
stata_max = stata_comp.loc[max_prob.dropna().values]

max_prob = data_out_long.groupby(['parcel_id','year'])['pprob_py'].idxmax()
py_max = data_out_long.loc[max_prob.dropna().values]
max_comp = stata_max.rename(columns={'crop_cat':'choice_stata'}).merge(py_max.rename(columns={'alternatives':'choice_py'}))

max_comp[max_comp.choice_stata != max_comp.choice_py]

# %%
prob_comp = stata_comp[['parcel_id','year','crop_cat','alternatives','pprob_stata']].merge(data_out_long)


# %%
prob_comp['error'] = prob_comp.pprob_stata - prob_comp.pprob_py
prob_comp['perc_error'] = prob_comp.error*100/((prob_comp.pprob_stata + prob_comp.pprob_py)/2)


prob_comp.to_csv(join(out_dir, 'stata_python_comparison.csv'))
# want to evaluate the error in probabilities individually to ensure not too significant

# more important test is to ensure that the randomly selected crop is the same

# %%
prob_comp.error.hist()
plt.xlabel('PP Error (Stata - Py)')
plt.ylabel('Count')
plt.savefig(join(out_dir, 'all_PP_error_histogram.png'), bbox_inches='tight')


# %%
# prob_comp

# %%
# look at parcels where the error is greater than 1% of the 100% distribution, below this it likely isn't
# making a big difference in random sampling or max selection
err_01 = prob_comp[prob_comp.error>0.01]

print('There are %.2f %% with greater than 1%% error' %(100*err_01.shape[0]/prob_comp.shape[0]))

err_01_count = err_01.groupby('alternatives')['error'].count()
err_01_count.plot(kind='bar')
plt.ylabel('Count of parcels across all years \nwith greater than 1% error')
plt.savefig(join(out_dir, 'count_1_perc_error_plot.png'), bbox_inches='tight')
# the error is pretty evenly distributed among the crop type alternatives with the lowest among other



# %%
sns.catplot(
    prob_comp,
    # err_01, 
            x='year',y='error', col='alternatives', col_wrap=3, 
            kind='box',#color='tab:blue',
            sharey=False)

plt.savefig(join(out_dir, 'frac_error_box_plot_by_crop_and_year.png'), bbox_inches='tight')


# %%
# compare the selected crop in each case


# %% [markdown]
# # Run the crop choice in python
#
# Comparing with simulated dtw

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

