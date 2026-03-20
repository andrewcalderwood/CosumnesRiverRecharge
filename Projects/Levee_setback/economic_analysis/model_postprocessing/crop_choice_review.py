# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Script to to check the crop choice model results against observations based on MODFLOW predictions
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
crop_choice_dir = '../parcelchoicemodelupdate'
add_path('../functions')
add_path(crop_choice_dir)

# import f_predict_landuse
# reload(f_predict_landuse)
from f_predict_landuse import predict_crops, predict_crops_rand

# %% [markdown]
# # load static data

# %%
fig_dir = join(proj_dir, "model_results","stata_comp", 'figures')
os.makedirs(fig_dir, exist_ok=True)

# %%
run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Economic/'

# %%
# load Sac Valley WYT
wyt = pd.read_csv(join(proj_dir, 'model_inputs', 'sacramento_WY_types.txt'))
# whether it is a critical or dry year or not
wyt['wy_dc'] = 0
# # for wet or above normal years the boolean will be 0 
wyt.loc[wyt['Yr-type'].isin(['C','D']),'wy_dc'] = 1

# %%
wyt[wyt.WY>2014][['WY','Yr-type','wy_dc']].to_csv(join(proj_dir, 'model_results','stata_comp','WYT_2015_2022.csv'))

# %%
data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
# expect updated column name for pod
# add column to POD that was available in previous dataset
# pod = data[['parcel_id','pod']].copy().rename(columns={'pod':'pod_bool'})
# pod['pod'] = 'No Point of Diversion on Parcel'
# pod.loc[pod.pod_bool==1, 'pod'] = 'Point of Diversion on Parcel'

# %%

# Read logit coefficients
logit_coefs = pd.read_csv(join(crop_choice_dir, 'logit_coefs.csv'))

# Import prior year revenue data by crop
rev_prior_yr_df = pd.read_csv(join(crop_choice_dir, "data_model/rev_prior_yr.csv"))

# this is supposed to vary by year but right now is static


# %%
# old input data and test
data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
# data

# %%
logit_coefs_adj = pd.read_csv(join(crop_choice_dir, 'logit_coefs_last_crop.csv'))
# logit_coefs_adj


# %%
crop_dict = logit_coefs_adj.set_index('Crop_Eq_new')['Crop_Eq'].to_dict()

crop_dict['Mixed pasture'] = 'Mixed_pasture'

# %%
# load results from stata to confirm python is doing as expected
stata_df = pd.read_csv(join(proj_dir, "model_results/stata_predicted_probabilities.csv"))

# replace crop name with original
stata_df = stata_df.replace(crop_dict)


# %%
# save revenue as a dataframe for reference
rev_out = stata_df[stata_df.chosen==1].drop_duplicates(['alternatives','year'])[['alternatives','year','rev_prior_yr']]
rev_out = rev_out.sort_values(['alternatives','year']).reset_index(drop=True).rename(columns={'alternatives':'crop_cat'})
rev_out = rev_out.replace({'Mixed pasture and miscellaneous grasses':'Mixed_pasture'})
rev_out.to_csv(join(crop_choice_dir, "data_model", 'rev_prior_yr_all.csv'), index=False)
# this is still missing 2016 information

# %%
# identify crops specified in LandIQ
stata_chosen = stata_df[stata_df.chosen==1]
# identify the full list of parcels
stata_chosen_init = stata_chosen.sort_values('year').drop_duplicates(['parcel_id']).reset_index()
stata_chosen_init.to_csv(join(proj_dir, 'model_inputs', 'stata_input_initial.csv'))

# %% [markdown]
# Review of crops that should fill in other

# %%
uzf_dir = join(gwfm_dir,'UZF_data')

ag_lu = gpd.read_file(join(uzf_dir,'county_landuse', 'domain_ag_lu_2018.shp'))
# 5 duplicates in irrigation efficiency
ag_irr_eff = ag_lu[['geom_id', 'name', 'irr_name', 'Avg_eff']].drop_duplicates()
# loaded to evaluate if worth using to specify Other, determined best to assign average conditions of most crops

# %%
ag_irr_eff.name.unique()
opt_crops = ['Mixed pasture','Grain and hay crops','Miscellaneous grain and hay', 
             'Corn (field & sweet)','Sudan','Vineyards','Alfalfa & alfalfa mixtures']
top_count = ag_irr_eff[~ag_irr_eff.name.isin(opt_crops)].groupby('name')['Avg_eff'].count().sort_values(ascending=False)[:12]
# adding up the top 12 categories is 395 fields
# the major types are melons, strawberries, rice, tomatoes, clover, walnuts (in grape category), flow nuresry, safflow, misc truck

# %%
# load in pre-processed array of ETc for all time, m/day
ETc_long = pd.read_hdf(join(uzf_dir, "dwr_ETc",'long_ETc_all_lu.hdf5'), key='variable')
# convert to column format
ETc_all = ETc_long.pivot(columns='variable', values='value')


# %%
# test to review top categories to see which are consistent
# rep_ETc = ETc_all[top_count.index[top_count.index.isin(ETc_all.columns)]]
# rep_ETc[rep_ETc.index>'2023-1-1'].plot()
# remove Mixed deciduous and safflower but the others all follow a fairly consistent trend/magnitude
other_et_cat = ['Strawberries','Rice','Tomatoes (processing)',
                'Melons, squahs and cucumbers (all types',
                'Miscellaneous truck','Flowers nursery & Christmas tree farms']
# use the subset categories to represent other
rep_ETc = ETc_all[other_et_cat]
# use average ET to represent
rep_ETc['Other'] = rep_ETc.mean(axis=1)

# %% [markdown]
# ## depth to water review
#  Give DTW in fall wasn't explicitly stated, we need to sample it or extract it again from a dataset

# %%

model_nam = 'input_write_2014_2022'
scen = 'R200'
# it is probably better to create a slightly different file name then to copy these over for a set scenario
econ_model_ws = join(loadpth, model_nam+'_'+scen)

all_run_dates = pd.read_csv(join(econ_model_ws, 'crop_modflow', 'all_run_dates.csv'),parse_dates=['date'])

years = all_run_dates[all_run_dates.use=='irrigation'].date.dt.year

swb_ws = join(econ_model_ws, 'rep_crop_soilbudget')


# %%
# crop_soilbudget/field_dtw/dtw_ft+crop data has dtw for re-running SWB to get actual costs
# issue is that it uses a ffill() for the last month before next model is run and it stops at end of irrigation season

# %%
from f_gw_dtw_extract import sample_dtw, avg_heads


# %%
m_model_ws = join(dirname(dirname(loadpth)), 'Regional', model_nam)
load_only=['DIS', 'BAS6']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= m_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only)
m_dim = (m.dis.nlay, m.dis.nrow,  m.dis.ncol)


# %%
parcel_wells = pd.read_csv(join(econ_model_ws, 'crop_modflow', 'parcel_wells_with_layer.csv'))

# %%
well_dtw_all = pd.DataFrame()
for m_per in np.arange(1, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date
    year = m_strt.year

    model_ws_last = join(econ_model_ws, 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date()))
    hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')

    # model output dates
    m_dates = all_run_dates.loc[m_per-1].date+np.array(hdobj.get_times()).astype('timedelta64[D]')
    m_dates = pd.DataFrame(m_dates, columns=['dates']).set_index('dates')
    m_dates['kstpkper'] = hdobj.get_kstpkper()
    # subset to 1 month of output
    # determine dates for fall sampling
    fall_dates = m_dates[m_dates.index.month==10]
    # determine dates for spring sampling
    spring_dates = m_dates[m_dates.index.month==3]

    # get head value from last 30 days to avoid using extreme single day value
    fall_heads = avg_heads(fall_dates.kstpkper.values, hdobj, m_dim)
    
    # the dtw conversion runs a little slow
    # get the DTW for the wels in the simulation from the last period
    well_dtw = sample_dtw(fall_heads, parcel_wells)
    # need to make integer for join with crop choice
    well_dtw.UniqueID = well_dtw.UniqueID.astype(int)
    well_dtw_all = pd.concat((well_dtw_all, well_dtw.assign(year=year)))

# %% [markdown]
# ## comparison of DTW
#
# The histogram of stata DTW shows a range of 10-40 with median around 20. This seems to suggest to me the units are in meters and not feet.

# %%
# gpd.read_file(join(proj_dir, 'parcel_zonalstats', 'gw_dtw_all.shp'))
# verifies that dtw was in meters in dataset sent to Yusuke

# %%
# get just DTW info used in stata
stata_dtw = stata_df[['parcel_id','year','dtwfa']]
# compare results in Stata with those in MODFLOW
dtw_comp = stata_dtw.merge(well_dtw_all.rename(columns={'UniqueID':'parcel_id'}))
# if the stata dtw is converted to ft then it aligns with the model predicted
dtw_comp.dtwfa /= 0.3048
# seems highly likely that stata was used dtw in meters, which is okay because it was consistent
# that just means that the use of stata in the model needs to do the same

# %%
fig,ax = plt.subplots(2,1)
dtw_comp.dtwfa.hist(ax=ax[0])
ax[0].set_title('STATA')

dtw_comp.dtw_ft.hist( ax=ax[1])
ax[1].set_title('MODFLOW')
ax[1].set_xlabel('DTW (ft)')
plt.savefig(join(fig_dir, 'dtw_histogram_stata_modflow.png'), bbox_inches='tight')

# %% [markdown]
# # Predict crops

# %% [markdown]
# ## with new test data
#
# compare against stata results that account for crop_cat lats

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

# %% [markdown]
# ## prediction based on observation data
#
# Could also re-run prediction with new dtw and previous crop to see how different

# %%
# data_out_all = pd.DataFrame()
# data_out_sel = pd.DataFrame()
# # 
# for year in stata_parcel.year.unique():
# # for year in [2018]:

#     # has issue in np.exp
#     stata_in = stata_parcel[stata_parcel.year==year]
#     data_out = predict_crops_rand(stata_in.reset_index(drop=True).copy(), rev_prior_yr_df, logit_coefs_adj, return_prob=True)
#     # old version runs without issue
#     # data_out = predict_crops(stata_in.reset_index().copy(), rev_prior_yr_df, logit_coefs, return_prob=True)
#     data_out_all = pd.concat((data_out_all, data_out.assign(year=year)))
#     # run prediction again to get just the selected crop
#     data_out = predict_crops_rand(stata_in.reset_index(drop=True).copy(), rev_prior_yr_df, logit_coefs_adj, return_prob=False)
#     data_out_sel = pd.concat((data_out_sel, data_out.assign(year=year)))


# %% [markdown]
# ## load crops predicted by modflow model

# %%
crop_all = pd.DataFrame()
for year in years:
    crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    crop_all = pd.concat((crop_all, crop_in.assign(year=year)))

# switch to underscores
crop_all = crop_all.replace(crop_dict)


# %% [markdown]
# ## compare predicted to observed

# %%
# columns relevant for comparison
stata_comp = stata_df[['parcel_id','year','crop_cat_last','crop_cat','crop','alternatives', 'chosen', 'dtwfa','wy_dc']]
# stata_comp = stata_df[['parcel_id','year','crop_cat','alternatives','pprob_stata']]
stata_comp.crop_cat = stata_comp.crop_cat.replace({'Mixed pasture and miscellaneous grasses':'Mixed_pasture'})
stata_comp.alternatives = stata_comp.alternatives.replace({'Mixed pasture and miscellaneous grasses':'Mixed_pasture'})


# %%
# join data on the selected crop from each for the parcel
# with the connected model we also want to check what crop choice last looked like
stata_sel = stata_comp[stata_comp.chosen==1][['parcel_id','year','crop_cat']].rename(columns={'crop_cat':'crop_choice_stata'})

comp_sel = stata_sel.merge(crop_all.rename(columns={'name':'crop_choice_MF'}))

# count how many where there is a difference
comp_sel['diff_crop'] = comp_sel.crop_choice_stata != comp_sel.crop_choice_MF
# comp_sel.to_csv(join(out_dir, 'selected_crop_stata_python.csv'))

# %%
# summarize crop count prediction by year and crop
comp_sel.groupby(['crop_choice_stata'])[['diff_crop']].sum().to_csv(join(fig_dir,'stata_modflow_crop_diff.csv'))

# %%
comp_sel.diff_crop.sum()/comp_sel.shape[0]


# %%
mf_count = crop_all.groupby(['year','name'])[['parcel_id']].count().reset_index()
stata_count = stata_sel.groupby(['year','crop_choice_stata'])[['parcel_id']].count().reset_index().rename(columns={'crop_choice_stata':'name'})
# mf_count = mf_count.rename(columns={'parcel_id':'mf_count'}).reset_index()
# stata_count = stata_count.rename(columns={'crop_choice_stata':'stata_count'}).reset_index()
# comp_count = mf_count.merge(stata_count)

# %%
count_comp = pd.concat((stata_count.assign(grp='stata'), mf_count.assign(grp='mf'))).reset_index()

# %%
# count_comp[count_comp.duplicated(['year','name','grp'])]
# # count_comp

# %%
sns.catplot(count_comp,x='year',y='parcel_id', col='name', col_wrap=3, hue='grp',
            kind='bar', #color='tab:blue',
            sharey=False
           # facet_kws={'sharey': False, 'sharex': True}
)

plt.savefig(join(fig_dir, 'parcel_count_by_crop_stata_modflow.png'), bbox_inches='tight')

