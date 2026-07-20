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
# Script to evaluate pre-processed SWB to check relationship of AW, DTW and profit make sense.

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


import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from importlib import reload

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

add_path('../')


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
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
swb_version = '_no_p_o'
swb_ws = join(proj_dir, 'model_inputs', 'swb_rep', 'version'+swb_version)
print(swb_ws)


# %%
out_dir = join(swb_ws, 'output')
os.makedirs(out_dir, exist_ok=True)

# %%
all_run_dates = pd.read_csv(join(swb_ws, 'all_run_dates.csv'), parse_dates=['date'])
# years to sample output
run_years = all_run_dates[all_run_dates.use=='irrigation'].date.dt.year.values
# run_years = run_years[:-1]


# %% [markdown]
# ## spring dtw
# I forgot to save the intermediate file but right now, the swb representative output is based on depth to water profiles going from 0 to 300 ft with a decline of 5 ft each season.  
#
# There is also no crop or parcel info involved because we allow equal chance.
#

# %% [markdown]
# # Process water budget
# At minimum need to start showing irrigation on a sub-annual scale for an average field or average of fields
#
# Percolation is really only important to the modflow side of things so probably don't need to plot.

# %%
df_all = pd.DataFrame()
for year in run_years:
# for year in [2015]:
    # load SWB folder
    # crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    print('\n', year, end=' - ')
    for var in ['percolation','GW_applied_water', 'SW_applied_water']:
        print(var, end=',')
        # need to call final crop_soilbudget not rep value
        name = join(swb_ws, 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        with h5py.File(name) as dset:
            finished_crops = list(dset['array'].keys())
            print(finished_crops, end='.')
        for crop in finished_crops:
        # for crop in finished_crops[2]:
            # need dates for time series water budget output
            var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop)
            yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
            yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
            # get the total extent of the irrigation season (calculation period)
            strt_date = yield_start.min()
            end_date = yield_end.max()
            dates = pd.date_range(strt_date, end_date, freq='D')
            # extract output and convert to dataframe with ID columns
            arr = read_crop_arr_h5(crop, name)
            df = pd.DataFrame(arr, columns=dates)
            df['dtw_id'] = np.arange(0,len(df))
            # add parcel information back
            # df = pd.concat((df,crop_in[crop_in.name==pred_dict[crop]].reset_index(drop=True)),axis=1)
            # melt to long format for easier appending
            df = df.melt(var_name='date', id_vars='dtw_id')
            df['date'] = pd.to_datetime(df['date'])
            df = df.assign(crop=crop, year=year, var=var)
            # specify where there were irrigation events
            irr_dates = pd.date_range(yield_start.min(), yield_end.max(), freq=str(int(var_crops.gap_irr))+'D')
            irr_dates = pd.DataFrame(irr_dates, columns=['date'])
            irr_dates['irr'] = True
            df = df.merge(irr_dates, how='left')
            # concat to existing data
            df_all = pd.concat((df_all, df))


# %%
# distinguishing column for rep with and without pod
df_all['pod_bool'] = False
df_all.loc[df_all.dtw_id>29, 'pod_bool'] = True

# %%
# rename as econ for plotting reference
df_all = df_all.rename(columns={'parcel_id':'UniqueID'})#.merge(parcels[['UniqueID','acres']])

# scale value rates (1/acre) into total volumes (m3)
# df_all['total_value'] = df_all['value']*df_all.acres*4046.86


# %%
df_all_out = df_all.copy()
# for applied water we we want to drop the days that are not irrigated
# remove AW where there is no irrigation event planned (not if 0 zero irrigation)
df_all_out = df_all_out.loc[~((df_all_out['var'].str.contains('applied_water'))&(df_all_out.irr!=True))]
# based on a test review, there is no impact of this on the mean

# aggregate results by crop type
# add column for acres
df_all_out = df_all_out.groupby(['crop','date','var'])[['value']].agg({ 'value':'mean'})
# calculate rate, averaged by acres 
# this is the same as value but mutliplied by the 4046 scaling factor
# df_all_out['value_per_area'] = df_all_out.total_value / (df_all_out.acres*4046.86)
# save data for Yusuke
df_all_out.to_csv(join(out_dir, 'daily_WB_long.csv'))

# convert to wide format so Yusuke can plot easier
df_all_out_wide = df_all_out.pivot_table(index=['crop','date'], values=['value'], columns=['var'])
df_all_out_wide.to_csv(join(out_dir, 'daily_WB_wide.csv'))


# %% [markdown]
# the irrigation time series are not the most helpful if we are averaging they are more relevant if plot the distribution or CI to show how the range of DTW impacts this.

# %%

var = 'GW_applied_water'
# for var in ['GW_applied_water']:
for var in ['GW_applied_water', 'SW_applied_water','percolation']:

    crops = df_all.crop.unique()
    for crop in crops:
        # 1, len(run_years)-1,
        fig,ax = plt.subplots( sharey=True, figsize=(12,3), layout='constrained', dpi=300)
        plt_df = df_all[(df_all.crop==crop)&(df_all['var']==var)&(df_all['pod_bool']==False)]
    
        # add in NA values to prevent line connection in dry season
        plt_df_na = plt_df.loc[plt_df.date.diff().dt.days>1].copy()
        plt_df_na.date -= pd.DateOffset(days=1)
        plt_df_na[['value']] = np.nan
        plt_df = pd.concat((plt_df, plt_df_na))
    
        plt_df.plot(x='date',y='value', ax=ax, legend=False)
        # sns to include standard deviation lines
        # sns.lineplot(plt_df, x='date',y='value', errorbar = ("sd",1), ax=ax)
        # for n,year in enumerate(run_years[:-1]):
        #         ax_n = ax[n]
        #         plt_df[plt_df.year==year].plot(x='date',y='value',ax=ax_n,legend=False)
        #         ax_n.set_title(year)
            
        fig.suptitle(crop)
        fig.supylabel(var.replace('_',' ')+'(m)')
        plt.xlabel(None)
        fig.supxlabel('Date')
        plt.savefig(join(out_dir, var+'_'+crop+'.png'))
        plt.close()

# %%
# df_all

# %%
df_all['dtw_ft'] = df_all.dtw_id*10
df_all.loc[df_all.pod_bool, 'dtw_ft'] -= 300
# df_all

# %%
# check the average budget that it meets constraints
var = 'GW_applied_water'
# df_all_plt = df_all['var']==var
df_annual_sum = df_all.groupby(['crop','year','dtw_id','var','pod_bool']).sum(numeric_only=True).reset_index()
df_annual_sum = df_annual_sum.groupby(['crop','pod_bool','year','var']).mean(numeric_only=True).reset_index()


sns.catplot(df_annual_sum[df_annual_sum['var']==var],x='year',y='value', col='crop',  
            kind='bar', color='tab:blue',
            sharey=False
           # facet_kws={'sharey': False, 'sharex': True}
)

plt.savefig(join(out_dir, var+'_annual_total_m.png'))


# %%
# check to explore if there is a clear relationsihp between
# applied water and DTW
crop = 'Grape'
var = 'GW_applied_water'
plt_df = df_all[(df_all.crop==crop)&(df_all['var']==var)]
plt_df = plt_df[plt_df['pod_bool']==False]
plt_df = plt_df[plt_df['year']==2019]

df_annual_sum = plt_df.groupby(['dtw_id','crop','year','var']).agg({'value':'sum','dtw_ft':'mean'}).reset_index()

# %%
# check on outliers for irrigation
# df_annual_sum.value.describe()
# df_annual_sum.value.quantile([0,0.05,0.25,0.5,0.75,0.95,1])
# the majority of fields 90% + are in the normal range
# but we see a value at 0 and almost double the median (at max value)
# suggesting the SWB with optimization might fail occasionally

# %%
# # in addition to grouping by crop, need to group by field on some level to confirm
# # field properties aren't impacting
# plt_df.plot(x='dtw_ft',y='value', kind='scatter')
df_annual_sum.plot(x='dtw_ft',y='value', kind='scatter')
# chk_id = df_annual_sum.UniqueID.unique()
# n=600
# df_chk = df_annual_sum[df_annual_sum.UniqueID==chk_id[n]]

# # df_chk.plot(x='dtw_ft',y='value', kind='scatter')

# # it seems like most fields have a small range of AW values
# # but no big connection to DTW in any case
# # a few seem to show errors with too big of differences
# plt.ylim(0.4,.8)


# %% [markdown]
# Plotting the annual total for irrigation seems problematic when comparing with DTW in ft for a couple reasons.
# 1. pod_bool seems like it's broken because for True there is hardly any GW AW for all DTW then for False it is mostly GW AW up to 150 ft then switches which indicates an issue because POD bool = False means no SW should be utilized.
# 2. There is still no clear negative trend with increasing DTW and less irrigation unless it is obscured between years. Plotting an individua year shows problematic results as the overall trend is decreasing but there are anomalies where irrigation increases with DTW or where it spikes up to crazy values.
#
# Overal it seems that the SWB optimization is what is being problematic.
