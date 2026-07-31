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
swb_version = '_no_p_o_2019'
# swb_version = '_no_p_o'

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

# %%
# only load one year as it is static
dtw_df = pd.read_csv(join(swb_ws, 'field_SWB', 'dtw_ft_WY2015.csv'), index_col=0)

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
            df['dtw_id_in'] = np.arange(0,len(df))
            # add parcel information back
            # df = pd.concat((df,crop_in[crop_in.name==pred_dict[crop]].reset_index(drop=True)),axis=1)
            # melt to long format for easier appending
            df = df.melt(var_name='date', id_vars='dtw_id_in')
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
df_all['dtw_id'] = df_all['dtw_id_in']
# it would nice to make the dtw_id part more automatic
max_dtw_id = dtw_df.columns.astype(int).max()
print('max dtw id', max_dtw_id)
df_all.loc[df_all.dtw_id>max_dtw_id, 'pod_bool'] = True
df_all.loc[df_all.dtw_id>max_dtw_id, 'dtw_id'] -= max_dtw_id

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
# create dataframe to identify dtw_id to dtw_ft at starting point
dtw_by_id = dtw_df.iloc[[0]].transpose().reset_index()
dtw_by_id.columns=['dtw_id','dtw_ft']
dtw_by_id = dtw_by_id.dtw_id.astype(int)
df_all = df_all.merge(dtw_by_id)

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


# %% [markdown]
# We don't need to show each year so much as showing the dtw profile. The importance of by year is that hotter years will require greater irrigation. Maybe we should simply have DTW 50-150 by 20 ft then average.
#
# Another issue is that pod_bool False or True lets SW irrigation occurred. At least with pod_bool = True there is typically much less irrigation from GW but that shows there is still an issue in the code for optimization.
#
# The choose_water_source is agnostic of pod_bool, it is purely based on cost of DTW vs cost of surface water. It is after that where the sw_con is set to 0 if no POD is available. The issue is that after this gw_con and sw_con are set to 0 for whatever the cheaper cost is which would result in a zero irrigation total if pod_bool is False and surface water is cheaper.
#
# Fixed the water_source issue with no pod by adding fix to make only gw. The other issue with water source conflicting with con is that the wrong values would be written in.
#
# Alfalfa shows stronger declines but is more problematic as there are break points then it levels off/rises indicating it is sensitive on a larger scale but not a smaller scale.
# Corn is totally flat. Pasture has a sig decrease from 0-60 ft then totally flat which means that it would be flat as most DTW is >60 ft (same pattern for misc grain and hay).
# **Takeaway is that below the actual depths in the basin the irrigation totals are constant**

# %%
# check to explore if there is a clear relationsihp between
# applied water and DTW
crop = 'Corn'
var = 'GW_applied_water'
plt_df = df_all[(df_all.crop==crop)&(df_all['var']==var)]
plt_df = plt_df[plt_df['pod_bool']==False]



# %%
df_annual_sum = plt_df.groupby(['dtw_id','crop','year','var']).agg({'value':'sum','dtw_ft':'mean'}).reset_index()
df_annual_sum = df_annual_sum.pivot(columns='year',values='value', index=['dtw_ft'])
# df_annual_sum

# %%
# plt_df[(plt_df.dtw_id==5)&(plt_df.year==2018)].plot(x='date',y='value')
plt_chk = plt_df[(plt_df.year==2018)].copy()
plt_chk = plt_chk[plt_chk.irr==True]
sns.relplot(plt_chk, x='date',y='value', hue='dtw_id')

# %%
# # in addition to grouping by crop, need to group by field on some level to confirm
# # field properties aren't impacting
# plt_df.plot(x='dtw_ft',y='value', kind='scatter')
fig, ax_n = plt.subplots()
# for y in run_years:
#     plt_df_y = plt_df[plt_df['year']==y]
    
#     df_annual_sum = plt_df_y.groupby(['dtw_id','crop','year','var']).agg({'value':'sum','dtw_ft':'mean'}).reset_index()
# df_annual_sum.plot(x='dtw_ft',y='value', ax=ax_n, kind='scatter')
df_annual_sum.plot( ax=ax_n, )
# chk_id = df_annual_sum.UniqueID.unique()
# n=600
# df_chk = df_annual_sum[df_annual_sum.UniqueID==chk_id[n]]

# # df_chk.plot(x='dtw_ft',y='value', kind='scatter')

# # it seems like most fields have a small range of AW values
# # but no big connection to DTW in any case
# # a few seem to show errors with too big of differences
# plt.ylim(0.4,.8)


# %% [markdown]
# After correction for POD sorting in 2019 version it becomes clear that the water year has a bigger impact on irrigation optimization than the actual DTW. Essentially farmers are likely to irrigate to maximize their profit and bigger changes will occur due to crop choices between years.

# %% [markdown]
# Notes from the 2019 version.  
# Misc grain and hay had worse error as in 2019 there is a massive drop where the optimizer goes to zero for several DTW before recovering. Either we need to force a range of irrigation total for the season then the optimizer is primarily adjusting rates between events. One way might be to simply look at the expected irrigation deficit ignoring soil moisture to set irrigation but this would honestly be just as far off as soil water is significant.

# %% [markdown]
# Notes from year variable version.
# Plotting the annual total for irrigation seems problematic when comparing with DTW in ft for a couple reasons.
# 1. pod_bool seems like it's broken because for True there is hardly any GW AW for all DTW then for False it is mostly GW AW up to 150 ft then switches which indicates an issue because POD bool = False means no SW should be utilized.
# 2. There is still no clear negative trend with increasing DTW and less irrigation unless it is obscured between years. Plotting an individua year shows problematic results as the overall trend is decreasing but there are anomalies where irrigation increases with DTW or where it spikes up to crazy values.
#
# Overal it seems that the SWB optimization is what is being problematic.

# %%
var = 'SW_applied_water'
for var in ['SW_applied_water','GW_applied_water']:
    plt_df1 = df_all[(df_all['var']==var)]
    for pod_bool in [True, False]:
        plt_df = plt_df1[plt_df1['pod_bool']==pod_bool]
        
        crops = plt_df.crop.unique()
        fig,ax = plt.subplots(1, len(crops), sharey=False, figsize=(12,3), layout='constrained', dpi=300)
        
        for n,crop in enumerate(crops):
            plt_df_c = plt_df[(plt_df.crop==crop)]
        
            df_annual_sum = plt_df_c.groupby(['dtw_id','year','var']).agg({'value':'sum','dtw_ft':'mean'}).reset_index()
            df_annual_sum = df_annual_sum.pivot(columns='year',values='value', index=['dtw_ft'])
            ax_n = ax[n]
            df_annual_sum.plot( ax=ax_n,legend=False )
            ax_n.set_title(crop)
        
        ax[0].legend()
        fig.supylabel(var.capitalize())
        
        plt.savefig(join(out_dir, 
                         var+'_year_total_by_crop_pod'+str(plt_df.pod_bool.iloc[0])+'.png'), 
                         bbox_inches='tight')
        plt.close()


# %% [markdown]
# Looking across the broader range of crops it's clearer that there are certain break points in each crop but generally they remain within a certain range that is about equal to WY impacts.
#
# These results to Yusuke might actually reverse the decision to keep in place the DTW variable adjustment since it appears that a code error was causing hte issue.
#
# Should update code to properly reference the pre-processed SWB profile in all cases to better complete this.

# %%
