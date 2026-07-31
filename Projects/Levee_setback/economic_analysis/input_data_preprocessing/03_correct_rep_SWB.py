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
# Script to load pre-processed SWB to correct relationship of AW, DTW and profit make sense.
#
# After discussing with Yusuke, we will hold on these corrections as we don't expect them to have a significant impact and we will wait to see if a review or co-author is concerned.
#
# TODO: move back code for reading in cost/profit into plt script and save intermediate csv if needed?

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
def add_pod_info(df_all, dtw_df):
    # distinguishing column for rep with and without pod
    df_all['pod_bool'] = False
    df_all['dtw_id'] = df_all['dtw_id_in']
    # it would nice to make the dtw_id part more automatic
    max_dtw_id = dtw_df.columns.astype(int).max()
    print('max dtw id', max_dtw_id)
    df_all.loc[df_all.dtw_id>max_dtw_id, 'pod_bool'] = True
    df_all.loc[df_all.dtw_id>max_dtw_id, 'dtw_id'] -= max_dtw_id
    
    # create dataframe to identify dtw_id to dtw_ft at starting point
    dtw_by_id = dtw_df.iloc[[0]].transpose().reset_index()
    dtw_by_id.columns=['dtw_id','dtw_ft']
    dtw_by_id.dtw_id = dtw_by_id.dtw_id.astype(int)
    df_all = df_all.merge(dtw_by_id)
    return df_all



# %%
df_all = add_pod_info(df_all, dtw_df)


# %%
# rename as econ for plotting reference
# df_all = df_all.rename(columns={'parcel_id':'UniqueID'})


# %% [markdown]
# the irrigation time series are not the most helpful if we are averaging they are more relevant if plot the distribution or CI to show how the range of DTW impacts this.

# %%
df_cost = pd.DataFrame()
for year in run_years:
# for year in [2015]:
    # load SWB folder
    # crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    print('\n', year, end=' - ')
    for var in ['cost', 'profit']:
        print(var, end=',')
        # need to call final crop_soilbudget not rep value
        name = join(swb_ws, 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        with h5py.File(name) as dset:
            finished_crops = list(dset['array'].keys())
            print(finished_crops, end='.')
        for crop in finished_crops:
            # extract output and convert to dataframe with ID columns
            arr = read_crop_arr_h5(crop, name)
            df = pd.DataFrame(arr, columns=['value'])
            df['dtw_id_in'] = np.arange(0,len(df))
            df = df.assign(crop=crop, year=year, var=var)
            df_cost = pd.concat((df_cost, df))

# %%
df_cost = add_pod_info(df_cost, dtw_df)
# dtw_df.iloc[0]
# df_cost.columns

# %%
var='profit'
pod_bool=False
sns.relplot(df_cost[(df_cost['var']==var)&(df_cost.pod_bool==pod_bool)],
            x='dtw_ft',y='value', hue='year', col='crop',kind='line', facet_kws={'sharey': False})

plt.savefig(join(out_dir, var+'_by_crop_year_dtw_pod'+str(pod_bool)+'.png'))
# out_dir
# reviewing profit shows that in most cases the profit decreases as expected with dtw, althought there are some exceptoins
# there is some levelign off or local decreases
# there is some odd behavior where misc. grain and hay irrig goes off (dry farming) and improves profit (unrealstic break in the line)

# for cost it makes sense for no pod, for pod true there is some odd behavior as if SW cost is not accounted

# %%
# check to explore if there is a clear relationsihp between
# applied water and DTW
crop = 'Corn'
var = 'GW_applied_water'
plt_df = df_all[(df_all.crop==crop)&(df_all['var']==var)]
plt_df = plt_df[plt_df['pod_bool']==False]



# %%
df_annual_sum_long = plt_df.groupby(['dtw_id','crop','year','var']).agg({'value':'sum','dtw_ft':'mean'}).reset_index()
df_annual_sum = df_annual_sum_long.pivot(columns='year',values='value', index=['dtw_ft'])
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

df_annual_sum.plot( ax=ax_n, )



# %%
var = 'SW_applied_water'
for var in ['SW_applied_water','GW_applied_water']:
    plt_df1 = df_all[(df_all['var']==var)]
    for pod_bool in [True, False]:
        plt_df = plt_df1[plt_df1['pod_bool']==pod_bool]
        
        crops = plt_df.crop.unique()
        # fig,ax = plt.subplots(1, len(crops), sharey=False, figsize=(12,3), layout='constrained', dpi=300)
        
        for n,crop in enumerate(crops):
            plt_df_c = plt_df[(plt_df.crop==crop)]
        
            df_annual_sum = plt_df_c.groupby(['dtw_id','year','var']).agg({'value':'sum','dtw_ft':'mean'}).reset_index()
            df_annual_sum = df_annual_sum.pivot(columns='year',values='value', index=['dtw_ft'])
        #     ax_n = ax[n]
        #     df_annual_sum.plot( ax=ax_n,legend=False )
        #     ax_n.set_title(crop)
        
        # ax[0].legend()
        # fig.supylabel(var.capitalize())
        
        # plt.savefig(join(out_dir, 
        #                  var+'_year_total_by_crop_pod'+str(plt_df.pod_bool.iloc[0])+'.png'), 
        #                  bbox_inches='tight')
        # plt.close()


# %% [markdown]
# Looking across the broader range of crops it's clearer that there are certain break points in each crop but generally they remain within a certain range that is about equal to WY impacts.
#
# These results to Yusuke might actually reverse the decision to keep in place the DTW variable adjustment since it appears that a code error was causing hte issue.
#
# Should update code to properly reference the pre-processed SWB profile in all cases to better complete this.

# %%
cost_var = pd.read_excel(join(data_dir, 'static_model_inputs_no_p_o.xlsx'), sheet_name = 'General')
# cost per acre-in
p_sw = cost_var[(cost_var.variable=='p_sw')&(cost_var.year==2019)].value.values[0]*12
# energy cost/kWh
p_e = cost_var[(cost_var.variable=='p_e')&(cost_var.year==2019)].value.values[0]
# phi is kWh per acre in per ft
p_gw = p_e*cost_var[cost_var.variable=='phi'].value.values[0]*12
# p_sw, p_gw*100

# %%
# # df_cost[(df_cost['var']=='profit')]
# # check whether cost captures SW AW
# cost_sw_chk = df_cost[(df_cost.pod_bool==True)&(df_cost['var']=='cost')]
# cost_sw_chk
# df_chk = df_all[df_all.pod_bool==True].groupby(['crop','year','var','dtw_ft'])[['value']].sum().reset_index()
# df_sw_chk = df_chk[(df_chk['var']=='SW_applied_water')].drop(columns=['var']).rename(columns={'value':'SW_AW'})
# df_gw_chk = df_chk[(df_chk['var']=='GW_applied_water')].drop(columns=['var']).rename(columns={'value':'GW_AW'})
# # the rep profiles have an assumed area of 1 acre so just need to convert to AF from m*acre
# df_sw_chk['SW_cost'] = (df_sw_chk.SW_AW/0.3048)*p_sw
# # # this doesn't account for the assumed 5 ft decline during irrigation season but should be minimal
# # # could subtract 2.5 ft to provide more accurate result
# df_gw_chk['GW_cost'] = (df_gw_chk.GW_AW/0.3048) *p_gw * df_gw_chk.dtw_ft
# cost_chk = cost_sw_chk.merge(df_sw_chk.merge(df_gw_chk))
# # cost_chk[cost_chk.value==0]
# # where cost is 0 the AW is 0 for both SW and GW
# # cost_chk['cost_chk'] = cost_chk.SW_cost+cost_chk.GW_cost
# # check that only one water source is used ine ach case comes clean
# cost_chk[cost_chk.SW_AW>0].GW_AW.sum(), cost_chk[cost_chk.GW_AW>0].SW_AW.sum()

# cost_chk[cost_chk.GW_AW>0]
# # something might be off in cost still, double check units

# %%
