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
# Goal: Convert the final output into a more useable format for post-processing

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import glob
import sys
import time
from importlib import reload
import h5py

import pandas as pd
import numpy as np

# standard geospatial python utilities
import shapely
import geopandas as gpd

import matplotlib.pyplot as plt

# %%
import seaborn as sns

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
from functions.output_processing import get_wb_by_parcel
from functions.f_gw_dtw_extract import sample_dtw, avg_heads
import functions.Basic_soil_budget_monthly as swb

# initialize HDF5 files for the year
from functions.data_functions import init_h5


# %%
# from functions.swb_functions import run_swb
# import f_rep_swb_profit_opt
# reload(f_rep_swb_profit_opt)

from f_rep_swb_profit_opt import load_run_swb


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 


# %%
def read_crop_arr_h5(crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "r") as f:
        grp = f['array'] 
        # print(grp.keys())
        dset = grp[crop]
        arr = dset[:]
    return(arr)


# %%
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic'

# update to different modflow models here, next step is using the 20 year model
# base_model_ws = loadpth + 'crop_soilbudget'
m_nam = 'historical_simple_geology_reconnection'
m_nam = 'input_write_2014_2020'
model_ws = join(loadpth, m_nam)


# %%
# provide representative soil water budget folder
swb_ws = join(model_ws, 'rep_crop_soilbudget')
# os.makedirs(join(swb_ws, 'output'), exist_ok=True)

out_dir = join(model_ws, 'output_clean')
os.makedirs(out_dir, exist_ok=True)


# %%
# load parcel data for reference as needed
parcels = gpd.read_file(join(proj_dir,'Parcels shapefile/parcels.shp'))
parcels['area_m2'] = parcels.geometry.area
parcels['acres'] = parcels.area_m2/(43560*0.3048**2)
parcels.UniqueID = parcels.UniqueID.astype(int)

# %%
all_run_dates = pd.read_csv(join(model_ws, 'crop_modflow', 'all_run_dates.csv'), parse_dates=['date'])
# years to sample output
run_years = all_run_dates[all_run_dates.use=='irrigation'].date.dt.year

# %% [markdown]
# # Review data available

# %%
for year in [2015]:
    crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    # crop_in[crop_in.name=='Vineyards']
    # load SWB folder
    for var in ['profit', 'yield', 'percolation','GW_applied_water', 'SW_applied_water']:
        print(var)
        name = join(model_ws, 'crop_soilbudget', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        with h5py.File(name) as dset:
            finished_crops = list(dset['array'].keys())
            print(finished_crops)

# %% [markdown]
# # Process economic indicators

# %%
df_all = pd.DataFrame()
for year in run_years:
# for year in [2015]:
    # load SWB folder
    crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    print('\n', year, end=' - ')
    # for var in ['profit', 'yield', 'percolation','GW_applied_water', 'SW_applied_water']:
    for var in ['profit', 'yield']:
        print(var, end=',')
        name = join(model_ws, 'crop_soilbudget', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        with h5py.File(name) as dset:
            finished_crops = list(dset['array'].keys())
            print(finished_crops, end='.')
        for crop in finished_crops:
        # for crop in finished_crops[2]:
            # need dates for time series water budget output
            var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
            yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
            yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
            # get the total extent of the irrigation season (calculation period)
            strt_date = yield_start.min()
            end_date = yield_end.max()
            dates = pd.date_range(strt_date, end_date, freq='D')
            # extract output and convert to dataframe with ID columns
            arr = read_crop_arr_h5(crop, name)
            df = pd.DataFrame(arr, columns=['value']).assign(crop=crop, year=year, var=var)
            # add parcel information back
            df = pd.concat((df,crop_in[crop_in.name==pred_dict[crop]].reset_index()),axis=1)

            df_all = pd.concat((df_all, df))

# correct profit from negative to positive
df_all.loc[df_all['var']=='profit','value'] *= -1
# fix name before ID join
df_all = df_all.rename(columns={'parcel_id':'UniqueID'})
# rename as econ for plotting reference
df_econ = df_all.merge(parcels[['UniqueID','acres']])

# %%
# scale value rates (1/acre) into totals 
df_econ['total_value'] = df_econ['value']*df_econ.acres
# we want to aggregate yield and profit by the profit/acre and yield/acre to the total
# look at average rate, and summed total (scaled by acreage)
df_econ_agg = df_econ.groupby(['crop','name','var','year'])[['total_value','value']].agg({'total_value':'sum', 'value':'mean'})
# df_econ.groupby(['crop','name','var','year'])['total_value'].agg(['sum', 'mean'])

df_econ_agg = df_econ_agg.reset_index()
df_econ_agg.year = df_econ_agg.year.astype(str)
# df_econ_agg['end_date'] = pd.to_datetime(df_econ_agg.year.astype(str)+'-9-30')

# %%
# save data for Yusuke
df_econ_agg.to_csv(join(out_dir, 'annual_profit_yield_long.csv'))

# convert to wide format so Yusuke can plot easier
df_econ_agg_wide = df_econ_agg.pivot_table(index=['name','year'], values=['total_value','value'], columns=['var'])
df_econ_agg_wide.to_csv(join(out_dir, 'annual_profit_yield_wide.csv'))

# %%
# plot the total profit and yield after scaling by acreage
# sns.relplot(df_econ_agg,x='year',y='total_value', col='crop', row='var', 
#            facet_kws={'sharey': False, 'sharex': True})

# sns.catplot(df_econ_agg,x='year',y='total_value', col='crop', row='var', 
#             kind='bar', color='tab:blue',
#             sharey=False
# )

# %%
# plot the average profit and yield (not-weighted by acreage) 
# g=sns.relplot(df_econ_agg,x='year',y='value', col='crop', row='var', 
#            facet_kws={'sharey': False, 'sharex': True})

sns.catplot(df_econ_agg,x='year',y='value', col='crop', row='var', 
            kind='bar', color='tab:blue',
            sharey=False
           # facet_kws={'sharey': False, 'sharex': True}
)


# %%




# %%
var = 'profit'
df_plt = df_econ_agg[df_econ_agg['var']==var].copy()
crops = df_plt.crop.unique()
fig,ax = plt.subplots(1, len(crops), sharey=False, figsize=(12,3), layout='constrained', dpi=300)

for n,crop in enumerate(crops):
    ax_n = ax[n]
    df_plt[df_plt.crop==crop].plot(x='year',y='value', ax=ax_n, kind='bar',legend=False)
    ax_n.set_title(crop)
    ax_n.set_xlabel('Year')

fig.supylabel(var.capitalize())
plt.savefig(join(out_dir, var+'_field_avg.png'), bbox_inches='tight')

# %%
crop='Alfalfa'
crops = df_econ_agg.crop.unique()
dtw_mean_all = pd.DataFrame()
for crop in crops:
    fig,ax = plt.subplots(1, len(run_years)-1, sharey=True, figsize=(12,3), layout='constrained', dpi=300)
    
    for n,year in enumerate(run_years[:-1]):
        name = join(model_ws,'crop_soilbudget','field_dtw', 'dtw_ft_'+crop+'_'+str(year)+'.csv')
        if exists(name):
            dtw_arr = pd.read_csv(name,index_col=0,parse_dates=[0])
            ax_n = ax[n]
            dtw_arr_mean = dtw_arr.mean()
            dtw_mean_all = pd.concat((dtw_mean_all, pd.DataFrame(dtw_arr_mean).assign(year=year)))
            dtw_arr_mean.hist(ax=ax_n)
            ax_n.set_title(year)
        
    fig.suptitle(crop)
    ax[0].set_ylabel('Number of fields')
    fig.supxlabel('Mean depth to water (ft)')
    plt.savefig(join(out_dir, 'dtw_ft_histogram_'+crop+'.png'))
    plt.close()

# %%
# dtw_mean_all

# %%
fig,ax = plt.subplots(1, len(run_years)-1, sharey=True, figsize=(12,3), layout='constrained', dpi=300)

for n,year in enumerate(run_years[:-1]):
    dtw_mean_all.loc[dtw_mean_all.year==year,0].hist(ax=ax[n])
    ax[n].set_title(year)
 
ax[0].set_ylabel('Number of fields')
fig.supxlabel('Mean depth to water (ft)')
plt.savefig(join(out_dir, 'dtw_ft_histogram_all.png'))
plt.close()

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
    crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    print('\n', year, end=' - ')
    for var in ['percolation','GW_applied_water', 'SW_applied_water']:
        print(var, end=',')
        name = join(model_ws, 'crop_soilbudget', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        with h5py.File(name) as dset:
            finished_crops = list(dset['array'].keys())
            print(finished_crops, end='.')
        for crop in finished_crops:
        # for crop in finished_crops[2]:
            # need dates for time series water budget output
            var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
            yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
            yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
            # get the total extent of the irrigation season (calculation period)
            strt_date = yield_start.min()
            end_date = yield_end.max()
            dates = pd.date_range(strt_date, end_date, freq='D')
            # extract output and convert to dataframe with ID columns
            arr = read_crop_arr_h5(crop, name)
            df = pd.DataFrame(arr, columns=dates)
            # add parcel information back
            df = pd.concat((df,crop_in[crop_in.name==pred_dict[crop]].reset_index(drop=True)),axis=1)
            # melt to long format for easier appending
            df = df.melt(var_name='date', id_vars=crop_in.columns)
            df = df.assign(crop=crop, year=year, var=var)
            # concat to existing data
            df_all = pd.concat((df_all, df))

# # correct profit from negative to positive
# df_all.loc[df_all['var']=='profit','values'] *= -1
# # fix name before ID join
# df_all = df_all.rename(columns={'parcel_id':'UniqueID'})
# # rename as econ for plotting reference
# df_econ = df_all.merge(parcels[['UniqueID','area_m2']])

# %%
# rename as econ for plotting reference
df_all = df_all.rename(columns={'parcel_id':'UniqueID'}).merge(parcels[['UniqueID','acres']])
# scale value rates (1/acre) into totals 
df_all['total_value'] = df_all['value']*df_all.acres

# %%
df_all_out = df_all.copy().drop(columns=['pod_bool','pod'])

df_all_out = df_all_out.groupby(['name','date','var'])[['total_value','value']].agg({'total_value':'sum', 'value':'mean'})
# save data for Yusuke
df_all_out.to_csv(join(out_dir, 'daily_WB_long.csv'))

# convert to wide format so Yusuke can plot easier
df_all_out_wide = df_all_out.pivot_table(index=['name','date'], values=['total_value','value'], columns=['var'])
df_all_out_wide.to_csv(join(out_dir, 'daily_WB_wide.csv'))

# %%
# this runs pretty slowly because there is daily data which in theory can be dropped since 
# we really only care about when there is action (irrigation events)
crop = 'Alfalfa'
# subset to specific crop and variable of interest for consistent axes/numbers
plt_df = df_all[(df_all.crop==crop)&(df_all['var']=='percolation')].copy()
plt_df = df_all[(df_all['var']=='GW_applied_water')].copy()

# sns.relplot(plt_df, x='date',y='value', col='year', row='crop',
#            facet_kws={'sharey': True, 'sharex': 'col'}, 
#             kind='line', err_style="bars"
#            )
# plt_df



# %%

var = 'GW_applied_water'
crops = df_econ_agg.crop.unique()
for crop in crops:
    # 1, len(run_years)-1,
    fig,ax = plt.subplots( sharey=True, figsize=(12,3), layout='constrained', dpi=300)
    plt_df = df_all[(df_all.crop==crop)&(df_all['var']==var)]

    # add in NA values to prevent line connection in dry season
    plt_df_na = plt_df.loc[plt_df.date.diff().dt.days>1].copy()
    plt_df_na.date -= pd.DateOffset(days=1)
    plt_df_na[['value','total_value']] = np.nan
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

# %% [markdown]
#
# ### applied water
# - for grape we have only SW applied which shows that something is wrong because realistically most vineyards use groundwater
# - alfalfa shows no groundwater used either
# - corn is also only surface water
# - misc grain and hay is only surface water
#
# **When going back to 03b_summarize_output.py there should only be gw applied water and no sw so they got mixed up**
#
# ### percolation
# - misc grain and hay has very little percolation except 2017
# - grape shows a bunch in 2015 then little in the rest except 2017
# - corn looks good with just a few late season irrigation event causing recharge
# - alfalfa has recharge in late irrigation events
