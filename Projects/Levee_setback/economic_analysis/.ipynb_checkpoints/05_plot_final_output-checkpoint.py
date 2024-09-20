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

# %%
# sns.relplot(df_econ_agg,x='year',y='total_value', col='crop', row='var', 
#            facet_kws={'sharey': False, 'sharex': True})

# plot the average profit and yield (not-weighted by acreage) 
g=sns.relplot(df_econ_agg,x='year',y='value', col='crop', row='var', 
           facet_kws={'sharey': False, 'sharex': True})


# %%

# %% [markdown]
# # Process water budget

# %%
df_all = pd.DataFrame()
# for year in run_years:
for year in [2015]:
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
