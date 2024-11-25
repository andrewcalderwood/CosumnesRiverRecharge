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
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# %%
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic'

# update to different modflow models here, next step is using the 20 year model



# %%
m_nam = 'input_write_2014_2020'
s_nam = 'input_write_2014_2020_R1'
# s_nam = 'input_write_2014_2020_R3'

model_ws = join(loadpth, m_nam)
scenario_ws = join(loadpth, s_nam)

model_out = join(model_ws, 'output_clean')
scenario_out = join(scenario_ws, 'output_clean')

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

# %%
# save data for Yusuke
df_econ  = pd.read_csv(join(model_out, 'annual_profit_yield_long.csv'),index_col=0)
df_econ_s  = pd.read_csv(join(scenario_out, 'annual_profit_yield_long.csv'),index_col=0)

# convert to wide format so Yusuke can plot easier
# df_econ_agg_wide = pd.read_csv(join(out_dir, 'annual_profit_yield_wide.csv'),index_col=0)

# %%
# merge data to co-plot or difference
diff_econ = df_econ.merge(df_econ_s, on=['crop','name','var','year'], suffixes=('_m','_s'))


# %%
var = 'profit'
# df_plt = df_econ_agg[df_econ_agg['var']==var].copy()
df_plt = diff_econ[diff_econ['var']==var].copy()
crops = df_plt.crop.unique()
fig,ax = plt.subplots(1, len(crops), sharey=False, figsize=(12,3), layout='constrained', dpi=300)

for n,crop in enumerate(crops):
    ax_n = ax[n]
    df_plt[df_plt.crop==crop].plot(x='year',y=['value_m','value_s'], ax=ax_n, kind='bar',legend=False)
    ax_n.set_title(crop)
    ax_n.set_xlabel('Year')
df_plt[df_plt.crop==crop].plot(x='year',y=['value_m','value_s'], ax=ax_n, kind='bar',legend=True)
plt.legend(['Baseline','Scenario'])
fig.supylabel(var.capitalize())

# %% [markdown]
# # Process water budget
# At minimum need to start showing irrigation on a sub-annual scale for an average field or average of fields
#
# Percolation is really only important to the modflow side of things so probably don't need to plot.

# %%
# save data for Yusuke
df_all = pd.read_csv(join(model_out, 'daily_WB_long.csv'), parse_dates=['date'])
df_all_s = pd.read_csv(join(scenario_out, 'daily_WB_long.csv'), parse_dates=['date'])

# convert to wide format so Yusuke can plot easier
# df_all_wide = pd.read_csv(join(out_dir, 'daily_WB_wide.csv'))

# %%
# merge data to co-plot or difference
diff_wb = df_all.merge(df_all_s, on=['name','date','var'], suffixes=('_m','_s'))
# diff_wb

# %%
# df_all.name.unique()

# %%
# 

# %%

var = 'GW_applied_water'
# for var in ['GW_applied_water']:
# for var in ['GW_applied_water', 'SW_applied_water','percolation']:

crops = df_all.name.unique()
# for crop in crops:
for crop in ['Vineyards']:
    fig,ax = plt.subplots( sharey=True, figsize=(12,3), layout='constrained', dpi=300)
    # plt_df = df_all[(df_all.name==crop)&(df_all['var']==var)]
    plt_df = diff_wb[(diff_wb.name==crop)&(diff_wb['var']==var)]

    # add in NA values to prevent line connection in dry season
    plt_df_na = plt_df.loc[plt_df.date.diff().dt.days>1].copy()
    plt_df_na.date -= pd.DateOffset(days=1)
    plt_df_na[['value','total_value']] = np.nan
    plt_df_na[['value_s','total_value_s','value_m','total_value_m']] = np.nan
    plt_df = pd.concat((plt_df, plt_df_na))

    # plt_df.plot(x='date',y='value', ax=ax, legend=False)
    # plt_df.plot(x='date',y=['value_m','value_s'], ax=ax, legend=False, )
    plt_df = plt_df.set_index('date').resample('MS').sum(numeric_only=True).reset_index()
    plt_df.plot(x='date',y=['value_m'], ax=ax, legend=False)
    plt_df.plot(x='date',y=['value_s'], ax=ax, legend=False, linestyle='-.' )
        
    fig.suptitle(crop)
    fig.supylabel(var.replace('_',' ')+'(m)')
    plt.xlabel(None)
    fig.supxlabel('Date')
# plt_df.plot(x='date',y=['value_m','value_s'], ax=ax, legend=True)
plt.legend(['Baseline','Scenario'])

    # plt.savefig(join(out_dir, var+'_'+crop+'.png'))
    # plt.close()
    

# %%
