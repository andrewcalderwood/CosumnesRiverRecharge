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
# for batch file case we want to ignore the consant h5py deprecation warning
import warnings

warnings.filterwarnings("ignore")
# warnings.filterwarnings("ignore", category=DeprecationWarning) 

# %%
# updated version specifies concept_name and copy_files here so it can be easily
# seen as these are the main update to make in a script
in_data = sys.argv

if 'ipykernel' in in_data[0]:
    # m_nam = 'input_write_2014_2020_R3'
    # m_nam = 'input_write_2014_2020'
    m_nam = 'input_write_2014_2025_R204'
    # m_nam = 'input_write_2016_2022_R200'

    # input_name = 'static_model_inputs.xlsx'
    input_name = 'static_model_inputs_no_p_o.xlsx'

else:
    m_nam = in_data[1]
    # added option to specify different static_model_inputs to allow easier adjustment of operating costs, revenue
    input_name = in_data[2]


print('sys.argv[1] (m_nam) is...')
print(m_nam)

print('sys.argv[2] (input_name) is...')
print(input_name)

t_start = time.time()

# %%
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'D://WRDAPP/GWFlowModel/Cosumnes/Economic'

# update to different modflow models here
# m_nam = 'input_write_2014_2020'
# m_nam = 'input_write_2014_2022'

# scenario = '_R20' # pumping constraint
# scenario = '_R4' # 90/20 floodplain
# # scenario = '_R3' # 6x existing diversion for MAR vineyard

# scenario = '_R200' # 200 represents no p_o
# scenario = '_R203' # no p_o and 6x MAR
# scenario = '_R204' # no p_o and 90/20 flloodplain

# scenario=''
# m_nam = m_nam+scenario


model_ws = join(loadpth, m_nam)


# %%
# provide representative soil water budget folder
# swb_ws = join(model_ws, 'rep_crop_soilbudget')
# we should be assessing the final clean output
swb_ws = join(model_ws, 'crop_soilbudget')

out_dir = join(model_ws, 'output_clean')
os.makedirs(out_dir, exist_ok=True)


# %%
# save log by date so we can see old versions
os.makedirs(join(out_dir, 'log'), exist_ok=True)
if 'ipykernel' in in_data[0]:
    print('print to jupyter')
else:
    print('printing remaining output to log file')
    sys.stdout = open(join(out_dir, 'log', 'plot_output_only_log_'+str(pd.to_datetime('today').date())+'.txt'), 'w')


# %%
# load parcel data for reference as needed
parcels = gpd.read_file(join(proj_dir,'Parcels shapefile/parcels.shp'))
parcels['area_m2'] = parcels.geometry.area
parcels['acres'] = parcels.area_m2/(43560*0.3048**2)
parcels.UniqueID = parcels.UniqueID.astype(int)

# %%
all_run_dates = pd.read_csv(join(model_ws, 'crop_modflow', 'all_run_dates.csv'), parse_dates=['date'])
# years to sample output
run_years = all_run_dates[all_run_dates.use=='irrigation'].date.dt.year.values
# run_years = run_years[:-1]

# %%
# check which files are available and shorten as needed
files = pd.DataFrame(glob.glob(join(out_dir, 'profit_yield_long_*.csv')), columns=['fn'])
files['name'] = files.fn.apply(basename)
files['year'] = files.name.str.extract(r'(\d{4})').astype(int)

# subset run years to those with output
run_years = run_years[pd.Series(run_years).isin(files.year)]
# should set up a way to delete old output to avoid mixing old/new in plots

# %%
# temp for while model is still running
# run_years = np.arange(2001, 2011)
# run_years = np.arange(2001, 2020)

# %% [markdown]
# # Plot economic indicators

# %%
# general limits for plotting
ncol_max = 4
# subtract 1 from run_years since last doesn't cover full period
nyears = len(run_years)-1
ncol = np.min((ncol_max, nyears))
nrow = int(np.ceil(nyears/ncol_max))

figsize= (4*ncol, 3*nrow)

# %%
df_econ_agg = pd.DataFrame()

for year in run_years:
    df_econ_in = pd.read_csv(join(out_dir, 'profit_yield_long_'+str(year)+'.csv'),index_col=0)
    df_econ_agg = pd.concat((df_econ_agg, df_econ_in))

df_econ_agg = df_econ_agg.reset_index()
# df_all

# %%
df_econ_agg.to_csv(join(out_dir, 'annual_profit_yield_long.csv'))


# %%
# plot the total profit and yield after scaling by acreage
# sns.relplot(df_econ_agg,x='year',y='total_value', col='crop', row='var', 
#            facet_kws={'sharey': False, 'sharex': True})

# sns.catplot(df_econ_agg,x='year',y='total_value', col='crop', row='var', 
#             kind='bar', color='tab:blue',
#             sharey=False
# )

# %%
# df_econ_agg

# %%
# plot the average profit and yield (not-weighted by acreage) 
# g=sns.relplot(df_econ_agg,x='year',y='value', col='crop', row='var', 
#            facet_kws={'sharey': False, 'sharex': True})

sns.catplot(df_econ_agg,x='year',y='total_value', col='crop', row='var', 
            kind='bar', color='tab:blue',
            sharey=False
           # facet_kws={'sharey': False, 'sharex': True}
)

plt.savefig(join(out_dir,'profit_yield_annual_total_m.png'))


# %%
# var = 'profit'
for var in ['profit', 'yield']:
# for var in ['profit']:
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
    plt.close()

# %% [markdown]
# # Depth to water
# There are two options for plotting:
# 1. The Spring DTW (march avg) which informs the crop choice
# 2. the seasonal average value for crops that are irrigated
#
# To include fallow it makes the most sense to use the spring value because that is what is driving the decision rather than the final seasonal values

# %%
os.makedirs(join(out_dir,'DTW'), exist_ok=True)

# %%
crop_name_dict = pd.read_excel(join(proj_dir,'model_inputs',input_name),sheet_name='Name_dict')

# %% [markdown]
# ## spring dtw

# %%
dtw_spring = pd.DataFrame()
crop_all = pd.DataFrame()
for n,year in enumerate(run_years):
    name = join(model_ws, 'rep_crop_soilbudget','field_SWB','modflow_spring_dtw_ft_WY'+str(year)+'.csv')
    dtw_spring_year = pd.read_csv(name,index_col=0)
    dtw_spring = pd.concat((dtw_spring, dtw_spring_year.assign(year=year)))
    
    name = join(model_ws, 'rep_crop_soilbudget','field_SWB','crop_parcels_'+str(year)+'.csv')
    crop_year = pd.read_csv(name,index_col=0)
    crop_all = pd.concat((crop_all, crop_year.assign(year=year)))


# %%
# add info on crops and acreage to spring dtw
# there is an issue with crop_all dropping all rows because I added a column to save the dtw_ft used in the crop choice
dtw_spring_ref = dtw_spring.merge(parcels[['UniqueID','acres']])
dtw_spring_ref = dtw_spring_ref.merge(crop_all.rename(columns={'parcel_id':'UniqueID','dtw_ft':'dtw_ft_fall'}))
# # # add simpler crop name
dtw_spring_ref = dtw_spring_ref.merge(crop_name_dict[['Crop','pred_name']].rename(columns={'pred_name':'name','Crop':'crop'}), how='left')
# where there is not a short name use the existing
dtw_spring_ref.loc[dtw_spring_ref.crop.isna(),'crop'] = dtw_spring_ref.loc[dtw_spring_ref.crop.isna(),'name'] 

dtw_spring_ref.to_csv(join(out_dir, 'dtw_ft_spring_all.csv')) 


# %%
# df_plt.hist('dtw_ft', weights = df_plt.acres, histtype=u'step',)


# %%
crops = dtw_spring_ref.crop.unique()

for crop in crops:
# for crop in crops[[1]]:
    fig,ax = plt.subplots(nrow, ncol, sharey=True, figsize=figsize, layout='constrained', dpi=300)
    
    for n,year in enumerate(run_years):
        if nrow>1:
            ax_n = ax[int(n/ncol), int(n%ncol)]
        elif nrow==1:
            ax_n = ax[int(n%ncol)]
        # keeps column of unique ID as index
        df_plt = dtw_spring_ref[(dtw_spring_ref.crop==crop)&(dtw_spring_ref.year==year)]
        if len(df_plt)>0:
            df_plt.hist('dtw_ft', ax=ax_n, weights = df_plt.acres)
            # dtw_arr_mean.plot.kde(ax=ax_n) # version with kernel density plot
            ax_n.set_title(year)
        
    fig.suptitle(crop)
    fig.supylabel('Area (acres)')
    # fig.supylabel('Density')
    fig.supxlabel('Mean depth to water (ft)')
    plt.savefig(join(out_dir, 'DTW', 'spring_dtw_ft_histogram_acres_'+crop+'.png'))
    plt.close()
# dtw_mean_all.index.name = 'UniqueID'

# %% [markdown]
# ## seasonal avg dtw

# %%
crop='Alfalfa'
crops = df_econ_agg.crop.unique()
dtw_mean_all = pd.DataFrame()

for crop in crops:
    fig,ax = plt.subplots(nrow, ncol, sharey=True, figsize=figsize, layout='constrained', dpi=300)
    
    for n,year in enumerate(run_years):
        name = join(model_ws,'crop_soilbudget','field_dtw', 'dtw_ft_'+crop+'_'+str(year)+'.csv')
        if exists(name):
            dtw_arr = pd.read_csv(name,index_col=0,parse_dates=[0])
            if nrow>1:
                ax_n = ax[int(n/ncol), int(n%ncol)]
            elif nrow==1:
                ax_n = ax[int(n%ncol)]
            # keeps column of unique ID as index
            dtw_arr_mean = dtw_arr.mean()
            dtw_mean_all = pd.concat((dtw_mean_all, pd.DataFrame(dtw_arr_mean, columns=['dtw_ft']).assign(year=year, crop=crop)))
            dtw_arr_mean.hist(ax=ax_n)
            # dtw_arr_mean.plot.kde(ax=ax_n) # version with kernel density plot
            ax_n.set_title(year)
        
    fig.suptitle(crop)
    fig.supylabel('Number of fields')
    # fig.supylabel('Density')
    fig.supxlabel('Mean depth to water (ft)')
    plt.savefig(join(out_dir, 'DTW','season_avg_dtw_ft_histogram_'+crop+'.png'))
    plt.close()
dtw_mean_all.index.name = 'UniqueID'

# %%
dtw_mean_all_ref = dtw_mean_all.copy()
dtw_mean_all_ref.index = dtw_mean_all_ref.index.astype(int)

# add info on acreage
dtw_mean_all_ref = dtw_mean_all_ref.join(parcels.set_index('UniqueID')[['acres']])


# %%
dtw_mean_all_ref.to_csv(join(out_dir, 'dtw_ft_mean_all.csv')) # new version with acreage
# dtw_mean_all.to_csv(join(out_dir, 'dtw_ft_mean_all.csv'))

# %%
crop='Alfalfa'
crops = df_econ_agg.crop.unique()

for crop in crops:
    fig,ax = plt.subplots(nrow, ncol, sharey=True, figsize=figsize, layout='constrained', dpi=300)
    
    for n,year in enumerate(run_years):
        if nrow>1:
            ax_n = ax[int(n/ncol), int(n%ncol)]
        elif nrow==1:
            ax_n = ax[int(n%ncol)]
        # keeps column of unique ID as index
        df_plt = dtw_mean_all_ref[(dtw_mean_all_ref.crop==crop)&(dtw_mean_all_ref.year==year)]
        if len(df_plt)>0:
            df_plt.hist('dtw_ft', ax=ax_n, weights = df_plt.acres)
            # dtw_arr_mean.plot.kde(ax=ax_n) # version with kernel density plot
            ax_n.set_title(year)
        
    fig.suptitle(crop)
    fig.supylabel('Area (acres)')
    # fig.supylabel('Density')
    fig.supxlabel('Mean depth to water (ft)')
    plt.savefig(join(out_dir, 'DTW', 'season_avg_dtw_ft_histogram_acres_'+crop+'.png'))
    plt.close()
dtw_mean_all.index.name = 'UniqueID'

# %%

fig,ax = plt.subplots(nrow,ncol, sharey=True, figsize=figsize, layout='constrained', dpi=300)

for n,year in enumerate(run_years[:-1]):
    if nrow>1:
        ax_n = ax[int(n/ncol), int(n%ncol)]
    elif nrow==1:
        ax_n = ax[int(n%ncol)]    
    dtw_mean_all.loc[dtw_mean_all.year==year,'dtw_ft'].hist(ax=ax_n)
    ax_n.set_title(year)
 
fig.supylabel('Number of fields')
fig.supxlabel('Mean depth to water (ft)')
plt.savefig(join(out_dir, 'dtw_ft_histogram_all.png'))
plt.close()

# %% [markdown]
# # Process water budget
# At minimum need to start showing irrigation on a sub-annual scale for an average field or average of fields
#
# Percolation is really only important to the modflow side of things so probably don't need to plot.

# %%

# up until june 2026 this code appears to have been referencing the rep_crop_soilbudget which
# is not correct, it should reference the clean final soil budget

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
            # add parcel information back
            df = pd.concat((df,crop_in[crop_in.name==pred_dict[crop]].reset_index(drop=True)),axis=1)
            # melt to long format for easier appending
            df = df.melt(var_name='date', id_vars=crop_in.columns)
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
# rename as econ for plotting reference
df_all = df_all.rename(columns={'parcel_id':'UniqueID'}).merge(parcels[['UniqueID','acres']])
# for applied water we we want to drop the days that are not irrigated

# scale value rates (1/acre) into total volumes (m3)
df_all['total_value'] = df_all['value']*df_all.acres*4046.86


# %%
df_all_out = df_all.copy().drop(columns=['pod_bool','pod'])
# aggregate results by crop type
# add column for acres
df_all_out = df_all_out.groupby(['name','date','var'])[['total_value','value','acres']].agg({'total_value':'sum', 'value':'mean','acres':'sum'})
# calculate rate, averaged by acres 
# this is the same as value but mutliplied by the 4046 scaling factor
# df_all_out['value_per_area'] = df_all_out.total_value / (df_all_out.acres*4046.86)
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

# %%

var = 'GW_applied_water'
# for var in ['GW_applied_water']:
for var in ['GW_applied_water', 'SW_applied_water','percolation']:

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

# %%
# check the average budget that it meets constraints
var = 'GW_applied_water'
# df_all_plt = df_all['var']==var
df_annual_sum = df_all.groupby(['crop','year','name','UniqueID','var','pod']).sum(numeric_only=True).reset_index()
df_annual_sum = df_annual_sum.groupby(['crop','name','year','var']).mean(numeric_only=True).reset_index()


sns.catplot(df_annual_sum[df_annual_sum['var']==var],x='year',y='value', col='crop',  
            kind='bar', color='tab:blue',
            sharey=False
           # facet_kws={'sharey': False, 'sharex': True}
)

plt.savefig(join(out_dir, var+'_annual_total_m.png'))


# %%
t_final = time.time()
print('Total time was %.2f hours' %((t_final-t_start)/3600))

sys.stdout.flush()

# %%
# check to explore if there is a clear relationsihp between
# applied water and DTW
# crop = 'Grape'
# var = 'GW_applied_water'
# # plt_df = df_all[(df_all.crop==crop)&(df_all['var']==var)]
# # df_all

# %%
# df_annual_sum = plt_df.groupby(['UniqueID','crop','year','var']).agg({'value':'sum','dtw_ft':'mean'}).reset_index()

# %%
# check on outliers for irrigation
df_annual_sum.value.describe()
df_annual_sum.value.quantile([0,0.05,0.25,0.5,0.75,0.95,1])
# the majority of fields 90% + are in the normal range
# but we see a value at 0 and almost double the median (at max value)
# suggesting the SWB with optimization might fail occasionally

# %%
# # in addition to grouping by crop, need to group by field on some level to confirm
# # field properties aren't impacting
# # plt_df.plot(x='dtw_ft',y='value', kind='scatter')
# # df_annual_sum.plot(x='dtw_ft',y='value', kind='scatter')
# chk_id = df_annual_sum.UniqueID.unique()
# n=600
# df_chk = df_annual_sum[df_annual_sum.UniqueID==chk_id[n]]

# # df_chk.plot(x='dtw_ft',y='value', kind='scatter')

# # it seems like most fields have a small range of AW values
# # but no big connection to DTW in any case
# # a few seem to show errors with too big of differences

