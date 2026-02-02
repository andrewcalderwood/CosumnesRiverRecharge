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
# loadpth = 'D://WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Economic'


# %% [markdown]
# The first few plots are scenario specific, the last plots are across all listed scenarios

# %%
# update to different modflow models here
m_nam = 'input_write_2014_2020'
m_nam = 'input_write_2014_2022'

scenario = '_R20' # pumping constraint
scenario = '_R4' # 90/20 floodplain
scenario = '_R3' # 6x existing diversion for MAR vineyard

scenario = '_R203' # no p_o, 6x MAR
scenario = '_R204' # no p_o, 90/20 floodplain

s_nam = m_nam+scenario

m_version = '' # original baseline with p_o
m_version = '_R200' # no p_o # specified later?
# need to specify base version for easier scenario switching
base_model_ws = join(loadpth, m_nam)

model_ws = join(loadpth, m_nam+m_version)
scenario_ws = join(loadpth, s_nam)

model_out = join(model_ws, 'output_clean')
scenario_out = join(scenario_ws, 'output_clean')

# %%
scenario_summary = pd.read_excel(join(data_dir, 'scenario_summary.xlsx'))
# # subset to scenarios of interest
scenario_summary = scenario_summary[scenario_summary.base_m_nam==m_nam]
# # to identify file paths
scenario_summary['extension'] = scenario_summary.m_nam.str.extract(r'(_R\d{1,2})')
scenario_summary = scenario_summary.set_index('m_nam')


# %%
scenario_dir = join(scenario_out, 'baseline_comparison')
os.makedirs(scenario_dir, exist_ok=True)

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
# general limits for plotting
ncol_max = 4
# subtract 1 from run_years since last doesn't cover full period
nyears = len(run_years)-1
ncol = np.min((ncol_max, nyears))
nrow = int(np.ceil(nyears/ncol_max))

figsize= (4*ncol, 3*nrow)

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
# df_crop = df_plt[df_plt.crop==crop]
# df_crop.set_index('year').reindex(df_plt.year.unique())

# %%
var = 'yield'
var='profit'
for val in ['total_value', 'value']:
    for var in ['yield','profit']:
        # df_plt = df_econ_agg[df_econ_agg['var']==var].copy()
        df_plt = diff_econ[diff_econ['var']==var].copy()
        crops = df_plt.crop.unique()
        fig,ax = plt.subplots(1, len(crops), sharey=False, figsize=(12,3), layout='constrained', dpi=300)
        
        for n,crop in enumerate(crops):
            ax_n = ax[n]
            df_plt[df_plt.crop==crop].plot(x='year',y=[val+'_m', val+'_s'], ax=ax_n, kind='bar',legend=False)
            ax_n.set_title(crop)
            ax_n.set_xlabel('Year')
            # ax_n.set_xticks(df_plt.year.unique().astype(str))
        
        
        # df_crop = df_plt[df_plt.crop==crop]
        # add legend to the last plot
        df_plt[df_plt.crop==crop].plot(x='year',y=[val+'_m',val+'_s'], ax=ax_n, kind='bar',legend=True)
        plt.legend(['Baseline','Scenario'])
        fig.supylabel(var.capitalize())
        plt.savefig(join(scenario_dir, var+'_field_avg_'+val+'.png'), bbox_inches='tight')
        plt.close()

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
df_all.name.unique()

# %%
crops = df_all.name.unique()

var = 'GW_applied_water'
# for var in ['GW_applied_water']:
for var in ['GW_applied_water', 'SW_applied_water','percolation']:

    for crop in crops:
    # for crop in ['Corn, sorghum or Sudan']:
    # for crop in ['Alfalfa and alfalfa mixtures']:
    # for crop in ['Vineyards']:
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
    
        plt.savefig(join(scenario_dir, var+'_'+crop+'.png'))
        plt.close()


# %% [markdown]
# ## plot DTW
# Switched to plotting difference in spring levels as this is available for fallow and irrigated crops  
# The part below here plots for all three scenarios

# %%
# load scenario data at once to simplify plotting 
scenarios = ['','_R3','_R4']
scenarios = ['_R200','_R203','_R204']

dtw_all_scenario = pd.DataFrame()
for scenario in scenarios:
    dtw_mean_all = pd.read_csv(join(base_model_ws+scenario,'output_clean', 'dtw_ft_spring_all.csv'),index_col=0)
    dtw_all_scenario = pd.concat((dtw_all_scenario, dtw_mean_all.assign(scenario=scenario, scen_name=scenario_summary.loc[basename(base_model_ws)+scenario, 'plot_name'])))

    # could improve with scenario summary sheet
scenario_names = dtw_all_scenario.scen_name.unique()


# %%
out_dir = join(model_ws, 'output_clean')
# index is UniqueID
# average dtw across time for each parcel by year 
# seasonal average by parcel
# dtw_mean_all = pd.read_csv(join(model_out, 'dtw_ft_mean_all.csv'),index_col=0)
# dtw_mean_all_s = pd.read_csv(join(scenario_out, 'dtw_ft_mean_all.csv'),index_col=0)
# spring average, includes fallow
dtw_mean_all = pd.read_csv(join(model_out, 'dtw_ft_spring_all.csv'),index_col=0)
dtw_mean_all_s = pd.read_csv(join(scenario_out, 'dtw_ft_spring_all.csv'),index_col=0)

dtw_mean_all_mean = dtw_mean_all.groupby(['year','crop']).mean(numeric_only=True).reset_index()
dtw_mean_all_s_mean = dtw_mean_all_s.groupby(['year','crop']).mean(numeric_only=True).reset_index()



# %%
for crop in dtw_mean_all.crop.unique():
    fig,ax = plt.subplots()

    dtw_mean_all_mean[dtw_mean_all_mean.crop==crop].plot(x='year',y='dtw_ft', legend=False, ax=ax)
    dtw_mean_all_s_mean[dtw_mean_all_s_mean.crop==crop].plot(x='year',y='dtw_ft', legend=False,ax=ax)

    plt.legend(['Baseline','Scenario'])
    plt.ylabel('Depth to Water (ft)')
    plt.title(crop)
    plt.savefig(join(scenario_dir, 'dtw_ft_yearly_'+crop+'.png'))
    plt.close()

# %%
# dtw_mean_all[(dtw_mean_all.crop==crop)&(dtw_mean_all.year==year)]
# dtw_mean_all.year.max()


# %%
from scipy.stats import gaussian_kde

# crop='Grape'
# df_plt = dtw_mean_all[(dtw_mean_all.crop==crop)&(dtw_mean_all.year==year)]
# # df_plt.dtw_ft.values
# kernel = gaussian_kde(df_plt.dtw_ft)
# # kernel = gaussian_kde(df_plt.dtw_ft, weights = df_plt.acres)
# eval_points = np.linspace(np.min(df_plt.dtw_ft), np.max(df_plt.dtw_ft))

# plt.plot(kernel.pdf(eval_points))

# %%
# the current density plot is basing density on the number of parcels, the actual density should use acreage

crops = dtw_mean_all.crop.unique()

for crop in crops:
# for crop in crops[[1]]:

    fig,ax = plt.subplots(nrow, ncol, sharey=True, figsize=figsize, layout='constrained', dpi=300)
    
    for n,year in enumerate(run_years):
            if nrow>1:
                ax_n = ax[int(n/ncol), int(n%ncol)]
            elif nrow==1:
                ax_n = ax[int(n%ncol)]
            df_plt = dtw_mean_all[(dtw_mean_all.crop==crop)&(dtw_mean_all.year==year)]
            ax_n.set_ylabel('')
            if len(df_plt)>0:
                # df_plt.dtw_ft.plot.kde(ax=ax_n, label='Baseline Unweighted') # version with kernel density plot
                df_plt_s = dtw_mean_all_s[(dtw_mean_all_s.crop==crop)&(dtw_mean_all_s.year==year)]
                # found a way to plot histogram with just line that helps show results since KDE doesn't show acreage well
                df_plt.dtw_ft.hist(ax=ax_n, label='Baseline', weights = df_plt.acres, histtype=u'step',) 
                df_plt_s.dtw_ft.hist(ax=ax_n, label='Scenario', weights = df_plt_s.acres, histtype=u'step',) 

                # df_plt_s.dtw_ft.plot.kde(ax=ax_n, label='Scenario') # version with kernel density plot
                # weighted kde for baseline
                eval_points = np.linspace(np.min(df_plt.dtw_ft), np.max(df_plt.dtw_ft))
                kernel = gaussian_kde(df_plt.dtw_ft, weights = df_plt.acres)
                # ax_n.plot(kernel.pdf(eval_points), label='Baseline')
                # weighted kde for scenario
                eval_points = np.linspace(np.min(df_plt_s.dtw_ft), np.max(df_plt_s.dtw_ft))
                kernel = gaussian_kde(df_plt_s.dtw_ft, weights = df_plt_s.acres)
                # ax_n.plot(kernel.pdf(eval_points), label='Scenario')
            ax_n.set_title(year)
    ax[0,0].legend()
        
    fig.suptitle(crop)
    fig.supylabel('Crop Area (acres)')
    # fig.supylabel('Density')
    fig.supxlabel('Mean depth to water (ft)')
    plt.savefig(join(scenario_dir, 'dtw_ft_kde_'+crop+'.png'))
    plt.savefig(join(scenario_dir, 'dtw_ft_hist_line_'+crop+'.png'))
    plt.close()


# %%
join(scenario_dir, 'dtw_ft_kde_'+crop+'.png')

# %% [markdown]
# Adding weighting definitely shifts the results, especially since they are for the actual DTW extent. What doesn't make sense is the PDf still has values with negative they just aren't evaluated.
# - need to better understand KDE so ask Yusuke -> Yusuke was not very familiar with this so suggested whatever I found made sense.

# %%
os.makedirs(join(model_out, 'scenario_comparison'), exist_ok=True)

# %%
fig,ax = plt.subplots(nrow, ncol, sharey=True, figsize=figsize, layout='constrained', dpi=300)
# fig,ax = plt.subplots()
for sn in scenario_names:
    df_dtw_plt = dtw_all_scenario[dtw_all_scenario.scen_name==sn].copy()
    for n,year in enumerate(run_years):
            if nrow>1:
                ax_n = ax[int(n/ncol), int(n%ncol)]
            elif nrow==1:
                ax_n = ax[int(n%ncol)]
            df_plt = df_dtw_plt[(df_dtw_plt.year==year)]
            df_plt.dtw_ft.hist(label=sn, ax=ax_n, weights = df_plt.acres, histtype=u'step',)

ax[0,0].legend()
fig.supylabel('Field area (acres)')

fig.supxlabel('DTW (ft)')
plt.savefig(join(model_out, 'scenario_comparison', 'dtw_ft_spring_line_hist_acres_all.png'))
plt.close()

# %%

# spring average, includes fallow
dtw_mean_all = pd.read_csv(join(model_out, 'dtw_ft_spring_all.csv'),index_col=0)
dtw_mean_all_s = pd.read_csv(join(scenario_out, 'dtw_ft_spring_all.csv'),index_col=0)

# also Yusuke would ultimately like to plot all scenarios at once
# get the field count by year to plot as a time series over the simulation period
# dtw_count = dtw_mean_all.groupby(['year','crop']).count()[['acres']].reset_index()
# dtw_count_s = dtw_mean_all_s.groupby(['year','crop']).count()[['acres']].reset_index()
# better to do the sum for the total number of acres
dtw_count = dtw_mean_all.groupby(['year','crop'])[['acres']].sum().reset_index()
dtw_count_s = dtw_mean_all_s.groupby(['year','crop'])[['acres']].sum().reset_index()
for crop in crops:
    fig,ax = plt.subplots()

    for sn in scenario_names:
        df_plt = dtw_all_scenario[dtw_all_scenario.scen_name==sn].groupby(['year','crop'])[['acres']].sum().reset_index()

        df_plt = df_plt[(df_plt.crop==crop)]
        df_plt.plot(x='year',y='acres', label=sn, ax=ax)
        # df_plt = dtw_count_s[(dtw_count_s.crop==crop)]
        # df_plt.plot(x='year',y='acres', label='Scenario', ax=ax)
        # ax.set_ylabel('Number of fields')
    ax.set_ylabel('Total Area (acres)')
    ax.set_title(crop)

    ax.set_xlabel('Year')
    plt.savefig(join(model_out, 'scenario_comparison', 'crop_acreage_by_year_'+crop+'.png'))
    plt.close()


# %%
