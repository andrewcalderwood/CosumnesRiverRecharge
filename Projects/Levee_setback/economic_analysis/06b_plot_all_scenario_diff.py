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
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# %%
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'D://WRDAPP/GWFlowModel/Cosumnes/Economic'



# %% [markdown]
# The first few plots are scenario specific, the last plots are across all listed scenarios

# %%
# update to different modflow models here
m_nam = 'input_write_2014_2025'

m_version = '' # original baseline with p_o
m_version = '_R200' # no p_o # specified later?
# need to specify base version for easier scenario switching
base_model_ws = join(loadpth, m_nam)

model_ws = join(loadpth, m_nam+m_version)

model_out = join(model_ws, 'output_clean')


# %%
scenario_summary = pd.read_excel(join(data_dir, 'scenario_summary.xlsx'))
# # subset to scenarios of interest
scenario_summary = scenario_summary[scenario_summary.base_m_nam==m_nam]
# # to identify file paths
scenario_summary['extension'] = scenario_summary.m_nam.str.extract(r'(_R\d{1,2})')
scenario_summary = scenario_summary.set_index('m_nam')


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

# %% [markdown]
# # Move to separate script
# The code below should mostly be set up to compare all three major scenarios so this should just be moved to a separate script.

# %% [markdown]
# # DTW - all scenarios

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
crops = dtw_all_scenario.crop.unique()
ny = 3
nx = int(np.ceil(len(crops)/ny))
fig,ax = plt.subplots(nx, ny, sharey=False, sharex=True, figsize=(12,6), layout='constrained', dpi=300)
        
for n,crop in enumerate(crops):
    if ny==1:
        ax_n = ax[n]
    else:
        ax_n = ax[n%nx, int(n/nx)]

    # df_plt[df_plt.crop==crop].plot(x='year',y=[val+'_m', val+'_s'], ax=ax_n, kind='bar',legend=False)
    # df_plt[df_plt.crop==crop]

    for sn in scenario_names:
        df_plt = dtw_all_scenario[dtw_all_scenario.scen_name==sn].groupby(['year','crop'])[['acres']].sum().reset_index()

        df_plt = df_plt[(df_plt.crop==crop)]
        df_plt.plot(x='year',y='acres', label=sn, ax=ax_n, legend=False)
        # df_plt = dtw_count_s[(dtw_count_s.crop==crop)]
        # df_plt.plot(x='year',y='acres', label='Scenario', ax=ax)
        # ax.set_ylabel('Number of fields')
    ax_n.set_ylabel('Total Area (acres)')
    ax_n.set_title(crop)

    ax_n.set_xlabel('Year')
ax[0,0].legend( loc='lower right')

plt.savefig(join(model_out, 'scenario_comparison', 'combined_crop_acreage_by_year.png'))
plt.close()

# %%


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
# dtw_all_scenario.groupby(['year','scen_name']).agg({'acres':'sum','UniqueID':'count'})
# dtw_all_scenario
# dtw_all_scenario.name.unique()

# %% [markdown]
# ## plot profit and yield across all scenarios

# %%
# load scenario data at once to simplify plotting 
scenarios = ['_R200','_R203','_R204']

df_econ_scenario = pd.DataFrame()
for scenario in scenarios:
    df_econ_s = pd.read_csv(join(base_model_ws+scenario,'output_clean', 'annual_profit_yield_long.csv'),index_col=0)
    df_econ_scenario = pd.concat((df_econ_scenario, df_econ_s.assign(scenario=scenario, scen_name=scenario_summary.loc[basename(base_model_ws)+scenario, 'plot_name'])))

# could improve with scenario summary sheet
scenario_names = df_econ_scenario.scen_name.unique()
df_econ_scenario = df_econ_scenario.drop(columns='index')

# %%
df_econ_scenario.to_csv(join(model_out, 'scenario_comparison', 'annual_profit_yield_crop_scenarios.csv'))

# %%
var = 'yield'
var='profit'
for val in ['total_value', 'value']:
    for var in ['yield','profit']:
        df_plt = df_econ_scenario[df_econ_scenario['var']==var].copy()
        crops = df_plt.crop.unique()
        fig,ax = plt.subplots(1, len(crops), sharey=False, figsize=(12,3), layout='constrained', dpi=300)
        
        for n,crop in enumerate(crops):
            ax_n = ax[n]
            # df_plt[df_plt.crop==crop].plot(x='year',y=[val+'_m', val+'_s'], ax=ax_n, kind='bar',legend=False)
            df_plt[df_plt.crop==crop].pivot(index='year',columns='scenario',values=val).plot(ax=ax_n, kind='bar', legend=False)
            ax_n.set_title(crop)
            ax_n.set_xlabel('Year')
            # ax_n.set_xticks(df_plt.year.unique().astype(str))
        
        
        # df_crop = df_plt[df_plt.crop==crop]
        # add legend to the last plot
        # df_plt[df_plt.crop==crop].plot(x='year',y=[val+'_m',val+'_s'], ax=ax_n, kind='bar',legend=True)
        df_plt[df_plt.crop==crop].pivot(index='year',columns='scenario',values=val).plot(ax=ax_n, kind='bar', legend=True)
        plt.legend(df_plt.scen_name.unique(), loc='lower right')
        fig.supylabel(var.capitalize())
        plt.savefig(join(model_out, 'scenario_comparison', var+'_field_avg_'+val+'.png'), bbox_inches='tight')
        plt.close()

# %% [markdown]
# ### plot total profit across all fields
# Cameron said this is a good check to verify that even if the scenarios reduce profit in some crops that overall the region has higher revenue.
# Need to plot total revenue, not profit.

# %%
df_econ_sum = df_econ_scenario.groupby(['year','var','scen_name','scenario'])[['total_value']].sum().reset_index()

# %%
# need to correct ordering for plot legend
df_econ_sum = df_econ_sum.sort_values('scenario')

df_econ_sum.to_csv(join(model_out, 'scenario_comparison', 'annual_profit_yield_summed_scenarios.csv'))

# %%
for val in ['total_value']:
    for var in ['profit']: # 'yield', doesn't work because of different units
        df_plt = df_econ_sum[df_econ_sum['var']==var].copy()
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(12,3), layout='constrained', dpi=300)
        
        ax_n = ax
        df_plt.pivot(index='year',columns='scenario',values=val).plot(ax=ax_n, kind='bar', legend=False)
        ax_n.set_title('Total across crops')
        ax_n.set_xlabel('Year')
        ax_n.ticklabel_format(style='plain', axis='y', useOffset=False)

        plt.legend(df_plt.scen_name.unique(), loc='lower left')

        fig.supylabel(var.capitalize())
        plt.savefig(join(model_out, 'scenario_comparison', var+'_summed_across_crops_'+val+'.png'), bbox_inches='tight')
        plt.close()

# %% [markdown]
# # Compare applied water

# %%
# load scenario data at once to simplify plotting 
scenarios = ['_R200','_R203','_R204']

df_wb_scenario = pd.DataFrame()
for scenario in scenarios:
    df_wb_s = pd.read_csv(join(base_model_ws+scenario,'output_clean', 'daily_WB_long.csv'),index_col=0)
    df_wb_scenario = pd.concat((df_wb_scenario, df_wb_s.assign(scenario=scenario, scen_name=scenario_summary.loc[basename(base_model_ws)+scenario, 'plot_name'])))

# could improve with scenario summary sheet
scenario_names = df_wb_scenario.scen_name.unique()
df_wb_scenario.date = pd.to_datetime(df_wb_scenario.date)
# df_wb_scenario = df_wb_scenario.drop(columns='index')

# %%
df_wb_scenario['year'] = df_wb_scenario.date.dt.year
wb_ann_sum = df_wb_scenario.groupby(['name','var','scenario', 'scen_name', 'year'])[['value','total_value']].sum()
wb_ann_sum = wb_ann_sum.reset_index()
plt_var = 'GW_applied_water'
wb_ann_sum = wb_ann_sum[wb_ann_sum['var']==plt_var]

# %%
wb_dir = join(model_out, 'scenario_comparison', 'wb_comp')
os.makedirs(wb_dir, exist_ok=True)

# %%
crop_names = wb_ann_sum.name.unique()
units = ['(meters)', '(cu. m.)']
for nc, col in enumerate(['value','total_value']):
    for crop in crop_names:
        fig,ax = plt.subplots()
    
        for sn in scenario_names:
            df_plt = wb_ann_sum[wb_ann_sum.scen_name==sn]
    
            df_plt = df_plt[(df_plt.name==crop)]
            df_plt.plot(x='year',y=col, label=sn, ax=ax)
            # df_plt = dtw_count_s[(dtw_count_s.crop==crop)]
            # df_plt.plot(x='year',y='acres', label='Scenario', ax=ax)
            # ax.set_ylabel('Number of fields')
        ax.set_ylabel(plt_var.replace('_',' ')+units[nc])
        ax.set_title(crop)
    
        ax.set_xlabel('Year')
        plt.savefig(join(wb_dir, plt_var+'_by_year_'+col+'_'+crop+'.png'))
        plt.close()

# %% [markdown]
# ## econ review
#

# %%
