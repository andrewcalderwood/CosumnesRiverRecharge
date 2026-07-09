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


# %%
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'D://WRDAPP/GWFlowModel/Cosumnes/Economic'

# update to different modflow models here
m_nam = 'input_write_2014_2025'

m_version = '' # original baseline with p_o
m_version = '_R200' # no p_o # specified later?
# need to specify base version for easier scenario switching
base_model_ws = join(loadpth, m_nam)

model_ws = join(loadpth, m_nam+m_version)

model_out = join(model_ws, 'output_clean')

input_name = 'static_model_inputs_no_p_o.xlsx'


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

# %%
scenario_summary = pd.read_excel(join(data_dir, 'scenario_summary.xlsx'))
# # subset to scenarios of interest
scenario_summary = scenario_summary[scenario_summary.base_m_nam==m_nam]
# # to identify file paths
scenario_summary['extension'] = scenario_summary.m_nam.str.extract(r'(_R\d{1,3})')
# scenario_summary = scenario_summary.set_index('m_nam')

scen_join = scenario_summary[['plot_name','extension']].rename(columns={'extension':'scenario'})

# %% [markdown]
# ## load econ data on parcel level for finer review
# Want to evaluate relationship of profit and cost of water

# %%
# load scenario data at once to simplify plotting 
scenarios = ['','_R3','_R4']
scenarios = ['_R200','_R203','_R204']


# %%
df_wb = pd.read_csv(join(out_dir, 'daily_WB_long.csv'),parse_dates=['date'])
df_wb['year'] = df_wb.date.dt.year

# %%
# check the average budget that it meets constraints
var = 'GW_applied_water'
# df_all_plt = df_all['var']==var
# df_annual_sum = df_wb.groupby(['crop','year','name','UniqueID','var','pod']).sum(numeric_only=True).reset_index()
df_annual_sum = df_wb.groupby(['name','year','var']).sum(numeric_only=True).reset_index()


# %%
df_annual_sum[df_annual_sum['var']==var].groupby('name').value.mean()
# df_annual_sum[df_annual_sum['var']==var]

# %%
scenario = scenarios[0]
join(base_model_ws+scenario)
# not sure why I was loading the rep version?

# %%
df_all_econ = pd.DataFrame()

for scenario in scenarios:
    print(scenario)
    swb_ws = join(base_model_ws+scenario, 'crop_soilbudget')
    swb_ws_rep = join(base_model_ws+scenario, 'rep_crop_soilbudget')

    for year in run_years:
    # for year in [2015]:
        # load SWB folder
        crop_in = pd.read_csv(join(swb_ws_rep, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
        print('\n', year, end=' - ')
        # cost is the cost of water
        for var in ['profit', 'cost']:
            print(var, end=',')
            name = join(model_ws, 'crop_soilbudget', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
            with h5py.File(name) as dset:
                finished_crops = list(dset['array'].keys())
                # print(finished_crops, end='.')
            for crop in finished_crops:
            # for crop in finished_crops[2]:
                # need dates for time series water budget output
                var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop)
                # extract output and convert to dataframe with ID columns
                arr = read_crop_arr_h5(crop, name)
                df = pd.DataFrame(arr, columns=['value'])
                # add parcel information back
                df = pd.concat((df,crop_in[crop_in.name==pred_dict[crop]].reset_index(drop=True)),axis=1)
                df = df.assign(crop=crop, year=year, var=var, scenario=scenario)
                # concat to existing data
                df_all_econ = pd.concat((df_all_econ, df))

# %%
# rename as econ for plotting reference
df_all_econ = df_all_econ.rename(columns={'parcel_id':'UniqueID'}).merge(parcels[['UniqueID','acres']])
# scale value rates (1/acre) into total volumes (m3)
# df_all['total_value'] = df_all['value']*df_all.acres*4046.86

# %%
df_wide = df_all_econ.pivot_table(index=['UniqueID','name', 'year','dtw_ft','pod_bool','acres','scenario'], columns='var', values='value')
df_wide = df_wide.reset_index()
df_wide.profit *= -1
# convert cost/profit to totals
df_wide[['cost','profit']] = df_wide[['cost','profit']].apply(lambda x: x*df_wide.acres)

# %%
phi = 0.12 # KWh/acre-in/ft
p_e = 0.31 # $/KWh
p_vol = 100*phi*p_e # 100 ft * KWh/acre-in/ft * $/Kwn = $/acre-in
p_vol *= 12 # convert to $/AF
p_vol * 0.8 /0.3048 # multiply by typical grape depth of 0.8 m and convert to ft
# the code is set up so the profit and cost is per acre

# %%
# df_wide

# %%
# # start with summary plot of cost of water to confirm it aligns with pumped water
df_wide_sum = df_wide.groupby(['year','name','scenario'])[['cost','profit']].sum()
df_wide_sum = df_wide_sum.reset_index()
df_wide_sum['revenue'] = df_wide_sum.cost+df_plt.profit
df_wide_sum_scen = df_wide_sum[df_wide_sum.scenario=='_R204']
# df_wide_sum
sns.catplot(df_wide_sum_scen,x='year',y='cost', col='name',
            kind='bar', color='tab:blue',
            sharey=False
           # facet_kws={'sharey': False, 'sharex': True}
)


# %%
# plot total revenue across years
# df_wide_sum.groupby('year')[['revenue']].sum().plot.bar()

# %%
val='cost'
# val='revenue'
# df_plt = df_econ_scenario[df_econ_scenario['var']==var].copy()
df_plt = df_wide_sum.copy()
df_plt = df_plt.merge(scen_join)
crops = df_plt.name.unique()
fig,ax = plt.subplots(1, len(crops), sharey=False, figsize=(12,3), layout='constrained', dpi=300)

for n,crop in enumerate(crops):
    ax_n = ax[n]
    # df_plt[df_plt.crop==crop].plot(x='year',y=[val+'_m', val+'_s'], ax=ax_n, kind='bar',legend=False)
    df_plt[df_plt.name==crop].pivot(index='year',columns='scenario',values=val).plot(ax=ax_n, kind='bar', legend=False)
    ax_n.set_title(crop)
    ax_n.set_xlabel('Year')
    # ax_n.set_xticks(df_plt.year.unique().astype(str))


# df_crop = df_plt[df_plt.crop==crop]
# add legend to the last plot
# df_plt[df_plt.crop==crop].plot(x='year',y=[val+'_m',val+'_s'], ax=ax_n, kind='bar',legend=True)
df_plt[df_plt.name==crop].pivot(index='year',columns='scenario',values=val).plot(ax=ax_n, kind='bar', legend=True)
plt.legend(df_plt.plot_name.unique(), loc='lower right')
# plt.legend([1,2,3], loc='lower right')
fig.supylabel(val.capitalize())
# plt.savefig(join(model_out, 'scenario_comparison', var+'_field_avg_'+val+'.png'), bbox_inches='tight')
# plt.close()

# %%
df_plt.scenario.unique()
