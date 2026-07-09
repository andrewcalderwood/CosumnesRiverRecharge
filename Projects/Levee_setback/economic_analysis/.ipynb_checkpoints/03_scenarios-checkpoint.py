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
# Script to write out scenarios to run under 03_copy_model_modflow and then 03_model_connect.py

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename
import sys
import glob
from importlib import reload
import shutil

import pandas as pd
import numpy as np
from scipy.stats import hmean, gmean

# import calendar
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import geopandas as gpd
# import rasterio

from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
gwfm_dir


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

from mf_utility import get_layer_from_elev
from map_cln import gdf_bnds



# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
run_dir = 'C:/WRDAPP/GWFlowModel'
run_dir = 'F://WRDAPP/GWFlowModel'
run_dir = 'D://WRDAPP/GWFlowModel'

# loadpth = run_dir +'/Cosumnes/levee_setback/streamflow/'
# # model_nam = 'setback_streamflow'
# model_nam = 'historical_streamflow'

loadpth = run_dir +'/Cosumnes/Regional/'

model_nam = 'historical_simple_geology_reconnection'
model_nam = 'input_write_2014_2020'

base_model_ws = loadpth+model_nam


# %%
load_only=['DIS','UPW','BAS6']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only=load_only,
                              )


nrow,ncol,nlay,delr,delc = (m.dis.nrow, m.dis.ncol, m.dis.nlay, m.dis.delr, m.dis.delc)

if 'LPF' in m.get_package_list():
    gel_nam = 'LPF'
else:
    gel_nam = 'UPW'
gel = m.__getattr__(gel_nam)

# %%
# loadpth = run_dir +'/Cosumnes/Economic/'

# m.model_ws = join(loadpth, model_nam, 'crop_modflow')
# # drop HOB since we don't want to update it
# m.remove_package('HOB')
# # re-write name file before copying
# # causes issue with file path for tab file
# m.write_name_file()

# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')


# %%
# dem_data = np.copy(m.dis.top.array)
# botm = np.copy(m.dis.botm.array)

# %%
# For the tab files the left column is time (in model units) and the right column is flow (model units)
# Time is days, flow is cubic meters per day
# USGS presents flow in cfs (cubic feet per second)
inflow = pd.read_csv(join(gwfm_dir, 'SFR_data', 'MB_daily_flow_cfs.csv'), index_col = 'datetime', parse_dates = True)
# covnert flow from cubic feet per second to cubic meters per day
inflow['flow_cmd'] = inflow.flow_cfs * (86400/(3.28**3))
inflow = inflow.dropna(how='all', axis=1)

# save tabfiles dict to variable for reference
# tabfiles_dict = m.sfr.tabfiles_dict
print('Inflow MB data', inflow.index.min(), inflow.index.max())

# %%
# deer creek doesn't flow in dry-season
# have minimum flow of 100/200 cfs to start deer creek flowing
# deer creek is approximated as about 10% of Cosumnes flow during wet season
dc_min_flow = 100*(86400/(3.28**3))
inflow_dc = inflow.copy()
# anything above the flow threshold is 10% of Cosumnes
inflow_dc.flow_cmd *= 0.1
# anything below the flow threshold is 0
inflow_dc.loc[inflow_dc.flow_cmd<dc_min_flow,'flow_cmd'] = 0

# %%
all_strt_date = pd.to_datetime(m.dis.start_datetime)
all_dates = all_strt_date + (m.dis.perlen.array.cumsum()-1).astype('timedelta64[D]')
all_end_date = all_dates[-1]
print(all_strt_date, all_end_date)
months = pd.date_range(all_strt_date, all_end_date, freq='MS')
years = pd.date_range(all_strt_date, all_end_date, freq='YS').year.values
# the model dates aren't actually needed here as we use the maximum range of streamflow data

# %% [markdown]
#
# ## update scenarios
# Create different recharge scenarios

# %%
div_flow = 300
mar = inflow.copy()
# the max permitted diversion rate is 15.6 cfs
div_rate_cfs = 15.6
print('Total estimated annual diversion %.1f' %(div_rate_cfs*86400/43560*90), 'AF at max diversion')

mar['rch_cfs'] = 0
mar.loc[mar.flow_cfs>div_flow, 'rch_cfs'] = div_rate_cfs
# diversion should only be Dec 15 to Mar 15 or more simply Jan - Mar
mar.loc[~mar.index.month.isin(np.arange(1,4)),'rch_cfs'] = 0
# calculate recharge recharge rate in m3/day
mar['rch'] = mar.rch_cfs*86400*0.3048**3
# need to assign recharge to cells

# %% [markdown]
# After creating a scenario you need to:
# 1. Make sure new model workspace and modflow run files are available
# 2. Run model_connect and make sure everything runs

# %% [markdown]
# ## vineyard recharge

# %%
# load kautz vineyards to assign to cells
teichert = gpd.read_file(join(gwfm_dir,'Mapping','Kautz_shapefiles', 'Kautz Property.shp'))
rooney = gpd.read_file(join(gwfm_dir,'Mapping','Kautz_shapefiles', 'Kautz Property2.shp'))
vineyards = pd.concat((teichert, rooney)).to_crs(grid_p.crs)
vineyards_grid = gpd.overlay(grid_p, vineyards).drop(columns=['Id','Name','Area'])
vineyards_grid = vineyards_grid[vineyards_grid.geometry.area>(200*200*0.5)]


# %%
# scale the total recharge flux by the area to get the average rate
mar['rch_rate'] = mar.rch/vineyards_grid.geometry.area.sum()
# convert mar to grid level
# need date as index, row, column, rch_rate (m/day)
mar_grid = mar[['rch_rate']].assign(id=0).reset_index().merge(vineyards_grid.assign(id=0)).drop(columns=['id'])
mar_grid_out = mar_grid[['datetime','row','column','rch_rate']].rename(columns={'datetime':'date'})
mar_grid_out.to_csv(join(proj_dir, 'scenarios', 'R1_MAR_max_diversion_for_available_flow.csv'),index=False)

# %%
# increase the rch_rate by the new maximum proposed rate, approx 3x greater (~38 cfs)
mar_grid_out_max = mar_grid_out.copy()
mar_grid_out_max.rch_rate *= (6000/2444)
mar_grid_out_max.to_csv(join(proj_dir, 'scenarios', 'R2_MAR_3x_diversion_for_available_flow.csv'),index=False)

# increase the rch_rate by the new maximum proposed rate, to 6x greater (~93.6 cfs
mar_grid_out_max = mar_grid_out.copy()
mar_grid_out_max.rch_rate *= (6)
mar_grid_out_max.to_csv(join(proj_dir, 'scenarios', 'R3_MAR_6x_diversion_for_available_flow.csv'),index=False)


# %%
def plt_rch(mar):
    fig,ax = plt.subplots(3,1, figsize=(8,8), sharex=True, dpi=300)
    mar.plot(y='flow_cfs', ax=ax[0], legend=False)
    ax[0].set_ylabel('Streamflow (cfs)')
    ax[0].set_yscale('log')
    mar.plot(y='rch_cfs', ax=ax[1], legend=False)
    ax[1].set_ylabel('Diversion \nrate (cfs)')
    
    mar.resample('AS-Oct')['rch_rate'].sum().plot(ax=ax[2])
    ax[2].set_ylabel('Recharge annual \nsum (m)')
    
    ax3a = ax[2].twinx()  # instantiate a second Axes that shares the same x-axis
    mar.resample('AS-Oct')['rch_cfs'].sum().multiply(86400/43560).plot(ax=ax3a)
    ax3a.set_ylabel('acre-feet')
    
    plt.xlim('2014-10-1','2025-9-30')
    plt.xlabel('Date')
    # mar


# %%
plt_rch(mar)

# %% [markdown]
# ## floodplain recharge

# %%
deer_ck = gpd.read_file(join(gwfm_dir, 'SFR_data','final_grid_sfr','Deer_Creek_sfr.shp'))
deer_ck['name'] = 'Deer Creek'
cosumnes_r = gpd.read_file(join(gwfm_dir, 'SFR_data','final_grid_sfr','Cosumnes_River_sfr.shp'))
cosumnes_r['name'] = 'Cosumnes River'

river = pd.concat((deer_ck, cosumnes_r))

# %%
buf=1E3
deer_ck_buf = deer_ck.copy()
deer_ck_buf.geometry = deer_ck.buffer(buf)

cr_buf = cosumnes_r.copy()
cr_buf.geometry = cr_buf.buffer(buf) #original

# original code used simple intersect which included cells in Wilton
floodplain = gpd.sjoin(cr_buf[['reach','geometry']], deer_ck_buf[['geometry']], predicate='intersects')

# simplify to one geometry
floodplain = gpd.GeoDataFrame([0], geometry = [floodplain.unary_union], crs=floodplain.crs)
floodplain_grid = gpd.overlay(grid_p, floodplain)
floodplain_grid = floodplain_grid.drop(columns=0)

# %% [markdown]
# Original method used 90th percentile across all time, I should verify but I believe this should be done on a daily historical basis as State Board says to use USGS website which does on a daily basis. So there is a higher flow requirement in winter which will mean less recharge. 90th percentile ranges from 40 in late summer to a few thousand in January.
#
#
# 1. Check the gage to determine if flows are above the 90th percentile at the point of diversion (POD). 2. If flows are above the 90th percentile, calculate their diversion rate as the lesser of: a. [Diversion rate] = [Actual flows] minus [the 90th percentile flow for that day] or b. [Diversion rate] = 0.2 multiply [Actual flows] 

# %%
# for floodplain
fp_mar = mar.copy()
# when flows exceed 5E3 cfs allow 90% flow then allow 20% diversion
flow_90 = fp_mar.flow_cfs.quantile([0.9]).values[0]
print('The 90th percentile is %.2e cfs' %flow_90)
fp_mar.loc[fp_mar.flow_cfs>flow_90, 'rch_cfs'] = fp_mar.loc[fp_mar.flow_cfs>flow_90, 'flow_cfs']*0.2
# calculate recharge recharge rate in m3/day
fp_mar['rch'] = fp_mar.rch_cfs*86400*0.3048**3

# scale the total recharge flux by the area to get the average rate
fp_mar['rch_rate'] = fp_mar.rch/floodplain_grid.geometry.area.sum()

fp_mar_old = fp_mar.copy()

# %%
# outdated?? this is overwritten below anyway
code below does the proper diversion criteria on daily basis

# convert mar to grid level
# need date as index, row, column, rch_rate (m/day)
fp_mar_grid = fp_mar[['rch_rate']].assign(id=0).reset_index().merge(floodplain_grid.assign(id=0)).drop(columns=['id'])
fp_mar_grid_out = fp_mar_grid[['datetime','row','column','rch_rate']].rename(columns={'datetime':'date'})
# fp_mar_grid_out.to_csv(join(proj_dir, 'scenarios', 'R4_floodplain_MAR_90_20_diversion.csv'),index=False)



# %%
# floodplain.drop(columns=0).to_file(join(proj_dir, 'scenarios', 'R4_floodplain_MAR_outline.shp'))

# %% [markdown]
# ## second floodplain scenario - more refined
# The actuall scenario should cover less of the area and follow the proper diversion criteria daily 90th percentile
#
# Also diversion should only technically be from December to March (when most water is around anyway)

# %%
# better approach would be to do an overlay
floodplain_clean = gpd.overlay(cr_buf[['reach','geometry']], deer_ck_buf[['geometry']])
# simplify to one geometry
floodplain_clean = gpd.GeoDataFrame([0], geometry = [floodplain_clean.unary_union], crs=floodplain_clean.crs)
floodplain_grid = gpd.overlay(grid_p, floodplain_clean)
floodplain_grid = floodplain_grid.drop(columns=0)

# %%
# for floodplain
fp_mar = mar.copy()
# when flows exceed 5E3 cfs allow 90% flow then allow 20% diversion
fp_mar['month'] = fp_mar.index.month
fp_mar['day'] = fp_mar.index.day
# the 90th percentile is calculated on a daily basis
flow_90 = fp_mar.groupby(['month','day'])[['flow_cfs']].quantile(0.9).reset_index().rename(columns={'flow_cfs':'flow_90'})
# print('The 90th percentile is %.2e cfs' %flow_90)
# add in 90th percentile flows to verify conditions
fp_mar = fp_mar.reset_index().merge(flow_90)
# calculate the recharge as 20% of daily flow when exceeding the 90th daily percentile
fp_mar.loc[fp_mar.flow_cfs>fp_mar.flow_90, 'rch_cfs'] = fp_mar.loc[fp_mar.flow_cfs>fp_mar.flow_90, 'flow_cfs']*0.2
# # calculate recharge recharge rate in m3/day
fp_mar['rch'] = fp_mar.rch_cfs*86400*0.3048**3

# # # scale the total recharge flux by the area to get the average rate
fp_mar['rch_rate'] = fp_mar.rch/floodplain_grid.geometry.area.sum()

# %%

# convert mar to grid level
# need date as index, row, column, rch_rate (m/day)
fp_mar_grid = fp_mar[['datetime','rch_rate']].assign(id=0).reset_index().merge(floodplain_grid.assign(id=0)).drop(columns=['id'])
fp_mar_grid_out = fp_mar_grid[['datetime','row','column','rch_rate']].rename(columns={'datetime':'date'})
fp_mar_grid_out.to_csv(join(proj_dir, 'scenarios', 'R4_floodplain_MAR_90_20_diversion.csv'),index=False)



# %%
floodplain_clean.drop(columns=0).to_file(join(proj_dir, 'scenarios', 'R5_floodplain_MAR_clean_outline.shp'))

# %%
# map review shows the floodplain scenario creeps into Wilton a bit so it is 
# a little impractical, would make more sense to force a right buffer on the Cosumnes River 
import contextily as ctx
fig,ax = plt.subplots(figsize=(8,8))
floodplain.plot(ax=ax, color='None')
floodplain_clean.plot(ax=ax, color='None', edgecolor='blue')
# deer_ck_buf.plot(ax=ax, color='None', edgecolor='green')
# cr_buf.plot(ax=ax, color='None', edgecolor='blue')
ctx.add_basemap(ax, source = ctx.providers.Esri.WorldImagery, crs='epsg:26910', alpha = 0.8, attribution=False)


# %%
fig,ax = plt.subplots()
fp_mar.set_index('datetime').plot(y='rch_rate', ax=ax)
fp_mar_old.plot(y='rch_rate', ax=ax,linestyle='--')
ax.set_xlim('2010-1-1','2026-1-1')
# not big differences between the approaches

# %% [markdown]
# ## plot climate

# %%
strt_date = '2000-10-1'
end_date = '2022-9-30'

# %%
uzf_dir = join(gwfm_dir, 'UZF_data')
## Potential ETo spatial interpolation from CIMIS
fn = glob.glob(join(uzf_dir,'CIMIS', 'Cosumnes_dailyET_precip*.csv'))
daily_data = pd.DataFrame()
for file in fn:
    new_data = pd.read_csv(file, index_col = ['Date'], parse_dates = True)
    daily_data = pd.concat((daily_data, new_data))
# units of mm
data_in = daily_data[daily_data['Stn Name']=='Fair Oaks']
# clean up data so columns are by location, units of Precip are in mm
rain_in = data_in.pivot_table(index = 'Date', columns = 'Stn Name', values = 'Precip (mm)')
rain_m = rain_in/1000
# subset to model and reindex to match streamflow dates
rain_plt = rain_m['Fair Oaks'].loc[strt_date:end_date]
rain_plt = rain_plt.reindex(inflow.index)



# %%
data_in.resample('AS-Oct').sum(numeric_only=True)

# %%
rain_plt_ft = rain_plt/0.3048

# %%
flw_lgd =  [
    # Patch(facecolor='none', edgecolor='black', linestyle = '--', alpha=1., label='Watershed Extent'),
    Line2D([0], [0],color='tab:blue', alpha=0.6, label='Rainfall'),
    Line2D([0], [0], color='brown', alpha=0.6,label='Streamflow'),
    # Line2D([0], [0], linestyle='-.', color='black', label='23 cms'),
    # Line2D([0], [0],  linestyle='--', color='black', label='71.6 cms'),
]

# %%
from map_cln import gdf_bnds, pnt_2_tup, lab_pnt, plt_cln, arr_lab, xy_lab, make_multi_scale, dir_arrow, plt_arrow


# %%
fig, ax = plt.subplots(figsize=(10,5))

dt = inflow.index.values
ax.plot(dt, inflow['flow_cfs'], color='brown', alpha=0.6)
# Create second axes, in order to get the bars from the top you can multiply by -1
ax2 = ax.twinx()
# ax2.bar(dt, -rain_plt.values, 0.9)
ax2.plot(dt, -rain_plt_ft, alpha=0.6)# regular line plot does better than a bar plot
# Now need to fix the axis labels
max_pre = np.max(rain_plt_ft)+0.01
# max_pre = 0.12
max_pre = 0.12/0.3048
# y2_ticks = np.arange(0, max_pre, 0.02)
y2_ticks = np.arange(0, max_pre, 0.06)
y2_ticklabels = [str(i) for i in y2_ticks]
ax2.set_yticks(-1 * y2_ticks)
ax2.set_yticklabels(y2_ticklabels)

ax2.set_ylabel('Rainfall (ft)')
ax.set_ylabel('Streamflow ($ft^3/s$)')
# ax.ticklabel_format(style='plain', axis='y') 
ax.set_yscale('log')
ax.set_ylim(1,1E5)
ax.set_xlabel('Date')

fig.legend(handles=flw_lgd, loc='outside upper center',ncol=2, bbox_to_anchor=(0.5, 1.0)) # with flow thresholds (0.5, 1.05)
# arr_lab(lak_extent, 'Reconnected\nFloodplain', ax, offset = (-100, 150), fontsize=fontsize)
# def flow_threshold_lines(ax):
#     ax.axhline(y=23, color='black', linestyle='--')
#     ax.axhline(y=71.6, color='black', linestyle='-.')
#     xy_lab((pd.to_datetime('2015-8-1'),23),'23 $m^3/s$', offset = (0,5), fontsize=10, bbox=False)
#     xy_lab((pd.to_datetime('2015-8-1'),71.6),'71.6 $m^3/s$', offset = (0,5), fontsize=10, bbox=False)
# flow_threshold_lines(ax=ax)

# fig.savefig(join(fig_dir, 'streamflow_hydrograph_with_rain.png'), bbox_inches='tight')
# plt.show()


# %% [markdown]
# # precip only

# %%
wyt_sac = pd.read_csv(join(gwfm_dir, 'GHB_data', 'sacramento_WY_types.txt'))
color_dict = {'C':'tab:purple','D':'tab:red','BN':'tab:orange','AN':'yellow','W':'tab:green'}
name_dict = {'C':'Critical','D':'Dry', 'BN':'Below Normal', 'AN':'Above Normal','W':'Wet'}
wyt_sac['color'] = [color_dict[yt] for yt in wyt_sac['Yr-type']]
wyt_sac['name'] = [name_dict[yt] for yt in wyt_sac['Yr-type']]

# %%

# %%
# wyt for plotting
# wy = dt_ref.wy.unique()
wy = inflow.index.year.unique()[1:]
wyt = wyt_sac[wyt_sac.WY.isin(wy)].copy()
wyt['plt_date'] = wyt.WY.astype(str)+'-1-1'
wyt = wyt.set_index('WY')

# %%
annual_rain = rain_m['Fair Oaks'].resample('AS-Oct').sum()
# look at data from last 20 years
annual_rain = annual_rain[annual_rain.index.year>annual_rain.index.year.max()-20]
# # summarize typical rainfall and devication
print('Avg rain %.3f m' %annual_rain.mean())
annual_rain[annual_rain.index.year.isin([2016,2018])] - annual_rain.mean()

# %%
np.where(rain_mon.index==wyt.loc[n, 'plt_date'])
n
# wyt.loc[2001]
rain_mon.index
 # rain_mon.quantile(0.95))
# wyt
wyt.loc[n, 'plt_date']

# %%
plt_strt = '2000-10-1'
plt_end = '2022-9-30'
fig, ax = plt.subplots(figsize=(6.5,3), dpi=300)
rain_mon = rain_plt.resample('MS').sum()
rain_mon = rain_plt.resample('AS-Oct').sum()
rain_mon = rain_mon.loc[plt_strt:plt_end]

ax.bar(np.arange(0, len(rain_mon)), rain_mon.values)
# ax.set_xticks(np.arange(0, len(rain_mon))[3::12], rain_mon.index[3::12].year);
ax.set_xticks(np.arange(0, len(rain_mon))[::3], rain_mon.index[::3].year+1);
for n in wy:
    xy_lab((wyt.loc[n, 'plt_date'], rain_mon.quantile(0.95)),
    # xy_lab((np.where(rain_mon.index==wyt.loc[n, 'plt_date'])[0][0], rain_mon.quantile(0.95)),
           wyt.loc[n, 'name'].replace(' ','\n'), ax=ax, offset = (0,2), fontsize=10, bbox=True)

ax.set_ylabel('Annual Rainfall Total (m)');
# ax.set_ylabel('Monthly Rainfall Total (m)');
# ax.set_xlabel('Date');
ax.set_xlabel('WY');
# plt.savefig(join(fig_dir, 'monthly_rainfall.png'), bbox_inches='tight')
# plt.close()
