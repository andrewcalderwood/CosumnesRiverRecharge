# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Review groundwater data seasonal trends
# Identify the typical seasonal trend for wells to create a generic depth to water profile to reduce the number
# of unique optimizations needed.

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

# # flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# # other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)
from mf_utility import get_layer_from_elev, param_load

# %%
# resampled ground surface elevation
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')

# %%
# parcel_wells = gpd.read_file(join(gwfm_dir, 'WEL_data', 'parcels_to_wells', 'parcels_to_wells.shp'))
mb_regional = gpd.read_file(join(gwfm_dir,"DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp"))


# %%
ghb_dir = join(gwfm_dir, 'GHB_data')

gwe = pd.read_csv(join(ghb_dir, 'domain_dwr_msmts.csv'))
gwe.msmt_date = pd.to_datetime(gwe.msmt_date, errors='coerce')
# assign year and season for grouping
gwe['year'] = gwe.msmt_date.dt.year
gwe['season'] = 'neither'
gwe.loc[gwe.msmt_date.dt.month.isin([3,4,5]), 'season'] = 'spring'
gwe.loc[gwe.msmt_date.dt.month.isin([9,10,11]), 'season'] = 'fall'

# drop measurements with known issues
wlm_qa_err = ['Nearby pump operating', 
       'Pumped recently', 'Pumping',
       'Recharge or surface water effects near well']
gwe = gwe[~gwe.wlm_qa_detail.isin(wlm_qa_err)]

# %%
# there is at least one well with measurements marked as questionable that show a drop from -20 to -200 ft
gwe[~gwe.wlm_qa_desc.isin(['Missing','Questionable'])]

# %%
strt_date= '2000-10-1'
gwe_new = gwe.loc[gwe.msmt_date> strt_date].copy()

gwe_chk = gwe_new.groupby(['year','site_code']).count()
# # only keep wells where 25th percentile year has at least 2 measurements
keep_sites= gwe_chk.groupby('site_code')['gwe'].quantile(0.25) >=2
keep_sites = keep_sites[keep_sites.values].index.values


# %%
# subset to useful sites
print(keep_sites.shape, gwe.site_code.unique().shape)
gwe_cln = gwe.loc[gwe.site_code.isin(keep_sites)]
gwe_cln = gwe_cln.loc[gwe_cln.msmt_date> strt_date].copy()
gwe_cln = gwe_cln.dropna(subset='gwe')

# create geodataframe of station locations
stns = gwe_cln.drop_duplicates('site_code')
stns_gdf = gpd.GeoDataFrame(stns, geometry = gpd.points_from_xy(stns.longitude, stns.latitude), crs='epsg:4326')
stns_gdf = stns_gdf.to_crs(mb_regional.crs)
# group wells by gse elevations
gse_quant = stns_gdf.wlm_gse.quantile([0, 0.25, 0.5, 0.75 ,1]).values
stns_gdf['gwe_group'] = np.nan
for n, gse in enumerate(gse_quant):
    stns_gdf.loc[stns_gdf.wlm_gse>gse, 'gwe_group'] = n+1

# find sites within the model domain
stns_gdf = gpd.overlay(stns_gdf, mb_regional)
# subset for sites in the domain
gwe_cln = gwe_cln[gwe_cln.site_code.isin(stns_gdf.site_code)]
# add elevation group to gwe data
gwe_cln = gwe_cln.merge(stns_gdf[['site_code','gwe_group']])


# %%
# stns_gdf.plot('gwe_group', legend=True)

# %%
# pattern is consistent between fall and spring
gwe_cln.hist(['year'], bins=len(gwe_cln.year.unique())) #.groupby('season')



# %%
# 
wyt_sac = pd.read_csv(join(ghb_dir, 'sacramento_WY_types.txt'))
wyt_sac['dry'] = 'tab:blue'
wyt_sac.loc[wyt_sac['Yr-type'].isin(['C','D','BN']), 'dry'] = 'tab:red' 
# wyt_sac
color_dict = {'C':'tab:purple','D':'tab:red','BN':'tab:orange','AN':'yellow','W':'tab:green'}
name_dict = {'C':'Critical','D':'Dry', 'BN':'Below Normal', 'AN':'Above Normal','W':'Wet'}
wyt_sac['color'] = [color_dict[yt] for yt in wyt_sac['Yr-type']]

# calculate the wyt group
wyt_sac['wy_group'] = (wyt_sac['Yr-type'] != wyt_sac['Yr-type'].shift(1)).cumsum()
wyt_sac['dry_group'] = (wyt_sac['dry'] != wyt_sac['dry'].shift(1)).cumsum()
wyt_sac['date'] = pd.to_datetime((wyt_sac.WY-1).astype(str)+'-10-1')

# %%
from matplotlib.patches import Patch
wy_lgd = []
for n in color_dict.keys():
    wy_lgd += [Patch(facecolor=color_dict[n], label=name_dict[n], alpha=0.5)]
# wy_lgd
dry_lgd = [Patch(facecolor='tab:red', label='Drier', alpha=0.7),
          Patch(facecolor='tab:blue', label='Wetter', alpha=0.7)]


# %%
def plt_wyt(wyt_sac, ax, alpha=0.5):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    for n in wyt_sac.wy_group.unique():
        wy = wyt_sac.loc[wyt_sac.wy_group==n]
        ax.fill_between([wy.date.min(), wy.date.max()+pd.DateOffset(months=11)], 
                         ylim[1], ylim[0], color=wy.color, alpha=alpha)
    ax.set_xlim(xlim[0],xlim[1])
    ax.legend(handles=wy_lgd, loc='lower right')
    
def plt_dry(wyt_sac, ax, alpha=0.5):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    for n in wyt_sac.dry_group.unique():
        wy = wyt_sac.loc[wyt_sac.dry_group==n]
        ax.fill_between([wy.date.min(), wy.date.max()+pd.DateOffset(months=11)], 
                         ylim[1], ylim[0], color=wy.dry, alpha=alpha)
    ax.set_xlim(xlim[0],xlim[1])
    # plt.legend(handles=dry_lgd, loc='lower right')


# %%
# use long-term datasets to identify trends
sns.relplot(gwe_cln, x='msmt_date', y='gwe', hue='site_code',row = 'gwe_group')


# %% [markdown]
# If we group by site code then the standard deviation is nearly equal to the mean for most. Grouping by water year type might improve this

# %%
wyt_join = wyt_sac[['WY','dry','Yr-type']].rename(columns={'WY':'year'})


# %%
chk_site = gwe_cln[gwe_cln.gwe<-200].site_code.unique()[0]
# chk_site = '385707N1211868W001'
gwe_cln[gwe_cln.site_code==chk_site].plot(x='msmt_date',y='gwe')
gwe_cln[gwe_cln.site_code==chk_site]

# %%
gwe_yr_plt = gwe_cln[gwe_cln.year==2020]
sns.relplot(gwe_yr_plt,x='msmt_date',y='gwe', hue='site_code',col = 'gwe_group', col_wrap= 2, kind='line')

# %%
# need to look at average decline for each year from spring to fall
# filtering to more recent years to avoid bad measurements in 2000-2010
gwe_plt = gwe_cln[gwe_cln.year>=2010].copy()
gwe_plt = gwe_plt.groupby(['site_code','year','season']).mean(numeric_only=True).reset_index()
gwe_plt = gwe_plt.merge(wyt_join)
# pivot to look at difference in fall and spring between each year
wy_var = 'dry'
wy_var='Yr-type'
gwe_plt_diff = gwe_plt.pivot(index=['site_code', 'year',wy_var], values='gwe', columns='season').reset_index()
gwe_plt_diff['diff'] = gwe_plt_diff.spring - gwe_plt_diff.fall
diff_stats = gwe_plt_diff.groupby(['site_code',wy_var])['diff'].agg(['count','mean','std'])
diff_stats

# %%
diff_stats.boxplot(by='Yr-type',column='mean')
plt.ylim(-10,20)
plt.title('Average GWE Decline (ft) from Spring to Fall by well')

# %%
# review soil data to confirm there are not
# significant differences to plot
field_ids = 'parcel' # 'ag'
# load cleaned soil data for ag fields
soil_path = join(uzf_dir,'clean_soil_data')
# soil data for each ag field
soil_ag_all = pd.read_csv(join(soil_path, 'soil_for_'+field_ids+'_fields.csv'), index_col=0)
soil_ag_all['Ks'] = np.where(soil_ag_all.Ksat_Low==0, soil_ag_all.Ksat/10, soil_ag_all.Ksat_Low)


# %%
# soil_ag_all[soil_ag_all.Texture=='clay']
# os.listdir(proj_dir)

# %%
soil_texture_summary = soil_ag_all.groupby(['Texture','HydGroup'])[['Ks','AWC', 'Porosity']].agg(['mean','std'])
soil_texture_summary.to_csv(join(proj_dir, 'model_inputs', 'soil_texture_summary.csv'))

# %%
# gwe_plt = gwe_riv_cln[gwe_riv_cln.msmt_date>'2000-1-1']
# npts = 40
# long_sites = gwe_plt.groupby('site_code').count()['gwe']
# long_sites = long_sites[long_sites>npts].index
# gwe_plt = gwe_plt[gwe_plt.site_code.isin(long_sites)]
# for n in np.arange(1,4):
#     # df_plt = gwe_riv_cln[gwe_riv_cln.gwe_group==n]
#     df_plt = gwe_plt[gwe_plt.gwe_group==n]
#     fig,ax = plt.subplots()
#     for s in df_plt.site_code.unique():
#         df_plt[df_plt.site_code==s].plot(x='msmt_date',y='gwe', color='black',
#                                          kind='line',ax=ax, legend=False)
#     # ax = sns.lineplot(df_plt,x='msmt_date',y='gwe',
#     #                   color='black',
#     #                   hue='site_code',
#     #                   legend=False
#     #                  )
#     plt_dry(wyt_sac, ax)
#     plt.ylabel('Groundwater Elevation (ft amsl)')
#     plt.xlabel('Date')
#     # plt.savefig(join(fig_dir, 'gw_hydrographs_near_river.png'), bbox_inches='tight')
#     plt.show()
#     plt.close()
