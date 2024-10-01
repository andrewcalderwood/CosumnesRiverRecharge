# ---
# jupyter:
#   jupytext:
#     formats: py:percent
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

# %%
# standard python utilities
import os
from os.path import basename, dirname, join, exists
import sys
from importlib import reload
import glob
import pandas as pd
import numpy as np
import calendar
import time

# standard python plotting utilities
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# standard geospatial python utilities
# import pyproj # for converting proj4string
# import shapely
import geopandas as gpd
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator

# import flopy
# import flopy.utils.binaryfile as bf
from importlib import reload
import h5py



# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

bas_dir = join(gwfm_dir, 'BAS6')
proj_dir = join(gwfm_dir, 'Projects','EcoFIP')
plt_dir = join(proj_dir,'figures/')
out_dir = join(proj_dir,'output')
gis_dir = join(proj_dir, 'GIS')


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

from map_cln import gdf_bnds, plt_cln, xy_lab
from report_cln import base_round

from mf_utility import get_dates, clean_hob


# %%
# method to set default plot parameters, no longer need to specify 300 dpi each time, may not need to specify dimensions either
plt.rcParams.update({"figure.dpi": 300})

# %% [markdown]
# # Load data
#
# 1. The quick method to estimate recharge rates for the areas is to do a nearest spatial join with the sfr grid cells 
# 2. While we only have stream seepage data for the existing stream channel, we could do a regression like Steve did to extrapolate to nearby cells what stream seepage rates we might expect based on the streambed conductance and/or stream depth. 
#

# %% [markdown]
# ## spatial reference

# %%
gdf_sfr = gpd.read_file(join(gis_dir, 'sfr_reach_reference.shp'))

# %%
grid_id = 'parcel'
grid_id = 'sq'
if grid_id=='parcel':
    gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'parcels_elevation.shp'))
elif grid_id == 'sq':
    gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'sq_10ac_elevation.shp'))

# simpler geodataframe to bring dataframe to geodataframe
gdf_id = gdf_elev[['Region','geometry']].copy()

# %%
gdf_elev.shape
# gdf_elev.Region.unique().shape

# %%
season='spring'
gdf_gw_long = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 
                         'GW_elevations_long_'+grid_id+'_'+season+'.shp'))
gdf_gw_long.year = gdf_gw_long.year.astype(int)

# %%
# gdf_sfr.plot()
# nearest spatial join to identify the sfr grid cells to reference for output in each
# shouldn't extrapolate too far away from the stream considering the mean lengths in the transverse are at 
# most 400 m
# the issue with joining by parcel is that is selects only one reach for a long parcel
# sfr_ref = gpd.sjoin_nearest(gdf_id, gdf_sfr, max_distance=500)
# find the nearest polygon to each stream reach (want to map each stream reach to a parcel)
# to make sure parcels aren't overly weighted by one stream reach
sfr_ref = gpd.sjoin_nearest(gdf_sfr, gdf_id, max_distance=500)
# but use the grid shapefile for the output 
sfr_ref = gdf_id.merge(sfr_ref.drop(columns=['geometry','index_right']))

for n in np.arange(0,4):
    # need to repeat the process to reach the polygons that weren't identified with stream reaches
    # identify grid polygons that weren't joined to a stream reach
    gdf_id_missing = gdf_id[~gdf_id.Region.isin(sfr_ref.Region.unique())].copy()
    sfr_ref2 = gpd.sjoin_nearest(gdf_sfr, gdf_id_missing, max_distance=500)
    # # but use the grid shapefile for the output 
    sfr_ref2 = gdf_id_missing.merge(sfr_ref2.drop(columns=['geometry','index_right']))
    sfr_ref = pd.concat((sfr_ref, sfr_ref2))

# %%
# plot shows that the grid polyogn aligns with the stream reach
# sfr_ref.plot('rch_order', legend=True)

# %%
gdf_sfr_plt = gdf_sfr.copy()
# manually got to 200 m round to have one of 10k to 50k by 5
gdf_sfr_plt['dist_m_round'] = base_round(gdf_sfr_plt.dist_m, 2E2)
sfr_dist_labels = gdf_sfr_plt[gdf_sfr_plt.dist_m_round.isin(np.arange(0, 55E3,5E3))]



# %%
# # xy_lab?
from map_cln import plt_arrow, make_multi_scale

# %%
fig,ax_n = plt.subplots(figsize=(6.5,5.5),dpi=300)

gdf_sfr.plot(ax=ax_n)
for n,d in enumerate(sfr_dist_labels.dist_m_round):
    d_lab = str(int(d/1E3))
    xy = list(sfr_dist_labels.iloc[n].geometry.centroid.coords)[0]
    # ax_n.annotate(d_lab, xy)
    xy_lab(xy, d_lab, ax=ax_n, fontsize=8)
    
ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
# gdf_bnds(gdf_id, ax=ax_n, buf=2E3)

ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
                crs='epsg:26910', attribution=False)

plt_arrow(ax_n, 0.925, 0.2)
make_multi_scale(ax_n, 0.65, 0.1, dist = 2E3, scales = [4,2,1], units='km')


# %%
(1/(30*180))*86400

# %% [markdown]
# ## time series

# %%
id_cols = ['Total distance (m)','segment','reach','rch_order']
val_cols = ['Qin','Qaquifer', 'Qout', 'Qaquifer_rate','width','depth']

# %%
# monthly means
sfrdf_mon = pd.read_hdf(join(out_dir, 'sfrdf_mon_all.hdf5'))

# %%
## it's helpful for the group to have a spatial reference to this
r = 5
# r = best10.realization.iloc[1]
# subset to one realization for spatial plots
sfr_r = sfrdf_mon[sfrdf_mon.realization==r].copy()
# take the simulation average for plotting (could also focus on wet-season)
sfr_r_avg = sfr_r.groupby(['Total distance (m)']).mean(numeric_only=True)

#

# %% [markdown]
# ### Load individual realization to compare updates against others

# %%
run_dir = 'C://WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'
model_nam = 'fine_tprogs'
base_model_ws = join(loadpth, model_nam)
load_only = ['DIS','BAS6','UPW','OC','SFR','LAK',
            ]
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only
                              )
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')
# round of the steady state period
dt_ref['dt'] = dt_ref.dt.dt.round('D')
dt_ref = dt_ref[~dt_ref.steady]

# %%
from mf_utility import clean_sfr_df
pd_sfr = pd.DataFrame(gdf_sfr.drop(columns=['geometry']))
pd_sfr = pd_sfr.set_index(['iseg','ireach']) #rename(columns={'iseg':'segment','ireach':'reach'})
sfrdf = clean_sfr_df(base_model_ws, dt_ref, pd_sfr)
# calculate the effective rate of seepage
sfrdf['Qaquifer_rate'] = sfrdf.Qaquifer/(sfrdf.rchlen*sfrdf.width)

# %%
sfrdf_r_avg = sfrdf.groupby(['dist_m']).mean(numeric_only=True)
fig,ax = plt.subplots()
sfrdf_r_avg.plot(y='Qaquifer_rate',ax=ax,  label='fine tprogs')
sfr_r_avg.plot(y='Qaquifer_rate', ax=ax, label='r5 old')
plt.yscale('log')
plt.legend()

# %% [markdown]
# Reviewing the mean and std deviation across realizations shows that there is more variability than originally thought. It's also important to note that the arithmetic mean overly weights the high seepage rates so the median is more appropriate. This also brings into question how to do the averaging across time, but perhaps the arithmetic mean is appropriate since the rates across time typically remained in a smaller range it appeared.
# - review standard error across time for characteristic reaches (grave/sand, sandy mud, mud)

# %%
# calculate average rates across all time for each reach
sfr_rch_avg = sfrdf_mon.groupby(['Total distance (m)', 'realization']).mean(numeric_only=True).reset_index('realization')
# across the 10 realizations summarize the results
# the mean should be representative if normally distributed with low std deviation
# using the mean across realizations leads to an overestimate of the peak values
sfr_mean = sfr_rch_avg.groupby(['Total distance (m)']).mean(numeric_only=True)
sfr_std = sfr_rch_avg.groupby(['Total distance (m)']).std(numeric_only=True)
# quantiles are an alternative to the 95% CI for identifying reasonable limits
sfr_q = sfr_rch_avg.groupby(['Total distance (m)']).quantile([0.1, 0.25,0.5,  0.75, 0.9], numeric_only=True)
sfr_q = sfr_q.reset_index().rename(columns={'level_1':'quant'})


# %%
sfr_q[sfr_q.quant==0.75].Qaquifer_rate.describe()

# %%
# calculate coefficient of variation
cv = (sfr_std/sfr_mean)
# cv.plot(y='Qaquifer_rate')

# %%
# # calculate the 95% confidence intervals
# shows pretty extreme values
# sfr_CI_up = sfr_mean.copy()
# sfr_CI_up[val_cols] = sfr_mean[val_cols] + 2*sfr_std[val_cols]
# sfr_CI_dn = sfr_mean.copy()
# sfr_CI_dn[val_cols] = sfr_mean[val_cols] - 2*sfr_std[val_cols]


# %%
# join averaged data to geodataframe to map
# sfr grid cells only
# gdf_sfr_r_avg = gdf_sfr.merge(sfr_r_avg[['rch_order','Qaquifer_rate']])
# extrapolate to EcoFIP analysis units
# gdf_sfr_r_avg = sfr_ref.merge(sfr_r_avg[['rch_order','Qaquifer_rate']])

def map_sfr_leak(sfr_ref, sfr_r_avg, ax_n, vscale=None):
    # extrapolate to EcoFIP analysis units
    gdf_sfr_r_avg = sfr_ref.merge(sfr_r_avg[['rch_order','Qaquifer_rate']])
    # need to solve rates for when a sfr cell was equally close to an analysis unit
    gdf_sfr_r_avg = gdf_sfr_r_avg.dissolve(np.mean)
    
    # plot just the loss_rate to avoid negative values with log-scale
    gdf_sfr_r_avg['loss_rate'] = gdf_sfr_r_avg.Qaquifer_rate.copy()
    gdf_sfr_r_avg.loc[gdf_sfr_r_avg.Qaquifer_rate<=0, 'loss_rate'] = 0
    if vscale is None:
        vmin = gdf_sfr_r_avg.loc[gdf_sfr_r_avg.loss_rate>0, 'loss_rate'].min()
        vmax = gdf_sfr_r_avg.loss_rate.max()
    elif vscale is not None:
        vmin, vmax = vscale
    
    gdf_sfr_r_avg.plot('loss_rate', ax=ax_n, legend=True,
                      norm = mpl.colors.LogNorm(vmin=vmin, vmax = vmax),
                      legend_kwds={"shrink":0.6, 'label':'Simulation Average\nLoss Rate (m/day)',
                                  'extend':'both'})
    # SymLogNorm(linthresh=1E-3) is a way to have a symmetric log scale
    ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)
    
    gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
    # gdf_bnds(gdf_id, ax=ax_n, buf=2E3)
    
    ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
                    crs='epsg:26910', attribution=False)
    return gdf_sfr_r_avg


# %%
# # sfr grid cells only
# gdf_sfr_r_avg = gdf_sfr.merge(sfr_r_avg[['rch_order','Qaquifer_rate']])
# map_sfr_leak(gdf_sfr, sfr_r_avg)

# %%

# %%
# it turns out that

# g = sns.lineplot(sfr_rch_avg, x='Total distance (m)', y='Qaquifer_rate', errorbar=('ci',95))
# g = sns.lineplot(sfr_rch_avg, x='Total distance (m)', y='Qaquifer_rate', errorbar=('se',2))
g = sns.lineplot(sfr_rch_avg, x='Total distance (m)', y='Qaquifer_rate', errorbar=('pi',95))
g.set(yscale='log', ylim=(1E-2, 1E2), title = '10 Realizations')
g.set( ylabel='Stream loss rate (m/day)')

# %%
# plt_val = 'Qaquifer_rate'
# fig,ax_n = plt.subplots()
# sfr_mean.plot(y=plt_val, ax=ax_n, legend=False)
# sfr_CI_up.plot(y=plt_val, ax=ax_n,legend=False)
# sfr_CI_dn.plot(y=plt_val,ax=ax_n, legend=False)
# ax_n.set_yscale('log')

# %%
# columns to save for reference
out_cols = ['Region','rch_order','dist_m','Qaquifer_rate','geometry']
# gdf_sfr_r_avg[gdf_sfr_r_avg.Qaquifer_rate !=gdf_sfr_r_avg.loss_rate]

# %%
# gdf_sfr_r_avg


# %%
# vscale = (1E-2, 1E0)

# fig,ax_n = plt.subplots(figsize=(6.5,5.5))
# q=50
# gdf_sfr_r_avg = map_sfr_leak(sfr_ref, sfr_q[sfr_q.quant==q/100],ax_n, vscale=vscale)
    

# %%
# sfr_q[sfr_q.quant==q/100].plot(x='Total distance (m)', y='Qaquifer_rate')

# %%
os.makedirs(join(out_dir, 'stream_loss_'+grid_id), exist_ok=True)


# %%
# 1538
# map_sfr_leak(gdf_sfr_r_avg)
vscale = (1E-2, 1E1) # better for mean which has high values
vscale = (1E-2, 1E0)
# map_sfr_leak(gdf_sfr_r_avg, vscale=vscale)
# map_sfr_leak(sfr_ref, sfr_r_avg, vscale=vscale)
# map_sfr_leak(sfr_ref, sfr_mean, vscale=vscale)
# fig, ax = plt.subplots(3,1,figsize=(6.5,10))
for n, q in enumerate([25, 50, 75]):
    # ax_n = ax[n]
    fig,ax_n = plt.subplots(figsize=(6.5,5.5))
    gdf_sfr_r_avg = map_sfr_leak(sfr_ref, sfr_q[sfr_q.quant==q/100],ax_n, vscale=vscale)
    # save quantiles to shapefile
    sfr_out_file = gdf_sfr_r_avg[out_cols].rename(columns={'Qaquifer_rate':'Qaq_rate'})
    sfr_out_file.to_file(join(out_dir, 'stream_loss_'+grid_id, 'stream_loss_p'+str(q)+'.shp'))
    plt.savefig(join(out_dir, 'stream_loss_'+grid_id, 'stream_loss_p'+str(q)+'.png'))
    plt.close()
# fig.tight_layout()
# gdf_sfr_r_avg = map_sfr_leak(sfr_ref, sfr_q[sfr_q.quant==0.5], vscale=vscale)
# gdf_sfr_r_avg = map_sfr_leak(sfr_ref, sfr_q[sfr_q.quant==0.75], vscale=vscale)


# %%
# sfr_out_file

# %% [markdown]
# The stream loss map shows even less contrast without log scale, and it helps to limit the log scale value extent to avoid the middle values looking the same.  
# The quantiles 0.25 and 0.75 show a much lower peak value than the mean which indicates high rates are overwhelming the lower values. 
# - This means the median or a geometric mean should be used.
# - the quantile values of 0.25, 0.5, 0.75 are still helpful as they show hot spots with more conservative values

# %%
# temporary code to fix the saved file format
# sfrdf_mon.to_hdf(join(out_dir, 'sfrdf_mon_all.hdf5'), 
#                                   key='all', complevel=4, index=False,
#                  data_columns = sfrdf_mon.columns, format='table', mode='w')

# %%
# # testing new hdf5 format
# with h5py.File(join(out_dir, 'sfrdf_mon_all.hdf5')) as f:
#     grp = f['monthly']
#     tbl = grp['table'][:]
#     print(grp.keys())
# # testing
# pd.DataFrame(tbl)

# %%
sfrdf_all = pd.read_hdf(join(out_dir, 'sfrdf_all.hdf5'))

# %%

sfrdf_all['wy'] = sfrdf_all.dt.dt.year
sfrdf_all.loc[sfrdf_all.dt.dt.month>=10,'wy']+=1

# %%
sfr_chk = sfrdf_all.copy()
# subset to upper reaches only to avoid excess recharge predicted in lower reaches
sfr_chk = sfr_chk[sfr_chk['Total distance (m)'] <40E3]
wy_mean = sfr_chk.groupby(['wy','realization']).mean(numeric_only=True)
wy_sum = sfr_chk.groupby(['wy','realization']).sum(numeric_only=True)

# %%
wy_Qaq_rate = wy_sum.Qaquifer/(wy_mean.width*sfr_ref.dist_m.max())
wy_Qaq_rate.groupby('wy').mean()
# wy_Qaq_rate.mean()

# %%
8E7/(52*(1E3)*30)
202/(0.3048**3)
