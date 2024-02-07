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
gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'parcels_elevation.shp'))
# gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'sq_10ac_elevation.shp'))

# simpler geodataframe to bring dataframe to geodataframe
gdf_id = gdf_elev[['Region','geometry']].copy()

# %%
grid_id= 'parcel'
season='spring'
gdf_gw_long = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 
                         'GW_elevations_long_'+grid_id+'_'+season+'.shp'))
gdf_gw_long.year = gdf_gw_long.year.astype(int)

# %%
# gdf_sfr.plot()
# nearest spatial join to identify the sfr grid cells to reference for output in each
# shouldn't extrapolate too far away from the stream considering the mean lengths in the transverse are at 
# most 400 m
sfr_ref = gpd.sjoin_nearest(gdf_id, gdf_sfr, max_distance=500)

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

# %%
# join averaged data to geodataframe to map
# sfr grid cells only
# gdf_sfr_r_avg = gdf_sfr.merge(sfr_r_avg[['rch_order','Qaquifer_rate']])
# extrapolate to EcoFIP analysis units
# gdf_sfr_r_avg = sfr_ref.merge(sfr_r_avg[['rch_order','Qaquifer_rate']])

def map_sfr_leak(gdf_sfr_r_avg):
    # plot just the loss_rate to avoid negative values with log-scale
    gdf_sfr_r_avg['loss_rate'] = gdf_sfr_r_avg.Qaquifer_rate.copy()
    gdf_sfr_r_avg.loc[gdf_sfr_r_avg.Qaquifer_rate<=0, 'loss_rate'] = 0
    vmin = gdf_sfr_r_avg.loc[gdf_sfr_r_avg.loss_rate>0, 'loss_rate'].min()
    vmax = gdf_sfr_r_avg.loss_rate.max()
    
    fig,ax_n = plt.subplots(figsize=(6.5,5.5),dpi=300)
    gdf_sfr_r_avg.plot('loss_rate', ax=ax_n, legend=True,
                      norm = mpl.colors.LogNorm(vmin=vmin, vmax = vmax),
                      legend_kwds={"shrink":0.6, 'label':'Simulation Average\nLoss Rate (m/day)'})
    # SymLogNorm(linthresh=1E-3) is a way to have a symmetric log scale
    ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)
    
    gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
    # gdf_bnds(gdf_id, ax=ax_n, buf=2E3)
    
    ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
                    crs='epsg:26910', attribution=False)
    return None


# %%
# # sfr grid cells only
# gdf_sfr_r_avg = gdf_sfr.merge(sfr_r_avg[['rch_order','Qaquifer_rate']])
# map_sfr_leak(gdf_sfr_r_avg)

# %%
# extrapolate to EcoFIP analysis units
gdf_sfr_r_avg = sfr_ref.merge(sfr_r_avg[['rch_order','Qaquifer_rate']])
# need to solve rates for when a sfr cell was equally close to an analysis unit
gdf_sfr_r_avg = gdf_sfr_r_avg.dissolve(np.mean)
# 1538
map_sfr_leak(gdf_sfr_r_avg)

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
