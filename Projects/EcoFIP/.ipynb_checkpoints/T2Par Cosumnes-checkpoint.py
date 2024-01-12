# ---
# jupyter:
#   jupytext:
#     formats: py:percent
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

# %%
# standard python utilities
import os
import sys
from os.path import basename, dirname, join, exists, expanduser
import glob
import time

import pandas as pd
import numpy as np
import numpy.ma as ma
from scipy.stats import gmean

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
# import pyproj # for converting proj4string
import geopandas as gpd
from osgeo import gdal
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study


sfr_dir = gwfm_dir+'/SFR_data/'
upw_dir = gwfm_dir+'/UPW_data/'

proj_dir = join(gwfm_dir, 'Projects', 'EcoFIP')


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

from importlib import reload

add_path(join(doc_dir,'GitHub', 'flopy'))
import flopy

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

from mf_utility import get_layer_from_elev
from report_cln import base_round



# %% [markdown]
# # Load data
# 1. Pre-processed AEM to array data
# 2. Well completion reports
# 3. Reference GIS

# %%
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Regional'
model_ws = join(loadpth, 'T2par_test')
os.listdir(model_ws)

# %%
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')

# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_elev = gpd.read_file(join(gwfm_dir,'DIS_data','grid_elevation_m_statistics.shp'))

# %%

# %%
wcr_path = join(dirname(gwfm_dir), 'Large_TPROGS_run')
# wcr_file = join(wcr_path, 'WCR_OD_99_plusAR_new_clean_elevation.csv') # only has 698 wells but covers most of domain, looks to be missing Alisha's wells
# id_col = 'WCR_No'
# files created by Alisha with final data
wcr_loc_file = join(wcr_path, 'All_locations_Large_Cosumnes_run.csv') # has all locations
wcr_file = join(wcr_path, 'TPROGS_update_full_run_12_2020_minus_NaN_and_duplicates.csv') # has all locations and depths

id_col = 'WCR Number'

# %%
wcr_df_in = pd.read_csv(wcr_file)
wcr_gdf = pd.read_csv(wcr_loc_file)

# %%
wcr_gdf_full = gpd.GeoDataFrame(wcr_gdf, geometry=gpd.points_from_xy(wcr_gdf.Easting, wcr_gdf.Northing), crs='epsg:3310')


# %%

domain_wcr =m_domain.to_crs(wcr_gdf_full.crs)
xll, yll = list(domain_wcr.geometry.values[0].exterior.coords)[0]


# %%
fig,ax = plt.subplots()
wcr_gdf_full.drop_duplicates(id_col, keep='first').plot(ax=ax)
domain_wcr.plot(ax=ax, color='None')
plt.plot([xll], [yll], marker='*', color='red')
print(xll, yll)

chk_wcr = 'WCR2005-003341'
wcr_gdf_full[wcr_gdf_full['WCR Number']==chk_wcr].plot(color='green', ax=ax)

# %% [markdown]
# Remove WCRs outside of the domain

# %%
# required extra buffer to remove well that Texture2PAR considered outside of the MODFLOW grid
domain_wcr_buf = domain_wcr.copy()
domain_wcr_buf.geometry = domain_wcr_buf.buffer(-100)
wcr_gdf = gpd.overlay(wcr_gdf_full, domain_wcr_buf)

print(wcr_gdf_full.shape[0]-wcr_gdf.shape[0], 'WCRs removed that were outside the boundary')

# %%
# wcr_gdf[wcr_gdf['WCR Number']==chk_wcr]

# %% [markdown]
# Wells were removed from the GeoDataFrame but the input dataframe wasn't reduced.

# %%
wcr_df_in[wcr_df_in['WCR Number'].isin(wcr_gdf['WCR Number'])].shape, wcr_df_in.shape


# %% [markdown]
# # t2py
# It's not clear how t2py is meant to run. I think it's expecting a csv with a set column naming.
# WellName - string well name  
# Well - integer well ID (1 to nwells)  
# Point - integer ID of the point in the well log (1 to ndepths)  
# PC - percent coarse3
#

# %% [markdown]
# `WellName	Well	Point	PC	X	Y	Zland	Depth`  
# `CVHMSJ_1856	1	1	1.0	0.5	496.5	25.0	0.701`  
# `CVHMSJ_1856	1	2	0.0	0.5	496.5	25.0	1.403`  
#

# %% [markdown]
# The tprogs input has from and to depths in feet. Texture2Par uses the depth from the to column.

# %% [markdown]
# ## Clean WCR data

# %%
def fill_wcr_df(wcr_df):
    wcr_df['d_from'] = np.round((wcr_df['Elevation (m)'] - wcr_df['Relative Elevation (m)'])/0.3048,1)
    wcr_df['d_to'] = wcr_df['d_from'] + 1
    print('Number of mismatch from depth to from:', wcr_df[(wcr_df['from'] - wcr_df['d_from'])>0].shape[0])
    # fill in the missing to and from data
    to_fill = wcr_df.to.isna()
    wcr_df.loc[to_fill, 'to'] = wcr_df.loc[to_fill, 'd_to']
    wcr_df.loc[to_fill, 'from'] = wcr_df.loc[to_fill, 'd_from']
    return wcr_df
wcr_df = fill_wcr_df(wcr_df_in)
# filter out to WCRs in the domain
wcr_df = wcr_df[wcr_df['WCR Number'].isin(wcr_gdf['WCR Number'])]

# %%
# has been split up into 0.5 m intervals which may not benefit classification
# I need to convert the split up well depth data to longer depths to simplify the input
wcr_nums = wcr_df[id_col].unique()
n=400
wcr_cln = pd.DataFrame()
t0 = time.time()

for n in np.arange(0, len(wcr_nums)):
    wcr_n = wcr_df.loc[wcr_df[id_col]==wcr_nums[n]].copy()
    for nc in wcr_n.index[:-1]: 
        # if the classification equals the next the update the from and drop the previous
        if wcr_n.loc[nc, 'classification']==wcr_n.loc[nc+1, 'classification']:
            wcr_n.loc[nc+1, 'from'] = wcr_n.loc[nc, 'from']
            wcr_n = wcr_n.drop(labels=wcr_n.loc[nc].name)
    # concat the simplified dataframe to the full dataframe
    wcr_n['Point'] = np.arange(1, len(wcr_n)+1)
    wcr_n['Well'] = n
    wcr_cln = pd.concat((wcr_cln, wcr_n))
t1 = time.time()
print('Total time %.1f min' %((t1-t0)/60))
# averages about 1-2  min

# %%
wcr_cln.to_csv(join(proj_dir, 'WCR_all_cleaned.csv'))

# %% [markdown]
# Need to decide how to convert facies classification to percent coarse.  
# The cutoffs I used for going from PC to facies was 0, 30, 40, 60 to match the Meirovtiz proportions and 0, 40, 50, 80 to match the new proportions from Alisha.  
# Mud:0.2, Sandy Mud:0.45, Sand:0.65, Gravel:0.9 (Alisha's proportions matching makes more sense intuitively).

# %% [markdown]
# ## Organize T2py input and write 

# %%
# take needed columns
wcr_out = wcr_cln.copy().reset_index()
# convert depth to meters
wcr_out['Depth'] = wcr_out.to*0.3048
pc_facies = pd.DataFrame([0.2, 0.45, 0.65, 0.9], [4,3,2,1], columns=['PC'])
wcr_out['PC'] = pc_facies.loc[wcr_out.classification].values
# subset columns for output
wcr_csv = wcr_out[['WCR Number', 'Well', 'Point', 'PC', 'Easting', 'Northing', 'Elevation (m)', 'Depth']]
# clean columns names for import
wcr_csv = wcr_csv.rename(columns={'WCR Number':'WellName','Easting':'X','Northing':'Y', 'Elevation (m)':'Zland'})
# the dash might causes issues with fortran code
wcr_csv['WellName'] =wcr_csv['WellName'].str.replace('-','_')

# %%
wcr_csv.to_csv(join(proj_dir, 'WCR_all_T2py_cleaned.csv'))

# %%
with open(join(model_ws, 'well_log.dat'), 'w') as f:
    f.write(wcr_csv.columns.str.ljust(12).str.cat(sep = ' '))
    f.write('\n')
    for n in np.arange(0, len(wcr_csv)):
    #         f.write(pdata_zone.iloc[n].str.cat())
        f.write(wcr_csv.iloc[n].astype(str).str.ljust(12).str.cat(sep=' '))
        f.write('\n')


# %% [markdown]
# ## Variograms

# %%
wcr_elev = wcr_csv.copy()
wcr_elev['elev'] = wcr_elev.Zland - wcr_elev.Depth

# %%
print(wcr_elev.elev.median())
# subset a layer for identifying variogram spatially
h = -10
h_diff = 2
def subset_wcr(wcr_elev, h, h_diff):
    wcr_elev_var = wcr_elev[(wcr_elev.elev>h-h_diff)&(wcr_elev.elev<h)].copy()
    # x, y, field =  wcr_elev_var.Easting, wcr_elev_var.Northing, wcr_elev_var.classification
    x, y, z, field =  wcr_elev_var.X, wcr_elev_var.Y, wcr_elev_var.elev, wcr_elev_var.PC
    return x, y, z, field
x, y, z, field = subset_wcr(wcr_elev, -10, 2)
fig, ax=  plt.subplots()
ax.scatter(x,y, c=field)
domain_wcr.plot(ax=ax, color='None')

# %%
import gstools as gs
from pykrige import OrdinaryKriging

# %%
pos = (x,y)
pos = (x,y,z)
bin_edges = gs.variogram.standard_bins(pos=pos, dim=len(pos), max_dist=4000)

# %%
bin_edges

# %%
# specifying direction angle creates two variogram outputs for each direction
angle = 52.9 *(np.pi / 180)
bin_center, vario = gs.vario_estimate(
    pos, field, 
    bin_edges = bin_edges, 
    direction=gs.rotated_main_axes(dim=len(pos), angles=angle),
    angles_tol = 5*(np.pi/180),
    # estimator='cressie'
)


# %%
# length scale can be in each direction or list with main and transversal directions
cov_model_fit = gs.Exponential(
    dim=len(pos), len_scale=[200,200,10],  var=0.01, angles=angle
)

cov_model_fit.fit_variogram(bin_center, vario, nugget=True);
# # model_fit.plot()
# ax = cov_model_fit.plot(x_max=max(bin_center))
# ax.scatter(bin_center, vario)

# %%
print('Sill %.4f,' %cov_model_fit.sill, 'nugget %.4f,' %cov_model_fit.nugget, 'length scale %.4f,' %cov_model_fit.len_scale)

# %% [markdown]
# If the squared difference decreases with higher distance then we should set a max distance. Alisha has a max of 4000 m in x, 2500 m in y and 20m in Z.
# - this fixed the issue
#
# Issue with some early points having large differences.
# - Going to a 3D variogram fitting provided the best results (exponential better than gaussian)
# - can't apply the 3D variogram to the z direction because the bins are set for the X and Y.
#
# Need to re-run with Z as a separate version

# %%
fig, ax =plt.subplots(3,1)
# ax = cov_model_fit.plot(x_max=max(bin_center), ax=ax)
# ax.scatter(bin_center, vario)
for n in np.arange(0,3):
    ax[n].scatter(bin_center, vario[n], label="emp. vario: 52.9")
# ax.scatter(bin_center, vario[1], label="emp. vario: 52.9 * 5?")
    cov_model_fit.plot("vario_axis", axis=n, ax=ax[n], x_max=max(bin_center), label="fit on axis ")
# cov_model_fit.plot("vario_axis", axis=1, ax=ax, x_max=max(bin_center), label="fit on axis 1")

# %% [markdown]
# ## T2py library

# %%
add_path(join(doc_dir, 'GitHub','t2py'))


# %%
import t2py

# %%
# t2py.WellLogFile(dataframe=wcr_csv)

# %% [markdown]
# # Output
# Review texture to PAR output

# %%
m = flopy.modflow.Modflow.load(join(model_ws, 'MF.nam'))

# %%
lay_bot = m.dis.botm.array[:,0,0]
n=1
lay_bot[n-1], lay_bot[n-1]-lay_bot[n]

# %%
upw = flopy.modflow.ModflowUpw.load(join(model_ws, 'MF.upw'), m)
hk = upw.hk.array
vmin, vmax = hk.min(), hk.max()
# vmin = 10**(np.floor(int(np.log10(vmin))))

# %%

# %%
nlay = m.dis.nlay-1
ny=3
nx=int(np.ceil(nlay/ny))
fig, ax = plt.subplots(nx, ny, sharex=True, sharey=True, figsize=(6.5, 8), dpi=300)
for n in np.arange(1,nlay):
    ax_n = ax[int((n-1)/ny), (n-1)%ny]
    # need to get grid xy for wcr first
    # x, y, z, field = subset_wcr(wcr_elev, lay_bot[n-1], lay_bot[n-1]-lay_bot[n])
    # ax_n.scatter(x,y, c=field)
    im = ax_n.imshow(hk[n], norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
# m.dis.nlay
fig.tight_layout()
fig.colorbar(im, ax=ax, location = 'top', orientation='horizontal', shrink=0.7, label='Kx')


# %% [raw]
# Using default parameters for kriging didn't work, it created a pilot point heavy result with odd wavy patterns. Likely the range was much too small or the fine texture facies are overweighting the system. 
# - In some cases fewer wells for kriging may help identify local paths while more might yield more regional results
# - increasing the range created larger pilot points but didn't improve any patterns ( had been using a spherical variogram)
#
#
# The model was running quickly with kriging with 16 nearest wells, and I bumped it to 160 which took several minutes
# - increasing the number of wells didn't seem to remove the pilot point nature of the simulation.
