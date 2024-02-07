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


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

from importlib import reload

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

from report_cln import base_round

# %% [markdown]
# # Load data and crop to model extent

# %%
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')

# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_elev = gpd.read_file(join(gwfm_dir,'DIS_data','grid_elevation_m_statistics.shp'))

# %%
grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')

# %%
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
rivers = gpd.read_file(join(sfr_dir,"Sac_valley_rivers/Sac_valley_rivers.shp"))
rivers = rivers.to_crs('EPSG:32610')

rivers_clip = gpd.clip(rivers, m_domain)
rivers_clip = rivers_clip[['GNIS_ID','GNIS_Name','LengthKM', 'geometry']]
cr = rivers_clip[rivers_clip.GNIS_Name=='Cosumnes River']
mr = rivers_clip[rivers_clip.GNIS_Name=='Mokelumne River']

# %%
aem_folder = 'statewide_aem_survey_coarse_fraction_depth_slices_and_average'
aem_domain_f = join(upw_dir, aem_folder, 'domain_aem_data.shp')
if not exists(aem_domain_f):
    aem_fp = glob.glob(join(upw_dir, aem_folder)+'/*.shp')[0]
    aem = gpd.read_file(aem_fp)
    # buffer domain by a few km to include nearby AEM data
    m_buf = m_domain.copy()
    m_buf.geometry = m_buf.geometry.buffer(5E3)
    # crop AEM data to the model domain area (runs pretty quickly
    aem_domain = gpd.overlay(aem, m_buf.to_crs(aem.crs)).to_crs(m_domain.crs)
    
    aem_domain.to_file(aem_domain_f)
else:
    aem_domain = gpd.read_file(aem_domain_f)

# %% [markdown]
# # Clean coarse fraction to plot

# %%
# calculate the average coarse fraction
aem_domain[['PC_low','PC_high']] = aem_domain.PercentCoa.str.split('-', expand=True).astype(float)
aem_domain['PC_avg'] = aem_domain[['PC_low','PC_high']].mean(axis=1)

# pull out the upper 15 m average for plotting
aem_avg = aem_domain[aem_domain.Interval_T=='Average, 0-15m'].copy()


# %%
fig,ax = plt.subplots(figsize=(6.5,6.5),dpi=300)
# Spectral_r or viridis work well for plotting
aem_avg.plot('PC_avg', ax=ax, cmap='viridis', legend=True)
m_domain.plot(color='none',ax=ax)
cr.plot(color='black', linewidth=1, ax=ax)

ctx.add_basemap(source=ctx.providers.Esri.WorldImagery, ax=ax, alpha = 0.6,
                crs=m_domain.crs, attribution=False)

# %% [markdown]
# ## Depth slices  
# The DWR dataset joins adjacent polygons that have the same percentcoa which creates larger polygons that can span multiple cross-sections.
# - To plot the cross-sections I should use overlay and I likely need to conver to a grid format to have a way to specify stream reach for the x-axis.
# Data Summary  
# - The online data shows the AEM goes to roughly 300 m below ground in the profile viewer. The Depth_Mid_ data is likely in meters because the maximum value is 330 m. The layer thickness ranges from 5 - 44 m.
# - the coarse and fine fraction (splitting at 50%) align with WCR estimates of aroud 75% fine and 25% coarse. Gravel-mud fractions are 0.057, 0.207, 0.151, 0.584

# %%
aem_depth = aem_domain[aem_domain.Interval_T=='Depth Slice'].copy()

# %%
# there are large polygons
# aem_depth[aem_depth.area>1E8].iloc[[0]]

# %%
fig,ax= plt.subplots()
gpd.sjoin(aem_depth, cr, how='inner').plot(ax=ax)
aem_depth.iloc[[831]].plot(ax=ax,color='red')
aem_depth[aem_depth.area>1E8].iloc[[0]].plot(ax=ax, color='red')
cr.plot(ax=ax)

# %%
thick = -(aem_depth.Depth_Top_ - aem_depth.Depth_Bott)
print('Min thick %.1f' %thick.min(), 'max thick %.1f' %thick.max())
print('Max depth %.1f' %aem_depth.Depth_Bott.max())

# %%
# define the cutoffs between facies
# none of the percent coarse actually equals zero, but some equals 5%
facies_frac = [0, 40, 50, 80,100] # matches new Cosumnes proportions
# facies_frac = [0, 30, 40, 60,100]
facies = ['mud','sandy mud', 'sand', 'gravel']
# default is mud as most common, and > lower means 0 must be included as default
aem_depth['facies'] = 'mud' 
for n, f in enumerate(facies):
    aem_depth.loc[(aem_depth.PC_avg>facies_frac[n]) & (aem_depth.PC_avg<=facies_frac[n+1]),'facies'] = f

# assign facies numbers
facies_num = pd.DataFrame(np.arange(4,0,-1),facies, columns=['facies_num'])
aem_depth['facies_num'] = facies_num.loc[aem_depth.facies].values

aem_depth.to_file(join(upw_dir, aem_folder, 'aem_depth.shp'))


# %%

# %%
# shallow data to be more representative of alluvium - need to decide alluvium bottom
aem_shallow = aem_depth[aem_depth.Depth_Bott<100]
total_area = aem_shallow.geometry.area.sum()
aem_fine = aem_shallow[aem_shallow.PC_avg<50]
fine_fraction = 100*aem_fine.geometry.area.sum()/total_area
aem_coarse = aem_shallow[aem_shallow.PC_avg>=50]
coarse_fraction = 100*aem_coarse.geometry.area.sum()/total_area
print('Coarse fraction %.2f %%' %coarse_fraction, 'fine fraction %.2f %%' %fine_fraction)

for f in facies:
    area_frac = 100*aem_shallow[aem_shallow.facies==f].geometry.area.sum()/total_area
    print(f,'fraction %.2f %%' %area_frac)

# %% [markdown]
# Using 25, 50, 75, 100 as the initial cutoffs yielded a slightly smaller than expected mud fraction. Increasing from 25 to 40 gave the closest match as 30 underestimated the fraction.
# - Keep in mind that the percent coarse are on the scale of 0-10 so the changes in fraction must be by 10%
#
# - After reviewing Casey Meirovtiz's thesis I realized that he reported a much higher portion of coarse vs fine with more for the american than the cosumnes.
# - Decreasing the percent coarse threshold for sand-mud, sand and gravel helped better delineate the two shallow coarse deposits under the Cosumnes River when previously the slightly deeper one () was defined as fine instead of coarse |

# %% [markdown]
# | Facies | American | Cosumnes |
# | ----------- | ----------- | ----------- |
# | Channel | 0.21 | 0.13 |
# | Near-channel | 0.13 | 0.16 |
# |Proximal-channel/levee| 0.26 | 0.22| 
# |Floodplain| 0.4 | 0.49 |    

# %% [markdown]
# # Grid transform
# Overlay the model grid onto the AEM data to identify the facies expected in each grid cell. The issue is that I don't want to apply the full 1 km buffer from the output of the results so I may want to do anti bufer of 500 m first.

# %%
# 3820 full, 3328 with -200 buffer, 
# anti-buffer to remove excess interpolation
buf = -200
aem_buf = aem_depth[aem_depth.buffer(buf).geometry.area!=0].copy()
aem_buf##.plot()
# overlay the grid onto the adjusted AEM polygons
aem_grid = gpd.overlay(aem_buf, grid_p)
aem_grid['elev_m'] = dem_data[aem_grid.row-1, aem_grid.column-1] # add elevation for depths

# %%
d_cols = ['Depth_Mid_','Depth_Bott', 'Depth_Top_']
# create new columns of elevation
for n, d in enumerate(d_cols):
    aem_grid[d.replace('Depth','Elev')] = aem_grid['elev_m'] - aem_grid[d] 


# %%
# slow to write
# could also save the grid join for just the top layer?
aem_grid.to_file(join(upw_dir, aem_folder, 'aem_grid.shp'))


# %%
aem_grid_sfr = gpd.overlay(aem_buf, grid_sfr)

# %%
# aem_grid_sfr = aem_grid.merge( grid_sfr, on=['row','column','node'])
aem_grid_sfr = aem_grid[aem_grid.node.isin(grid_sfr.node)].reset_index()
sfr_cols= ['z','slope','z_min','reach','length_m']
sfr_cols = grid_sfr.set_index(['node']).loc[aem_grid_sfr.node.astype(int), sfr_cols].reset_index()
aem_grid_sfr = pd.concat((aem_grid_sfr, sfr_cols), axis=1)

# %%
import matplotlib as mpl
norm = mpl.colors.Normalize(vmin=0, vmax=100)
# n=50
# mpl.cm.viridis(norm(n))
# colors = mpl.cm.viridis(norm(aem_depth.PC_avg.values))

# %%
# grid_sfr.plot('reach', legend=True)

# %%
# aem_grid_sfr.plot(x='reach',y='Depth_Mid_')
colors = mpl.cm.viridis(norm(aem_grid_sfr.PC_avg.values))
plt.scatter(aem_grid_sfr.reach.values, aem_grid_sfr.Elev_Mid_.values, color=colors)

# %% [markdown]
# ## convert data to 3D arrays
# TPROGS arrays are 80 to -80 m amsl so the array shape should match that with 0.5 m steps.
#
# The ranges in the AEM data are larger so it might make more sense to use them to identify the upper and lower extents and assign values rather than start from the tprogs layers

# %%
nlay = int((80 - -80)/0.5)
nrow = 100
ncol = 230
# nlay=640 # just for plot testing to check inclusion of deeper facies

aem_array = np.zeros((nlay, nrow, ncol))
aem_cf_array = np.zeros((nlay, nrow, ncol))


# %%
# aem_grid.Elev_Mid_.min()
# for b in np.arange(80, -80, -0.5):
#     b

# %%
def elev_2_tprogs_lay(elev):
    """
    Elevation (m) to tprogs layers
    """
    elev = base_round(elev, 0.5)
    lay = int((80 - elev)/0.5)
    return lay



# %%
# a little slow but not bad (10 seconds or so)
for n in aem_grid.node.unique():
    aem_n = aem_grid[aem_grid.node==n]
    row_n = aem_n.row.iloc[0]-1
    col_n = aem_n.column.iloc[0]-1
    for k in np.arange(0, len(aem_n)):
        lay_b = elev_2_tprogs_lay(aem_n.Elev_Bott.iloc[k])
        lay_t = elev_2_tprogs_lay(aem_n.Elev_Top_.iloc[k])
        # assign the facies from the aem data
        aem_array[lay_t:lay_b, row_n, col_n] = aem_n.facies_num.iloc[k]
        aem_cf_array[lay_t:lay_b, row_n, col_n] = aem_n.PC_avg.iloc[k]
# base_round(aem_n.Elev_Top_.iloc[n], 0.5)


# %%
aem_mask = np.ma.masked_where(aem_array==0, aem_array)
aem_cf_mask = np.ma.masked_where(aem_cf_array==0, aem_cf_array)

# %%
fig,ax = plt.subplots(2,1, sharex=True)
k=160
# k=100
ax[0].imshow(aem_mask[k], cmap='viridis_r')
ax[1].imshow(aem_cf_mask[k], cmap='viridis')
fig.tight_layout(h_pad=0.1)

# %%
sfr_profile = aem_mask[:,grid_sfr.row.astype(int)-1, grid_sfr.column.astype(int)-1]
sfr_profile_cf = aem_cf_mask[:,grid_sfr.row.astype(int)-1, grid_sfr.column.astype(int)-1]
fig,ax = plt.subplots(2,1,sharex=True, figsize=(6.5,6))
ax[0].imshow(sfr_profile, cmap='viridis_r')
ax[1].imshow(sfr_profile_cf, cmap='viridis')
ax[0].set_aspect(1/5)
ax[1].set_aspect(1/5)
fig.tight_layout(h_pad=-4)


# %%
import h5py

# %%
with h5py.File(join(upw_dir, aem_folder, 'aem_array.hdf5'), 'w') as f:
    grp = f.require_group("facies")
    dset = grp.require_dataset("facies_array", aem_array.shape, dtype='i')
    dset[:] = aem_array
    grp = f.require_group("percent_coarse")
    dset = grp.require_dataset("pc_array", aem_cf_array.shape, dtype='i')
    dset[:] = aem_cf_array

# %%
with h5py.File(join(upw_dir, aem_folder, 'aem_array.hdf5'), 'r') as f:
    arr = f['facies']['facies_array'][:]
    arr_cf = f['percent_coarse']['pc_array'][:]
# arr_ma = np.ma.masked_where(arr==0,arr)
# plt.imshow(arr_ma[100])
