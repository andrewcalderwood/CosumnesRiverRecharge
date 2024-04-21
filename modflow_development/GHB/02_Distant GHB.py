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

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename
import sys
import glob

import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio

import matplotlib.pyplot as plt

# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
gwfm_dir

# %%

# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')


# %%
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')


# %%
ghb_dir = join(gwfm_dir, 'GHB_data')
# kriged water surface elevation to use
tif_path = join(ghb_dir, 'interpolated_data')
# kriged_tifs = glob.glob(join(tif_path,'*_kriged.tif'))


# %%
# get exterior polyline of model grid
grid_bnd = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [grid_p.unary_union.exterior], crs=grid_p.crs)
# find cells that construct the model boundary
bnd_cells_df = gpd.sjoin(grid_p, grid_bnd)
bnd_cells_df = bnd_cells_df.drop(columns=['index_right',0])
bnd_cells = bnd_cells_df[['row','column']] - 1
bnd_cells['grid_id'] = np.arange(0,len(bnd_cells))
bnd_rows, bnd_cols = bnd_cells.row.values, bnd_cells.column.values

# %% [markdown]
# Before 2012 there is not as many WSE points so the raster is even smaller so part of the domain is even cutoff.

# %%
raster_name = join(tif_path, 'fall2012_kriged.tif')
dem_f = rasterio.open(raster_name)
# dem = dem_f.read(1)
xmin, xmax = dem_f.bounds[0], dem_f.bounds[2]
ymin, ymax = dem_f.bounds[1], dem_f.bounds[3]
print(dem_f.meta)
from shapely.geometry import box
rio_box = gpd.GeoDataFrame([0], columns=['id'], geometry = [box(*dem_f.bounds)], crs=dem_f.crs)
rio_box.geometry = rio_box.geometry.exterior

dem_f.close()

# %%
# fig, ax = plt.subplots()
# rio_box.plot(ax=ax)
# bnd_cells_df.plot(ax=ax)

# %%
dist = 5000

# for cells buffered that would be outside of the raster use the furthest distance possible
bnd_max_dist = gpd.sjoin_nearest(bnd_cells_df, rio_box, distance_col='dist_m')
# subtract 100 m from distance to ensure a kriged cell is sampled
bnd_max_dist['ghb_dist'] = np.round(bnd_max_dist.dist_m-100,-2)
bnd_max_dist.loc[bnd_max_dist.dist_m>dist, 'ghb_dist'] = dist
bnd_max_dist = bnd_max_dist.drop(columns=['index_right','id'])

# %%
nrow = bnd_max_dist.row.max()

# %%
bnd_nw =  bnd_max_dist[bnd_max_dist.row==1].copy()
bnd_se =  bnd_max_dist[bnd_max_dist.row==nrow].copy()

# %%
rot = 52.9*np.pi/180 # rotation of the model grid
# change in x and y off of northwest and southeast
x= np.sin(rot)*dist
y= np.cos(rot)*dist


# %%
# northwest goes negative x and postive y
bnd_nw_dist = bnd_nw.copy()
bnd_nw_dist.geometry = bnd_nw.translate(-x, y)
# southeast goes positive x and negative y
bnd_se_dist = bnd_se.copy()
bnd_se_dist.geometry = bnd_se.translate(x, -y)

# %%
# only the southeast boundary has an overlap
bnd_fix = bnd_se[bnd_se.ghb_dist<dist].copy()
# subset for where ghb_dist < dist
for n in np.arange(0, len(bnd_fix)):
    x = np.sin(rot)*bnd_fix.iloc[n].ghb_dist
    y = np.cos(rot)*bnd_fix.iloc[n].ghb_dist
    bnd_fix.geometry.iloc[n] = bnd_fix.geometry.translate(x,-y).iloc[n]

# %%
bnd_se_dist_full = pd.concat((bnd_se_dist[~bnd_se_dist.node.isin(bnd_fix.node)], bnd_fix))

# %%
fig, ax = plt.subplots()
bnd_nw_dist.plot(ax=ax,)
bnd_nw.plot(ax=ax,color='red')

bnd_se_dist.plot(ax=ax, color='green')
bnd_se.plot(ax=ax,color='red')
bnd_se_dist_full.plot(ax=ax)
# plt.legend(labels=['NW Offset', 'Old','SE Offset','Old'])

# %%
# bnd_dist = pd.concat((bnd_nw_dist, bnd_se_dist)) # with straight lines that miss data
bnd_dist = pd.concat((bnd_nw_dist, bnd_se_dist_full))
bnd_dist['easting'] = bnd_dist.geometry.centroid.x
bnd_dist['northing'] = bnd_dist.geometry.centroid.y
bnd_dist

# %%
# kriged_tifs
season='fall'
y=2012
raster_name = join(tif_path, season+str(y)+'_kriged.tif')


# %%
# np.repeat(bnd_dist, 2, axis=0)
bnd_dist['season'] = season
bnd_dist['year'] = y
bnd_dist

# %%


# dem_f = rasterio.open(raster_name)
# dem = dem_f.read(1)
# xmin, xmax = dem_f.bounds[0], dem_f.bounds[2]
# ymin, ymax = dem_f.bounds[1], dem_f.bounds[3]
bnd_dist_all = pd.DataFrame()

point = bnd_dist.loc[:,['easting','northing']].values
for season in ['spring', 'fall']:
    for y in np.arange(2012,2021):
#     for y in [2012]:
        raster_name = join(tif_path, season+str(y)+'_kriged.tif')

        with rasterio.open(raster_name) as src:
            bnd_dist['wse_ft'] = [sample[0] for sample in src.sample(point)]
            bnd_dist.loc[bnd_dist.wse_ft==src.nodata,'wse_ft'] = np.nan
            # save information for all years and seasons
            bnd_dist['season'] = season
            bnd_dist['year'] = y
            bnd_dist_all = pd.concat((bnd_dist_all, bnd_dist))


# %%

# %%
import seaborn as sns

# %%
sns.relplot(bnd_dist_all[bnd_dist_all.row==1], 
            x='column', y='wse_ft', hue='year', col='season')

# %%
sns.relplot(bnd_dist_all[bnd_dist_all.row==nrow], 
            x='column', y='wse_ft', hue='year', col='season')

# %% [markdown]
# A piece of the northeast boundary in the mountains is missing in spring 2012

# %%
pd.DataFrame(bnd_dist_all.drop(columns=['geometry'])).to_csv(join(ghb_dir, 'boundary_distance_heads.csv'), index=False)

# %%
bnd_dist.plot('wse_m', legend=True)


# %% [markdown]
# # Correct for elevations about the boundary land surface
# Heads in the foothills are more uncertain in the contouring both in the domain and beyond.

# %%

# %%
