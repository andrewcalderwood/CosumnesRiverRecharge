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
from os.path import basename, dirname, join, exists
import glob

import pandas as pd
import numpy as np
import numpy.ma as ma
from scipy.stats import gmean

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# standard geospatial python utilities
import geopandas as gpd
# mapping utilities
# import contextily as ctx

# %%
doc_dir = os.getcwd()
while os.path.basename(doc_dir) != 'Documents':
    doc_dir = os.path.dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = os.path.dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study
proj_dir = gwfm_dir + '/Oneto_Denier/'
dat_dir = proj_dir+'Stream_level_data/'

sfr_dir = gwfm_dir+'/SFR_data/'

# %%
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')
# load parent model grid
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
nrow = int(grid_p.row.max())
ncol = int(grid_p.column.max())


# %%
# get exterior polyline of model grid
grid_bnd = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [grid_p.unary_union.exterior], crs=grid_p.crs)
# find cells that construct the model boundary
bnd_gdf = gpd.sjoin(grid_p, grid_bnd)
bnd_cells = bnd_gdf[['row','column']] - 1
bnd_cells['grid_id'] = np.arange(0,len(bnd_cells))
rows,cols = bnd_cells.row.values, bnd_cells.column.values

# %%
bnd_gdf

# %%
strtyear = 2011 # 2010 isn't a good year yet
endyear = 2021

# %%
# raster cropping will be done in outside script so the only part read in will be the final array
ghb_dir = gwfm_dir+'/GHB_data'

# %%

kriged_fall = np.zeros((int(endyear-strtyear),nrow,ncol))
kriged_spring = np.zeros((int(endyear-strtyear),nrow,ncol))

# keep track of which place in array matches to year
year_to_int = np.zeros((endyear-strtyear,2))

for t, year in enumerate(np.arange(strtyear,endyear)):
    # load and place spring kriged data in np array, load spring first
    filename = glob.glob(ghb_dir+'/final_WSEL_arrays/spring'+str(year)+'_kriged_WSEL.tsv')[0]
    kriged_spring[t,:,:] = np.loadtxt(filename)
    # load and place fall kriged data in np array
    filename = glob.glob(ghb_dir+'/final_WSEL_arrays/fall'+str(year)+'_kriged_WSEL.tsv')[0]
    kriged_fall[t,:,:] = np.loadtxt(filename)
    
    year_to_int[t,0] = t
    year_to_int[t,1] = year


# %%
# ceate multi, index to stack fall and spring data
sy_ind = np.repeat(['Apr','Oct'],(endyear-strtyear)),np.tile(np.arange(strtyear,endyear),2)
sy_ind = pd.MultiIndex.from_arrays(sy_ind, names=['month','year'])


# %%
# Set kriged water table elevations that are above land surface to land surface minus 15 ft (based on historical levels)
# no observed water elevations show depth closer than 15 ft even in foothills
dem_offset = 15*0.3048

# stack fall and spring before resampling
kriged = np.vstack((kriged_spring[:, rows, cols],kriged_fall[:, rows, cols]))
kriged = np.where(kriged>dem_data[rows, cols], dem_data[rows, cols]- dem_offset, kriged)


# %%

# # should bring back monthly interpolate along defined boundaries
kriged_df = pd.DataFrame(np.transpose(kriged),columns=sy_ind)
# long format for easier resampling and create datetime column
df_long = kriged_df.melt(ignore_index=False).reset_index(names='grid_id') # keep index it row or col number
df_long['date'] = pd.to_datetime(df_long.year.astype(str)+'-'+df_long.month)
# linearly interpolate between fall and spring measurements for each row,col id
df_mon = df_long.set_index('date').groupby(['grid_id']).resample('MS').interpolate('linear')
df_mon = df_mon.reset_index('grid_id', drop=True)
df_mon['year'] = df_mon.index.year
df_mon['month'] = df_mon.index.month

df_mon = df_mon.join(bnd_cells.set_index('grid_id'),on='grid_id')

# %%
sns.relplot(df_mon[df_mon.row==99], x='column',y='value', col='year', hue='month', col_wrap=4, kind='line')

# %%
