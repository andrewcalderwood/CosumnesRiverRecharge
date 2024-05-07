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
# Cosumnes Model 
#

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import sys
import glob
from importlib import reload

import pandas as pd
import numpy as np
from scipy.stats import hmean, gmean

# import calendar
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

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



# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
print(gwfm_dir)
sfr_dir = gwfm_dir+'/SFR_data/'


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')

import flopy 

from importlib import reload
# importlib.reload
# reload(flopy)

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

from mf_utility import get_layer_from_elev, param_load

# %% [Markdown] [markdown]
# # Correct for elevations about the boundary land surface
# Heads in the foothills are more uncertain in the contouring both in the domain and beyond.

# %%
ghb_dir = join(gwfm_dir, "GHB_data")

dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
## heads from regional kriging at a distance from the domain (5 km with some smaller due to raster constraints)
bnd_dist = pd.read_csv(join(ghb_dir, 'boundary_distance_heads.csv'))
# add month column for easier date creation
bnd_dist['month'] = 'Oct'
bnd_dist.loc[bnd_dist.season=='spring','month'] = 'Apr'
# convert from ft to meters and match variable with existing dataset
bnd_dist['wse_m'] = bnd_dist.wse_ft*0.3048
bnd_dist['value'] = bnd_dist.wse_ft*0.3048
bnd_dist[['row','column']] -=1
# remove column that won't work with resampling
bnd_dist = bnd_dist.drop(columns=['season'])

# %%
print('Span:',bnd_dist.year.min(), '-',bnd_dist.year.max())

# %%
# may need to reconsider boundary GWE that are driving inflow when it should be more static in the foothills
bnd_dist['dem_m'] = dem_data[bnd_dist.row, bnd_dist.column]
# alternatively we could switch to using a long-term average GHB head

# %%
# bnd_dist_avg = bnd_dist.groupby(['row','column']).mean(numeric_only=True).reset_index()
# the average of heads still show values above land surface

# %%
# bnd_dist[bnd_dist.value > bnd_dist.dem_m - 5]
# use the deep geology to identify region to correct boundary groundwater flow
# adj_min_col = np.where(deep_geology[0])[1].min()
# adjust 140 and above but only keep +10
adj_min_col = 140
bnd_dist_adj = bnd_dist[bnd_dist.column>adj_min_col].copy()
# groundwater elevations shouldn't be above land surface
adj_bool = bnd_dist_adj.value > bnd_dist_adj.dem_m - 5
bnd_dist_adj.loc[adj_bool, 'value'] = bnd_dist_adj.loc[adj_bool,'dem_m'] -5


# %% [markdown]
# # iterate over each year, season and row

# %%
# find unique year, season, row combos
date_combo = bnd_dist_adj[['year','month','row']].drop_duplicates().reset_index(drop=True)
for n in np.arange(0, len(date_combo)):
    y = date_combo.loc[n, 'year']
    month = date_combo.loc[n, 'month']
    row = date_combo.loc[n, 'row']
    temp_plt = (bnd_dist_adj.year==y)&(bnd_dist_adj.month==month)&(bnd_dist_adj.row==row)
    rper = 15
    temp = bnd_dist_adj[temp_plt].rolling(rper, min_periods=0,on='column', center=True)['value'].mean(numeric_only=True)
    bnd_dist_adj.loc[temp_plt, 'value'] = temp

# filter out to only the upper columns to avoid increasing in elevation 
bnd_dist_adj = bnd_dist_adj[bnd_dist_adj.column>= adj_min_col +rper]

# %%
temp_plt = (bnd_dist_adj.year==2019)&(bnd_dist_adj.month=='Apr')&(bnd_dist_adj.row==99)
fig,ax = plt.subplots()
bnd_dist_adj[temp_plt].plot(x='column', y=['value','wse_m','dem_m'], ax=ax)


# %%
bnd_dist['keep'] = True
bnd_dist = bnd_dist.set_index('node')
bnd_dist.loc[bnd_dist_adj.node.values, 'keep'] = False
bnd_dist = bnd_dist.reset_index()
bnd_out = pd.concat((bnd_dist[bnd_dist['keep']], bnd_dist_adj))
bnd_out.to_csv(join(ghb_dir, 'boundary_distance_heads_adjusted.csv'), index=False)
