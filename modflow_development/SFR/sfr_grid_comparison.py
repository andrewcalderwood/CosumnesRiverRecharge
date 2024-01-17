# ---
# jupyter:
#   jupytext:
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

# %% [markdown]
# Quick script to compare the local and regional sfr reach definition, specifically the stream top elevations. This is to confirm whether it is impacting the potential for baseflow by the lower floodplains.

# %%
# standard python utilities
import os
from os.path import join, basename, dirname, exists,expanduser
import sys
import glob
import pandas as pd
import numpy as np
import time

import geopandas as gpd

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')

# dir of all gwfm data
gwfm_dir = join(usr_dir,'Box/research_cosumnes/GWFlowModel')
print(gwfm_dir)
sfr_dir = join(gwfm_dir,'SFR_data')


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy

py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

from mf_utility import get_dates, get_layer_from_elev, clean_wb
from flopy_utilities import reach_data_gdf
from map_cln import gdf_bnds, plt_cln

# %%
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes'


# %% [markdown]
# ## Load data

# %%
# regional grid_sfr
grid_sfr = gpd.read_file(join(sfr_dir, 'final_grid_sfr/grid_sfr.shp'))

# %%

model_nam = 'Stream_seepage/oneto_denier_upscale4x_2014_2020'
model_ws = join(loadpth, model_nam)
grid_sfr_l = pd.read_csv(join(model_ws, 'grid_sfr.csv'),index_col=0)

# %%
m = flopy.modflow.Modflow.load(model_ws+'/MF.nam',  load_only=['DIS','SFR'])

# %%
model_grp = 'inset_oneto_denier'
grid_dir = join(gwfm_dir, 'DIS_data/streambed_seepage/grid')
grid_fn = join(grid_dir, model_grp,'rm_only_grid.shp')
grid_p = gpd.read_file(grid_fn)
grid_p.crs='epsg:32610'
m_domain = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [grid_p.unary_union], crs=grid_p.crs)

# %%
grid_sfr_l = reach_data_gdf(m.sfr, grid_p)
grid_sfr_l = grid_sfr_l[grid_sfr_l.strhc1!=0]
grid_sfr_l['local_reach'] = np.arange(0, len(grid_sfr_l))

# %% [markdown]
# ## merge regional and local data to compare stream top elevation

# %%
cols = ['geometry','iseg','ireach','rchlen','strtop','slope']
cols = ['geometry','strtop','local_reach']
grid_sfr_join = gpd.sjoin_nearest(grid_sfr, grid_sfr_l[cols], how='right')

# %%
grid_sfr_join.sort_values('local_reach').plot(x='local_reach',y=['strtop','z_min'],kind='line')

# %%
diff = grid_sfr_join.z_min - grid_sfr_join.strtop

print('Regional grid sfr is on avg %.2f m higher than local model' %diff.mean())
print('Std dev of %.2fm ' %diff.std())

# %% [markdown]
# This generally explains why the regoinal model can't show baseflow because the streamtop elevations are on average 2 meters higher and there also is generally more drawdown via pumping.
