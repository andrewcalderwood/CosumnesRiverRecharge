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

# %% [markdown]
# # Identify the spatial extent of cooarse pathways
# The Cosumnes model is based only on geologic data in the Cosumnes River area and doesn't include the incised valley fill (IVF) from the American River fan. The AEM data shows there is a strong appearnce of the IVF below the middle lower Cosumnes River. This means it should be included in the geologic model, but rather than produce a TProGS model of the American and crop out the channel deposits we can use the AEM data combined with the cross-sections and maps from Meirovitz (2010) to identify the true channel extents.

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

from mf_utility import get_layer_from_elev
from report_cln import base_round

# %% [markdown]
# # Load data

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

# %%
aem_depth = gpd.read_file(join(upw_dir, aem_folder, 'aem_depth.shp'))


# %%
# slow to read because of high cell count
aem_grid = gpd.read_file(join(upw_dir, aem_folder, 'aem_grid.shp'))


# %%
# default of 10 bins+1 has 40 m thickness
interval = pd.IntervalIndex(pd.cut(aem_grid.Elev_Mid_, bins=20))
aem_grid['bin_lower'] = interval.left.values
aem_grid['bin_higher'] = interval.right.values

# %% [markdown]
# # Visualize

# %%
vmin = aem_depth.PC_avg.min()
vmax = aem_depth.PC_avg.max()
# depths =np.sort(aem_depth.Depth_Mid_.unique())[:-2:2]
depths =np.sort(aem_grid.bin_higher.unique())[6:-2]
nx = 3
ny = int(np.ceil(len(depths)/nx))
fig,ax = plt.subplots(ny, nx, sharex=True, sharey=True, figsize=(6.5, 6.5*(ny/nx)), dpi=300)
for n, d in enumerate(depths):
    ax_n = ax[int(n/nx), n%nx]
    aem_plt = aem_grid[aem_grid.bin_higher==d]
    # aem_plt = aem_depth[aem_depth.Depth_Mid_==d]
    aem_plt.plot('PC_avg', ax=ax_n, vmin = vmin, vmax=vmax)
fig.tight_layout()

# %%
fig,ax= plt.subplots()
aem_depth.plot('PC_avg', ax=ax, legend=True)
cr.plot(ax=ax)

# %%
