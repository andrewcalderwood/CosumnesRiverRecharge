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
# Goal: Run scripts needed to extend the model into the future.
#
# Data downloads:
# - CIMIS potential evapotranspiration and rainfall (soil water budget)
# - USGS Streamflow at Michigan Bar (stream inflow)
# - DWR groundwater elevation data base (lateral boundaries and observation matching)
# - DWR stream stage at McConnell (observation matching)
#
# Scripts to run:
# - SFR/ add script to download USGS data, I have an example of this somewhere
# - UZF/Basic soil budget - fields, check if gridded is still needed for natives?
#     - need to download another year of CIMIS precip/rainfall
# - GHB/ need to run scripts 1 to 3 but these need updates to move testing code to separate notebooks

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


# %% [markdown]
# # Load data files

# %%

# %% [markdown]
# # GHB Update
#
# All the scripts in the GHB folder should be run in numeric order.  
#
# - `00_Preliminary DWR GW data analysis.py` loads the DWR periodic measurements and crops to the buffered extent of the model domain. Then saves the points as shapefiles by season and year. The DWR periodic measurements file has to be manually downloaded currently, should update to scripted.  
# - `01_Kriging_for_GHB_boundary.py` takes the point shapefiles and applies kriging to each year/season combination to create rasters of groundwater elevation. It takes specific years overwhich to run for the contours since we don't want to have to always re-create past years.
#
# - `02_Distance GHB.py` and `03_GHB head correction.py` don't depend on the years as they run on the entire datasets available  
