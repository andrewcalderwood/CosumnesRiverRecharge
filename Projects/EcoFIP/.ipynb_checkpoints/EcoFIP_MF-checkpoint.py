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
# Load the existing historical modflow model and adapt only the packages necessary for the EcoFIP anlaysis. 
# - If restoration includes existing farmland then WEL/RCH need to have those cells removed
# - Convert SFR to RIV for initial testing at project locations
# - Research further the use of MODFLOW-API to connect to HEC-RAS. This won't be fully done, because EcoFIP doesn't run HEC-RAS in a continuous series but assumes a water surface given a flow time series so it would be that we could try to have an adjustment factor in the model by estimating seepage so that downstream segments have lower flow so lower water surface elevations.

# %%
# standard python utilities
import os
from os.path import basename, dirname, join, exists, expanduser
import sys
from importlib import reload
import glob

import pandas as pd
import numpy as np
import time

# standard python plotting utilities
# import seaborn as sns
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates

# standard geospatial python utilities
import geopandas as gpd
import rasterio

# mapping utilities
# import contextily as ctx
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as fm
# from matplotlib.ticker import MaxNLocator


# %%
doc_dir = join(expanduser('~'), 'Documents')

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

bas_dir = join(gwfm_dir, 'BAS6')
proj_dir = join(gwfm_dir,'Mapping')
plt_dir = join(proj_dir,'figures/')


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

# from mf_utility import get_layer_from_elev
from map_cln import gdf_bnds, plt_cln


# %%
run_dir = 'C://WRDAPP/GWFlowModel'
run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'


model_ws = loadpth+model_nam


# %%
load_only = ['DIS','BAS6','UPW','OC','SFR','LAK',
            'RCH', 'WEL'
            ]
m = flopy.modflow.Modflow.load('MF.nam', model_ws= model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only
                              )

if 'LPF' in m.get_package_list():
    gel_nam = 'LPF'
else:
    gel_nam = 'UPW'
gel = m.__getattr__(gel_nam)

# %% [markdown]
# # WEL/RCH Updates
# - remove WEL/RCH where new floodplain interaction is

# %%

# %% [markdown]
# # SFR Updates
# Two options:
# 1. Totally replace the streamflow routing with RIV cells with stage based on output from HEC-RAS
# 2. Replace cells at proposed restoration sites with RIV cells
#     - downstream segment needs to divert flow from upstream to skip over the RIV cells
#     - if the RIV is removed in the dry season then the segment data could be redefined to allow flow passing through the segments
#     - likely need to define strhc1 by segment to be able to assign a value of 0 when RIV cells overlap.

# %%

# %%
