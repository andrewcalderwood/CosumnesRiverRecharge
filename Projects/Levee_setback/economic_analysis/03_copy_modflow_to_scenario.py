# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Goal: take already copied/split up modflow model and copy to the scenario directory

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename
import sys
import glob
from importlib import reload
import shutil

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
# import rasterio

# # mapping utilities
# import contextily as ctx
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as fm

# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
gwfm_dir


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

from mf_utility import get_layer_from_elev
from map_cln import gdf_bnds



# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
run_dir = 'C:/WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
# run_dir = 'D://WRDAPP/GWFlowModel'

# loadpth = run_dir +'/Cosumnes/levee_setback/streamflow/'
# # model_nam = 'setback_streamflow'
# model_nam = 'historical_streamflow'

loadpth = run_dir +'/Cosumnes/Regional/'

model_nam = 'historical_simple_geology_reconnection'
model_nam = 'input_write_2014_2020'
model_nam = 'input_write_2014_2022'
# model_nam = 'input_write_2000_2022'

base_model_ws = loadpth+model_nam


# %%
# # for files that don't need updates
# # they can be diretcly copied to save file formatting write time
# # # copy over basic package to enable flopy to read in the model, will update start heads later
# # ,'MF.evt',
# files_copy = ['MF.lak','MF.bath', 'MF.gage', 'MF.upw','MF.nwt', 'MF.bas']

# for m_per in np.arange(0, all_run_dates.shape[0]-1):
#     m_strt = all_run_dates.iloc[m_per].date
#     model_ws = join(m.model_ws, str(m_strt.date()))
#     for f in files_copy:
#         shutil.copy(join(base_model_ws,f), join(model_ws,f))

# %%

# %%
# need to identify the other files to copy

# %%
# # name file had HOB manually removed
# files_copy = ['MF.nam']

# for m_per in np.arange(0, all_run_dates.shape[0]-1):
#     m_strt = all_run_dates.iloc[m_per].date
#     model_ws = join(m.model_ws, str(m_strt.date()))
#     for f in files_copy:
#         shutil.copy(join(m.model_ws,f), join(model_ws,f))
