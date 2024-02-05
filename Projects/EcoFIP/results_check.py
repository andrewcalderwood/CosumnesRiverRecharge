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
from os.path import basename, dirname, join, exists
import sys
from importlib import reload
import glob
import pandas as pd
import numpy as np
import calendar
import time

# standard python plotting utilities
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

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
from matplotlib.ticker import MaxNLocator

# import flopy
# import flopy.utils.binaryfile as bf
from importlib import reload
import h5py



# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

bas_dir = join(gwfm_dir, 'BAS6')
proj_dir = join(gwfm_dir,'EcoFIP')
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

from map_cln import gdf_bnds, plt_cln
from report_cln import base_round

# %%
# the .to_hdf format for format='fixed' means that the name of each column is saved separate from the values
# e.g., block0_items is column name and block0_values has the values
out_dir = join(proj_dir,'output')
with h5py.File(join(out_dir, 'sfrdf_mon_all.hdf5')) as f:
    grp = f['monthly']
    cols = grp['axis0'][:]
    rows = grp['axis1'][:]
    n1 = grp['block0_items'][:]
    v1 = grp['block0_values'][:]
    print(grp.keys())

# %%
# n1, v1

# %%
sfrdf_mon = pd.read_hdf(join(out_dir, 'sfrdf_mon_all.hdf5'))

# %%
sfrdf_mon.to_hdf(join(doc_dir, 'sfrdf_test.hdf5'), 
                                  key='monthly', complevel=4, index=False,
                 data_columns = sfrdf_mon.columns, format='table', mode='w')

# %%
with h5py.File(join(doc_dir, 'sfrdf_test.hdf5')) as f:
    grp = f['monthly']
    tbl = grp['table'][:]
    print(grp.keys())

# %%
pd.DataFrame(tbl)

# %%
