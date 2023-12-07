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

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import glob
import sys
import time
from importlib import reload
import h5py

import pandas as pd
import numpy as np

# standard geospatial python utilities
import shapely
import geopandas as gpd

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Regional/'
model_ws = loadpth+'historical_simple_geology_reconnection'


# %%

# %%

# %%

def read_crop_arr_h5(crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "r") as f:
        grp = f['array'] 
        print(grp.keys())
        dset = grp[crop]
        arr = dset[:]
    return(arr)
        
    

# %%
pc_all.shape

# %%
crop='Alfalfa'
crop='Pasture'
crop='Corn'
year=2015
# need separte hdf5 for each year because total is 300MB, group by crop in array
fn = join(model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
pc_all = read_crop_arr_h5(crop, fn)

# # applied water (GW and SW are separate)
# fn = join(model_ws, 'field_SWB', "GW_applied_water_WY"+str(year)+".hdf5")
# irr_gw = read_crop_arr_h5(crop, fn)

# fn = join(model_ws, 'field_SWB', "SW_applied_water_WY"+str(year)+".hdf5")
# irr_sw = read_crop_arr_h5(crop, fn)

# %%
pc_all.shape

# %%
pc_all.shape
