# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Goal: Convert the final output into a more useable format for post-processing

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

import matplotlib.pyplot as plt

# %%
import seaborn as sns

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# %%
add_path('../')


# %%
from functions.output_processing import get_wb_by_parcel
from functions.f_gw_dtw_extract import sample_dtw, avg_heads
import functions.Basic_soil_budget_monthly as swb

# initialize HDF5 files for the year
from functions.data_functions import init_h5


# %%
# from functions.swb_functions import run_swb
# import f_rep_swb_profit_opt
# reload(f_rep_swb_profit_opt)

from f_rep_swb_profit_opt import load_run_swb


# %%
def read_crop_arr_h5(crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "r") as f:
        grp = f['array'] 
        # print(grp.keys())
        dset = grp[crop]
        arr = dset[:]
    return(arr)


# %%
# for batch file case we want to ignore the consant h5py deprecation warning
import warnings

warnings.filterwarnings("ignore")
# warnings.filterwarnings("ignore", category=DeprecationWarning) 

# %%
# updated version specifies concept_name and copy_files here so it can be easily
# seen as these are the main update to make in a script
in_data = sys.argv

if 'ipykernel' in in_data[0]:
    # m_nam = 'input_write_2014_2020_R3'
    # m_nam = 'input_write_2014_2020'
    m_nam = 'input_write_2014_2025_R204'
    # m_nam = 'input_write_2016_2022_R200'

    # input_name = 'static_model_inputs.xlsx'
    input_name = 'static_model_inputs_no_p_o.xlsx'

else:
    m_nam = in_data[1]
    # added option to specify different static_model_inputs to allow easier adjustment of operating costs, revenue
    input_name = in_data[2]


print('sys.argv[1] (m_nam) is...')
print(m_nam)

print('sys.argv[2] (input_name) is...')
print(input_name)

t_start = time.time()

# %%
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Economic'
loadpth = 'D://WRDAPP/GWFlowModel/Cosumnes/Economic'

# update to different modflow models here
# m_nam = 'input_write_2014_2020'
# m_nam = 'input_write_2014_2022'

# scenario = '_R20' # pumping constraint
# scenario = '_R4' # 90/20 floodplain
# # scenario = '_R3' # 6x existing diversion for MAR vineyard

# scenario = '_R200' # 200 represents no p_o
# scenario = '_R203' # no p_o and 6x MAR
# scenario = '_R204' # no p_o and 90/20 flloodplain

# scenario=''
# m_nam = m_nam+scenario


model_ws = join(loadpth, m_nam)


# %%
# provide representative soil water budget folder
swb_ws = join(model_ws, 'rep_crop_soilbudget')
# os.makedirs(join(swb_ws, 'output'), exist_ok=True)

out_dir = join(model_ws, 'output_clean')
os.makedirs(out_dir, exist_ok=True)


# %%
# save log by date so we can see old versions
os.makedirs(join(out_dir, 'log'), exist_ok=True)
if 'ipykernel' in in_data[0]:
    print('print to jupyter')
else:
    print('printing remaining output to log file')
    sys.stdout = open(join(out_dir, 'log', 'plot_output_only_log_'+str(pd.to_datetime('today').date())+'.txt'), 'w')


# %%
# load parcel data for reference as needed
parcels = gpd.read_file(join(proj_dir,'Parcels shapefile/parcels.shp'))
parcels['area_m2'] = parcels.geometry.area
parcels['acres'] = parcels.area_m2/(43560*0.3048**2)
parcels.UniqueID = parcels.UniqueID.astype(int)

# %%
all_run_dates = pd.read_csv(join(model_ws, 'crop_modflow', 'all_run_dates.csv'), parse_dates=['date'])
# years to sample output
run_years = all_run_dates[all_run_dates.use=='irrigation'].date.dt.year.values
# run_years = run_years[:-1]

# %%
# check which files are available and shorten as needed
files = pd.DataFrame(glob.glob(join(out_dir, 'profit_yield_long_*.csv')), columns=['fn'])
files['name'] = files.fn.apply(basename)
files['year'] = files.name.str.extract(r'(\d{4})').astype(int)

# subset run years to those with output
run_years = run_years[pd.Series(run_years).isin(files.year)]
# should set up a way to delete old output to avoid mixing old/new in plots

# %%
# temp for while model is still running
# run_years = np.arange(2001, 2011)
# run_years = np.arange(2001, 2020)

# %% [markdown]
# # Plot SWB results by dtw


# %%
y=2015
irr_gw = pd.read_csv(join(swb_ws, 'output','irr_gw_all'+str(y)+'.csv'),index_col=0)
