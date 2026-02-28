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
# Goal: take already copied/split up modflow model and fix files that are identified as an issue

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
import geopandas as gpd

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
run_dir = 'F://WRDAPP/GWFlowModel'
# run_dir = 'D://WRDAPP/GWFlowModel'

# loadpth = run_dir +'/Cosumnes/levee_setback/streamflow/'
# # model_nam = 'setback_streamflow'
# model_nam = 'historical_streamflow'

loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology_reconnection'
# model_nam = 'input_write_2014_2020'
model_nam = 'input_write_2014_2022'

base_model_ws = loadpth+model_nam

# %%
load_only = ['DIS','BAS6','SFR']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only,
                              )

nrow,ncol,nlay,delr,delc = (m.dis.nrow, m.dis.ncol, m.dis.nlay, m.dis.delr, m.dis.delc)

# if 'LPF' in m.get_package_list():
#     gel_nam = 'LPF'
# else:
#     gel_nam = 'UPW'
# gel = m.__getattr__(gel_nam)

# %%
from mf_utility import get_dates, get_layer_from_elev, clean_wb
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')


# %% [markdown]
# # Iterate over different economic runs to fix issue file
# There isn't output for the flow budget from OWHM or hob out because I turned them off thinking I wouldn't need to check them.  
# - turns out it was just because the 2014-2022 base run didn't have the manual edit to the bas file. Should try to add this for future runs.
#

# %%
loadpth = run_dir +'/Cosumnes/Economic/'
model_nam = 'input_write_2014_2022'
scen = 'R200'
# it is probably better to create a slightly different file name then to copy these over for a set scenario
econ_model_ws = join(loadpth, model_nam+'_'+scen, 'crop_modflow')

all_run_dates = pd.read_csv(join(econ_model_ws, 'all_run_dates.csv'))



# %%
# for scen in ['R200']:
#     econ_model_ws = join(loadpth, model_nam+'_'+scen, 'crop_modflow')
    
#     # iterate over all economic model ws and years
#     for n, d in enumerate(all_run_dates.date[:-1]):
#         econ_ws_yr =join(econ_model_ws, d)
    
#         # file to change
#         input_file_path = join(econ_ws_yr, m.name+'.bas')
#         output_file_path = join(econ_ws_yr, m.name+'.bas.temp')
#         # text to change/update
#         target_line_content = "FREE"
#         replacement_content = "FREE BUDGETDB flow_budget.txt"
#         # open existing file and create temporary output file to write to
#         with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
#             for n, line in enumerate(infile):
#                 # only make changes in top rows where the specificed text is present
#                 if (n<10)&(target_line_content in line)&(replacement_content not in line):
#                     # Modify the line
#                     outfile.write(line.replace(target_line_content, replacement_content))
#                 else:
#                     # Keep the original line
#                     outfile.write(line)
        
#         # Replace the original file with the new one
#         os.replace(output_file_path, input_file_path)


