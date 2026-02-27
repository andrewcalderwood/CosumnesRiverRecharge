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

model_nam = 'historical_simple_geology_reconnection'
model_nam = 'input_write_2014_2020'
# model_nam = 'input_write_2014_2022'

base_model_ws = loadpth+model_nam

# %%
load_only = ['DIS','BAS6','SFR']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              # load_only = load_only,
                              )

nrow,ncol,nlay,delr,delc = (m.dis.nrow, m.dis.ncol, m.dis.nlay, m.dis.delr, m.dis.delc)

if 'LPF' in m.get_package_list():
    gel_nam = 'LPF'
else:
    gel_nam = 'UPW'
gel = m.__getattr__(gel_nam)

# %%
from mf_utility import get_dates, get_layer_from_elev, clean_wb
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')


# %%
# I've only ran the 2014-2020 model now, switch to 2014-2022 later

wb, out_cols, in_cols = clean_wb(base_model_ws, dt_ref)
# # clean_wb?

# %%
fig,ax= plt.subplots(3,1, sharex=True)
ax2 = ax[0].twinx()
wb.plot(y='PERCENT_ERROR', ax=ax2)
ax2.set_ylabel('Percent Error')
wb.plot(y=['dSTORAGE_sum'], legend=False, ax=ax[0], color='black')
ax[0].set_ylabel('Cumulative Storage \nChange ($m^3$)')

wb.plot(y=out_cols, ax=ax[1], legend=True)
wb.plot(y=in_cols, ax=ax[2], legend=True)


# %% [markdown]
# # Iterate over different economic runs
# There isn't output for the flow budget from OWHM or hob out because I turned them off thinking I wouldn't need to check them. 
# The models take too long to to turn back on and re-run so I should plan to just run zonebudget with the cbc file for periods available and check heads at select locations. Maybe start by sampling heads because this might be faster if I choose just a few representative points. -> look at code that extracts simulated time series for HOB wells in oneto denier stream seepage

# %%
loadpth = run_dir +'/Cosumnes/Economic/'
model_nam = 'input_write_2014_2022'
scen = 'R203'
# it is probably better to create a slightly different file name then to copy these over for a set scenario
econ_model_ws = join(loadpth, model_nam+'_'+scen, 'crop_modflow')

all_run_dates = pd.read_csv(join(econ_model_ws, 'all_run_dates.csv'))



# %%
# load_only = ['DIS']
# m_per = flopy.modflow.Modflow.load('MF.nam', model_ws= join(econ_model_ws,d), 
#                                 exe_name='mf-owhm', version='mfnwt',
#                               load_only = load_only,
#                               )
# it looks like the start date time isn't correct in the economic model ws so can't really use this
# strt_date_per, end_date_per, dt_ref_per2 = get_dates(m_per.dis, ref='strt')


# %%
# for n, d in enumerate(all_run_dates.date[:-3]):
#     # create locally referenced dt_ref file to avoid reloading dis
#     dt_ref_per = dt_ref[(dt_ref.dt<=all_run_dates.date[n+1])&(dt_ref.dt>=d)].copy()
#     dt_ref_per.kstpkper = list(zip(np.zeros(len(dt_ref_per),dtype=int), np.arange(len(dt_ref_per))))

#     econ_ws_yr =join(econ_model_ws, d)
#     wb_per, out_cols_per, in_cols_per  = clean_wb(econ_ws_yr, dt_ref_per)

# econ_ws_yr

# %%
# wb, out_cols, in_cols = clean_wb(base_model_ws, dt_ref)
# clean_wb?
