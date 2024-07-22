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
# # Run soil water budget with no irrigation for all fields as a reference.

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
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

# # other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

import Basic_soil_budget_monthly as swb
reload(swb)

# %%
import swb_functions
reload(swb_functions)
# from swb_functions import prep_soil_dict, calc_S, calc_pc
# from swb_functions import calc_yield, calc_profit, choose_water_source
# from swb_functions import run_swb, mak_irr_con, mak_irr_con_adj
from swb_functions import run_swb_model

# %%
from functions.data_functions import init_h5, crop_arr_to_h5

# %%


# so if you have a dictionary d and want to access (read) its values with the syntax x.foo instead of the clumsier d['foo'], just do
# convert a dictionary to an object with object style referencing
class cost_variables(object):
  def __init__(self, adict):
    self.__dict__.update(adict)


# %%

# %%
field_ids = 'parcel'
field_ids = 'native'

soil_path = join(uzf_dir,'clean_soil_data')
# soil data for each ag field
soil_ag_all = pd.read_csv(join(soil_path, 'soil_for_'+field_ids+'_fields.csv'), index_col=0)

# curve numbers
CN = pd.read_csv(join(soil_path, field_ids+'_field_CN.csv'),index_col = 0)


# %%
def no_irr_swb(soil, rain, ETc)
    """ 
    Run simple SWB with no irrigation
    """
    nfield = soil.nfield
    nper = ETc.shape[0]
    
    irr_sw = np.zeros((nper,nfield))
    irr_gw = np.zeros((nper,nfield))
    # run the time series behind the soil water budget model
    pc, K_S = run_swb_model(soil, rain, ETc, irr_gw, irr_sw)
    return(pc)


# %%
soil_crop = swb.load_soil(pred_dict[crop], crop_in)
soil_ag = soil_crop.iloc[[ns]] #keep as dataframe for consistency 
soil_dict = prep_soil_dict(soil_ag, etc_arr, var_crops)
soil = cost_variables(soil_dict)

# %% [markdown]
# We should be able to run all fields at once since there is no interdependency or optimization required.

# %%
no_irr_swb(soil, rain, ETc)
