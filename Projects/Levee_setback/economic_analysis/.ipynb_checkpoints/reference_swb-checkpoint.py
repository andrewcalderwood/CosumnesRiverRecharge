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

import functions.Basic_soil_budget_monthly as swb
reload(swb)

# %%
import functions.swb_functions
reload(functions.swb_functions)
# from swb_functions import prep_soil_dict, calc_S, calc_pc
# from swb_functions import calc_yield, calc_profit, choose_water_source
# from swb_functions import run_swb, mak_irr_con, mak_irr_con_adj
from functions.swb_functions import run_swb_model, base_soil_dict

# %%
from functions.data_functions import init_h5, crop_arr_to_h5


# %%
# so if you have a dictionary d and want to access (read) its values with the syntax x.foo instead of the clumsier d['foo'], just do
# convert a dictionary to an object with object style referencing
class cost_variables(object):
  def __init__(self, adict):
    self.__dict__.update(adict)


# %%
field_ids = 'parcel'
field_ids = 'native'

soil_path = join(uzf_dir,'clean_soil_data')
# soil data for each ag field
soil_ag_in = pd.read_csv(join(soil_path, 'soil_for_'+field_ids+'_fields.csv'), index_col=0)

# curve numbers
CN = pd.read_csv(join(soil_path, field_ids+'_field_CN.csv'),index_col = 0)

# %%
# need to load dataframe which defines native land use type
nat_lu = pd.read_csv(join(proj_dir, 'native_parcel_zonalstats','native_land_use.csv'))

# %%
soil_ag_all = pd.merge(soil_ag_in, CN)

# add crop data to fields, NAs are introduced for when a crop isn't specified
soil_ag_all = soil_ag_all.merge(nat_lu)
# for now drop fields without crop specified 
soil_ag_all = soil_ag_all[~soil_ag_all.name.isna()]

# %%
# reference for field to grid cell
grid_soil = pd.read_csv(join(soil_path, field_ids+'_field_to_cell.csv'), index_col=0)


# %%
def no_irr_swb(soil, rain, ETc):
    """ 
    Run simple SWB with no irrigation
    """
    nfield = soil.nfield
    nper = ETc.shape[0]
    print('nfield', nfield)
    print('nper', nper)

    irr_sw = np.ones((nper,nfield))*100
    irr_gw = np.zeros((nper,nfield))
    # run the time series behind the soil water budget model
    pc, K_S = run_swb_model(soil, rain, ETc, irr_gw, irr_sw, calc_Ks=False)
    return(pc)


# %% [markdown]
# We should be able to run all fields at once since there is no interdependency or optimization required.

# %%
strt_date = '2014-10-1'
end_date = '2021-9-30'
dates = pd.date_range(strt_date, end_date)

# %%
# load CIMIS rain, ETo data m/day
rain, ETo_df = swb.load_hyd(dates)

# %%
plt.plot(rain)

# %%
# load in pre-processed array of ETc for all time, m/day
ETc_long = pd.read_hdf(join(uzf_dir, "dwr_ETc",'long_ETc_all_lu.hdf5'), key='variable')
# convert to column format
ETc_all = ETc_long.pivot(columns='variable', values='value')

# subset for model period
ETc_all = ETc_all[strt_date:end_date]
# fill in empty dates with interpolation (temporary since corn dates are slightly different)
ETc_all = ETc_all.reindex(dates)

# %%
# pull out the crop ET for the fields simulated
# it may be necessary to add ET in the winter months representing bare soil ET?
ETc = ETc_all.loc[:,soil_ag_all.name.values].values
ETc[np.isnan(ETc)] = 0


# %%
ETc.shape

# %%
# create soil dictionary
soil_dict = base_soil_dict(soil_ag_all)
# convert to class for cleaner referencing
soil = cost_variables(soil_dict)
# assume baseline depth of 2 meters
soil.depth = 2
soil.depth = soil_ag_all.SoilDepth.values

# %%
soil_dict['Ks']

# %%
# soil_crop = swb.load_soil(pred_dict[crop], crop_in)
# soil_ag = soil_crop.iloc[[ns]] #keep as dataframe for consistency 
# soil_dict = prep_soil_dict(soil_ag, ETc, var_crops)


# %%
from functions.swb_functions import calc_S, calc_pc

# %%
nfield = soil.nfield
nper = ETc.shape[0]

# m2_ac = (1/0.3048**2)/43560 # convert from m2 to acres

wc = np.zeros((nper+1, nfield)) # water content, add initial conditions with +1
pc = np.zeros((nper, nfield)) # percolation
rp = np.zeros((nper, nfield)) # runoff 
ETa = np.zeros((nper, nfield)) # actual ET
irr_sw = np.zeros((nper,nfield))
irr_gw = np.zeros((nper,nfield))
soildepth=2
soildepth = soil_ag_all.SoilDepth.values

# %%
# plt.plot(water_in)

# %%
for ns, n in enumerate(np.arange(-1, nper-1)):
# for ns, n in enumerate([-1, 0, 1]):
    ## Runoff ##
    S = calc_S(wc[ns+1], soil.Smax, soil.wc_f, soil.por)
    water_in = rain[n+1] 
    # calculate runoff only when there is rain, and rain is uniform
    if (water_in>0).any():
        rp[n+1] = ((water_in - 0.2*S)**2)/(water_in + 0.8*S)
    # where rainfall is less than initial abstraction (0.2S) set runoff as 0
    rp[n+1] = np.where(water_in<0.2*S, 0, rp[n+1])
    # add in irrigation after runoff (assume farm is set up to avoid runoff for irrigation season)
    water_in = water_in + irr_sw[n+1] + irr_gw[n+1]
    ## explicit percolation ##
    pc[n+1] = calc_pc(wc[ns], soil.por, soil.Ks, soil.m)
    # stepwise water budget, explicit to avoid iteration
    # add rain and take away runoff first
    wc[ns+1] = (wc[ns]*soildepth + (water_in - rp[n+1]))/soildepth
    # take away ET, add term to prevent going to zero
    ETa[n+1] = np.where(ETc[n+1] <= wc[ns+1]*soildepth, ETc[n+1], wc[ns+1]*soildepth - 1E-9)
    wc[ns+1] = wc[ns+1] + (-ETa[n+1])/soildepth
    # take away percolation
    pc[n+1] = np.where(pc[n+1] <= wc[ns+1]*soildepth, pc[n+1], wc[ns+1]*soildepth - 1E-9)
    wc[ns+1] = wc[ns+1] + (-pc[n+1])/soildepth

# %%
wc[ns]

# %%
ns=850
calc_pc(.4, soil.por, soil.Ks, soil.m)
# ns

# %%
import matplotlib.pyplot as plt

# %%

# %%
fig,ax = plt.subplots(3,1)
# plt.plot(water_in, label='water in')
ax[0].plot(rain, label='rain')
ax[0].plot(rp.mean(axis=1), label='rp')
ax[0].plot(ETa.mean(axis=1), label='ETc')
ax[1].plot(pc.mean(axis=1), label='pc')
ax[-1].plot(wc.mean(axis=1), label='wc')

# runoff potential follows rain at a lesser degree
ax[0].legend()

# %% [markdown]
# The model results seem to suggest water content is simply not filling enough to support percolation. A quick check shows once water content hits 0.4 then peroclation is quite high.
# - i had defaulted soildepth to 2 m which appears to have held excessive soil water storage capacity to allow percolation

# %% [markdown]
# From regional model with irrigation:
# - Rainfall+irrigation is usually 0.02 to 0.05
# - Runoff os around 0.01 to 0.02
# - ET is around 0.005 in summer
# - percolation hoavers around 0.001 or a little below
# - VWC is usually 0.2 in winter and 0.3 in wet years

# %%
# this versoin ends up with zero pc
# pc = no_irr_swb(soil, rain, ETc)

# %%

# %%
# convert percolation to dataframe for output
pc_df = pd.DataFrame(pc, columns=soil_ag_all.UniqueID.values)
pc_df.index = dates

pc_df.to_csv( join(proj_dir,'native_parcel_zonalstats', 
                   'swb_percolation.csv'))

# %% [markdown]
# # check results
# Something is broken because we see zero percolation even when I set irrigation to extreme values

# %%
import matplotlib as mpl
import matplotlib.pyplot as plt

# %%
fig,ax = plt.subplots(figsize=(6,3))
ax.plot(rain)
# plt.plot(ETc)
ax.plot(ETc.mean(axis=1))


# %%

plt.plot(pc.mean(axis=1))
