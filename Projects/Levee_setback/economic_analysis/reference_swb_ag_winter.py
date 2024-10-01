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
# # Goal
# Run soil water budget for agricultural fields to estimate the winter recharge due to rainfall this script doesn't have to be a function if we update it to read in the crop choice csv then have it save the swb output so the other script can read it in. This script is also important as it will be the sole calculator for deep percolation for fallowed fields which means we need to decide which ETc rates to assign, this should be less important as ET is less important in winter when percolation is happening.

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


# year = 2015

# %%
def run_swb_ag_winter(year):
    ''' 
    function to run the SWB for the agricultural fields in the winter
    needs some updates to better represent initial conditions (sample
    water content from last timestep of optimized model) and
    to use Kc from the economic model
    INPUT:
    year which is interpreted as water year
    OUTPUT:
    A csv with the deep percolation time series for each
    UniqueID agricultural field
    '''
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
    # field_ids = 'native'
    
    soil_path = join(uzf_dir,'clean_soil_data')
    # soil data for each ag field
    soil_ag_in = pd.read_csv(join(soil_path, 'soil_for_'+field_ids+'_fields.csv'), index_col=0)
    
    # curve numbers
    CN = pd.read_csv(join(soil_path, field_ids+'_field_CN.csv'),index_col = 0)

    # %%
    # need to load dataframe which defines native land use type
    # nat_lu = pd.read_csv(join(proj_dir, 'native_parcel_zonalstats','native_land_use.csv'),index_col=0)

    # %%
    loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic'
    
    m_nam = 'historical_simple_geology_reconnection'
    
    model_ws = join(loadpth, m_nam)
    
    swb_ws = join(model_ws, 'rep_crop_soilbudget')
    
    crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    crop_in = crop_in.rename(columns={'parcel_id':'UniqueID'})
    # temporarily replace Other and 0 with unclassified fallow
    crop_in.loc[crop_in.name.isin(['Other','0']),'name'] = 'Unclassified fallow'


    # %%
    soil_ag_all = pd.merge(soil_ag_in, CN)
    
    # add crop data to fields, NAs are introduced for when a crop isn't specified
    soil_ag_all = soil_ag_all.merge(crop_in)
    # for now drop fields without crop specified 
    soil_ag_all = soil_ag_all[~soil_ag_all.name.isna()]

    # %%
    # reference for field to grid cell
    grid_soil = pd.read_csv(join(soil_path, field_ids+'_field_to_cell.csv'), index_col=0)

# %% [markdown]
#     # We should be able to run all fields at once since there is no interdependency or optimization required.
#     # - for the irrigated fields in the winter we need to run a code like this that updates based on crops.

    # %%
    strt_date = str(year-1)+'-10-1'
    end_date = str(year)+'-9-30'
    dates = pd.date_range(strt_date, end_date)

    # %%
    # load CIMIS rain, ETo data m/day
    rain, ETo_df = swb.load_hyd(dates)

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
    # for the crop model we can use standard Kc for default
    # should update with specific Kc except for Fallow pasture
    crop_match = pd.read_csv(join(proj_dir,'data','crop_name_matching.csv'))
    crop_dict = crop_match.set_index('Crop_original')['Crop_model'].to_dict()
    # rename crops for the model
    ETc_model = ETc_all[crop_match.Crop_original.values].rename(columns = crop_dict)


    # %%
    # for fallow lands consider 
    # native pasture, mixed pasture and native vegetation are all the same values (this is likely)
    # an artifact from me ID joining to the closest match
    # grassland, light brush, medium brush have the same pattern as well with decreasing ETc
    # in the summer, light has high Kc than medium/heavy likely because less shade means more ET
    # the question is then is pasture or grassland/brush more representative of fallowed fields
    # fallow_cols1 = ['Native pasture',  'Mixed pasture']
    # fallow_cols2 = ['Native vegetation','Grassland', 'Light brush', 'Medium brush']
    
    # ETc_all[fallow_cols1].plot()
    # ETc_all[fallow_cols2].plot()
    
    # best to just assume pasture as this will most conservatively predict recharge by
    # reducing soil water storage in the summer. Also the Kc of .95 from the Guerra report closely matches the 0.9 from CIMIS 1980s
    
    

    # %%
    # pull out the crop ET for the fields simulated
    # it may be necessary to add ET in the winter months representing bare soil ET?
    ETc = ETc_model.loc[:,soil_ag_all.name.values].values
    ETc[np.isnan(ETc)] = 0


    # %%
    print(ETc_model.columns.values)

    # %%
    # create soil dictionary
    soil_dict = base_soil_dict(soil_ag_all)
    # convert to class for cleaner referencing
    soil = cost_variables(soil_dict)
    # assume baseline depth of 2 meters
    soil.depth = soil_ag_all.SoilDepth.values

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
    # initial conditions of wilting point water content after irrigation season
    wc[0,:] = soil_ag_all.w15bar/100
    # the initial water content should be sampled from the final day of the optimized SWB
    pc = np.zeros((nper, nfield)) # percolation
    rp = np.zeros((nper, nfield)) # runoff 
    ETa = np.zeros((nper, nfield)) # actual ET
    irr_sw = np.zeros((nper,nfield))
    irr_gw = np.zeros((nper,nfield))
    soildepth = soil_ag_all.SoilDepth.values

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
    import matplotlib.pyplot as plt

    # %%
    fig,ax = plt.subplots(3,1)
    # plt.plot(water_in, label='water in')
    ax[0].plot(rain, label='rain')
    ax[0].plot(rp.mean(axis=1), label='rp')
    ax[0].plot(ETa.mean(axis=1), label='ETc')
    ax[1].plot(pc.mean(axis=1), label='pc')
    ax[1].set_ylabel('Percolation')
    ax[-1].plot(wc.mean(axis=1), label='wc')
    ax[-1].set_ylabel('WC')
    # runoff potential follows rain at a lesser degree
    ax[0].legend()
    plt.close()

    # %%
    # convert percolation to dataframe for output
    pc_df = pd.DataFrame(pc, columns=soil_ag_all.UniqueID.values)
    pc_df.index = dates
    
    

    # %%
    pc_df.to_csv(join(model_ws,
                      'ag_non_irr_swb_percolation_'+str(year)+'.csv'))
    return pc_df
