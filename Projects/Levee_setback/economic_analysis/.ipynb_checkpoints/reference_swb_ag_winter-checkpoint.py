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
# # Goal
# Run soil water budget for agricultural fields to estimate the winter recharge due to rainfall this script doesn't have to be a function if we update it to read in the crop choice csv then have it save the swb output so the other script can read it in. This script is also important as it will be the sole calculator for deep percolation for fallowed fields which means we need to decide which ETc rates to assign, this should be less important as ET is less important in winter when percolation is happening. This script can also serve the purpose of calculating SWB for Other fields (to get pumping and eep percolation) if desired or this can be imported from the original model.

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
# year = 2018
# # m_nam = 'input_write_2014_2020'
# m_nam = 'input_write_2014_2022_R203'
# loadpth = 'F:/WRDAPP/GWFlowModel/Cosumnes/Economic'
# init_wc = None
# wc_df = pd.read_csv(join(swb_out, 'swb_water_content_'+str(year-1)+'.csv'))

# plot_bool = False

# %%
# init_wc = pd.read_csv('temp_wc.csv')

# %%
def run_swb_ag_winter(year, 
                      m_nam = 'input_write_2014_2020', loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Economic',
                     init_wc = None, plot_bool=False):
    ''' 
    function to run the SWB for the agricultural fields in the winter
    needs some updates to better represent initial conditions (sample
    water content from last timestep of optimized model) and
    to use Kc from the economic model
    INPUT:
    year which is interpreted as water year
    m_nam: basename of model workspace that is used
    loadpth: the complete file path of the directory that holds the model workspace
    init_wc: defines intial water content conditions, should be an array the same shape as soil_ag_all (UniqueID)
            it also specifies the date that it occurs on
            if init_wc is None then the wilting point is assumed as initial conditions
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
    # import functions.swb_functions
    # reload(functions.swb_functions)
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
    
    model_ws = join(loadpth, m_nam)
    # identifies which fields are which crop
    swb_ws = join(model_ws, 'rep_crop_soilbudget')
    # should hold the output from this SWB
    swb_out = join(model_ws, 'crop_soilbudget','ag_non_irr')
    os.makedirs(swb_out, exist_ok=True)
    
    crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    crop_in = crop_in.rename(columns={'parcel_id':'UniqueID'})



    # %%
    soil_ag_all = pd.merge(soil_ag_in, CN)
    
    # add crop data to fields, NAs are introduced for when a crop isn't specified
    soil_ag_all = soil_ag_all.merge(crop_in)
    # for now drop fields without crop specified 
    soil_ag_all = soil_ag_all[~soil_ag_all.name.isna()]

    # %%
    avg_irr_eff = pd.read_csv(join(proj_dir, 'model_inputs', 'avg_irr_eff_by_crop.csv'), index_col=0)


    # %%
    # reference for field to grid cell
    grid_soil = pd.read_csv(join(soil_path, field_ids+'_field_to_cell.csv'), index_col=0)

# %% [markdown]
#     # We should be able to run all fields at once since there is no interdependency or optimization required.
#     # - for the irrigated fields in the winter we need to run a code like this that updates based on crops.

# %% [markdown]
# The start date and end date are uniform across all fields but we'll need to make a bit more flexible to account for varying end dates of irrigation when passing the final water content to start this. Misc grain and hay ends in July, corn ends in September. How did I choose the WY range for this? It should run for the same period of the model from (year)-4-1 to (year+1)-3-31

    # %%
    # strt_date = str(year-1)+'-10-1'
    # end_date = str(year)+'-9-30'
    strt_date = str(year)+'-4-1'
    end_date = str(year+1)+'-3-31'
    dates = pd.date_range(strt_date, end_date)
    print('\n\n', 'Running simple SWB for ag fields to provide gap filling from', strt_date, end_date,
          'data is only used to fill in gaps around irrigation SWB, i.e., mostly winter')

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
    # remove Mixed deciduous and safflower but the others all follow a fairly consistent trend/magnitude
    other_et_cat = ['Strawberries','Rice','Tomatoes (processing)',
                    'Melons, squahs and cucumbers (all types',
                    'Miscellaneous truck','Flowers nursery & Christmas tree farms']
    # use the subset categories to represent other
    rep_ETc = ETc_all[other_et_cat].copy()
    # use average ET to represent
    ETc_all['Other'] = rep_ETc.mean(axis=1)
    # repeat the same for irrigation efficiency
    avg_irr_eff.loc['Other'] = avg_irr_eff.loc[other_et_cat].mean()

    # %%
    # for the crop model we can use standard Kc for default
    # should update with specific Kc except for Fallow pasture
    crop_match = pd.read_csv(join(proj_dir,'data','crop_name_matching.csv'))
    crop_dict = crop_match.set_index('Crop_original')['Crop_model'].to_dict()
    # add value for Other that is made on the fly 
    crop_match.loc[crop_match.index.max()+1] = ['Other','Other']
    # extract relevant ET
    ETc_model = ETc_all[crop_match.Crop_original.values]
    # rename crops for the model
    ETc_model.columns = crop_match.Crop_model.values
    # grab the irrigation efficiency multiple for crop of interest
    irr_eff_mult_df = 100/avg_irr_eff.loc[crop_match.Crop_original.values].Avg_eff
    # add value for fallow to ensure no irrigation
    irr_eff_mult_df.index =  crop_match.Crop_model.values
    irr_eff_mult_df.loc['Unclassified fallow'] = 0



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

    # sample the irrigation efficiency to the fields
    irr_eff_mult = irr_eff_mult_df.loc[soil_ag_all.name.values]


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
    # if the file doesn't exist then assign value of None which function assumes is field water content
    # load water content from winter SWB module 
    wc_fn = join(model_ws, 'crop_soilbudget', 'ag_non_irr',
                                  'swb_water_content_'+str(year-1)+'.csv')
    if exists(wc_fn):
        wc_df = pd.read_csv(wc_fn,index_col=0, parse_dates=True)
        wc_df.columns = wc_df.columns.astype(int)
        # sample the initial water content for the fields corresponding to the crop
        crop_wc = wc_df.loc[:,soil_ag_all['UniqueID'].values].copy()
    #     # sample the water content for the date that irrigation run starts
        crop_wc_init = crop_wc.iloc[-1].copy().values
    else:
        # the first year that uses the fall period from MODFLOW needs to assume initial conditions of field water content
        print('Assuming field water content for initial conditions')
        crop_wc_init = None

    # %%
    # # initial conditions should technically be a mix of last years ag winter + irrig
    # # identify the day of the run period to index and insert the final water content
    # # from the irrigation SWB resutls
    init_wc['wc_index'] = (pd.to_datetime(init_wc.end_date)-pd.to_datetime(strt_date)).dt.days
    # create dataframe that has the correct index in the parcels and in the timestep
    init_wc_ref = soil_ag_all[['UniqueID']].merge(init_wc[['parcel_id','init_wc','wc_index']].rename(columns={'parcel_id':'UniqueID'}), how='left')
    init_wc_ref['swb_index'] = init_wc_ref.index
    init_wc_ref = init_wc_ref.dropna(subset='init_wc')
    init_wc_ref.wc_index = init_wc_ref.wc_index.astype(int)



    # %%
    nfield = soil.nfield
    nper = ETc.shape[0]
    
    # m2_ac = (1/0.3048**2)/43560 # convert from m2 to acres
    
    wc = np.zeros((nper+1, nfield)) # water content, add initial conditions with +1
    if crop_wc_init is None:
        # initial conditions of wilting point water content after irrigation season
        # wc[0,:] = soil_ag_all.w15bar/100
        # might try field content again as recharge estimates are a bit low, as a test to see sensitivity
        # if start date is winter then it should be field capacity
        wc[0,:] = soil_ag_all.w3rdbar/100

    elif crop_wc_init is not None:
        wc[0,:] = crop_wc_init
        
    # the initial water content should be sampled from the final day of the optimized SWB
    pc = np.zeros((nper, nfield)) # percolation
    rp = np.zeros((nper, nfield)) # runoff 
    ETa = np.zeros((nper, nfield)) # actual ET
    AW = np.zeros((nper,nfield)) # really only need one simple value for this case as we don't select irrig source
    soildepth = soil_ag_all.SoilDepth.values

    # %%
    for ns, n in enumerate(np.arange(-1, nper-1)):
    # for ns, n in enumerate([-1, 0, 1]):
        # update the water content with the final value from the optimized irrigation when it is available
        # should not be overwritten because it references ns rather than ns+1 and calculation of ns+1 wc references ns
        fill_wc = init_wc_ref[init_wc_ref.wc_index==ns]
        # we shouldn't need to specify the index here for time since it is determined above
        wc[fill_wc.wc_index, fill_wc.swb_index] = fill_wc.init_wc
        # after end date of opt, no irrigation should be applied either
        stop_aw = init_wc_ref[init_wc_ref.wc_index<=ns]

        ## Runoff ##
        S = calc_S(wc[ns+1], soil.Smax, soil.wc_f, soil.por)
        water_in = rain[n+1] 
        # calculate runoff only when there is rain, and rain is uniform
        if (water_in>0).any():
            rp[n+1] = ((water_in - 0.2*S)**2)/(water_in + 0.8*S)
        # where rainfall is less than initial abstraction (0.2S) set runoff as 0
        rp[n+1] = np.where(water_in<0.2*S, 0, rp[n+1])
        # use effective preciptation to calculate irrigation storage as
        # farmers often do not actually check soil water storage in reality
        irr_stor = soildepth*np.where(water_in - rp[n+1] <0, 0, water_in - rp[n+1])
        # add irrigation with the assumption it is designed to limit runoff
        irr_req = np.where(ETc[n+1] - irr_stor<0, 0, ETc[n+1] - irr_stor)
        AW[n+1] = (irr_req)*irr_eff_mult
        # applied water is zero after irrigation season ends for optimized crops
        AW[n+1, stop_aw.swb_index] = 0
        # add applied water to soil zone
        wc[ns+1] = wc[ns+1] + AW[n+1]/soildepth
        water_in = water_in + AW[n+1]
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
    if plot_bool:
        import matplotlib.pyplot as plt
        # identify subset of fields to use for averaging in the plot
        crop = 'Vineyards'
        # crop = 'Mixed pasture'
        # crop = 'Unclassified fallow'
        plt_field = soil_ag_all.name==crop
        
        fig,ax = plt.subplots(3,1)
        # plt.plot(water_in, label='water in')
        fig.suptitle(crop)
        ax[0].plot(rain, label='rain')
        ax[0].plot(rp[:,plt_field].mean(axis=1), label='rp')
        ax[0].plot(ETa[:,plt_field].mean(axis=1), label='ETc')
        ax[1].plot(pc[:,plt_field].mean(axis=1), label='pc')
        
        # ax[1].set_ylabel('Percolation')
        ax[1].plot(AW[:,plt_field].mean(axis=1), label='AW')
        ax[-1].plot(wc[:,plt_field].mean(axis=1), label='wc')
        ax[-1].set_ylabel('WC')
        # runoff potential follows rain at a lesser degree
        ax[0].legend()
        ax[1].legend()
        
        # plt.close()

    # %%
    # convert percolation to dataframe for output
    pc_df = pd.DataFrame(pc, columns=soil_ag_all.UniqueID.values)
    pc_df.index = dates
    # need to export water content to define initial conditions for next year
    wc_df = pd.DataFrame(wc[1:], columns=soil_ag_all.UniqueID.values)
    wc_df.index = dates
    # runoff is just in case we decide to use later
    rp_df = pd.DataFrame(rp, columns=soil_ag_all.UniqueID.values)
    rp_df.index = dates

    # %%
    # save components to csv for reference if needed
    pc_df.to_csv(join(swb_out, 'ag_non_irr_swb_percolation_'+str(year)+'.csv'))
    wc_df.to_csv(join(swb_out, 'swb_water_content_'+str(year)+'.csv'))
    rp_df.to_csv(join(swb_out, 'swb_runoff_'+str(year)+'.csv'))
    return pc_df
