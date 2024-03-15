# ---
# jupyter:
#   jupytext:
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

# %% [markdown]
# # Iterate the irrigation optimizer over crops and years
#
# The first step is running all crops for each year
# For each crop in this version, a series of simulations will be run for different depth to water series
# then the optimal irrigation can be interpolated for each parcel
# Then stepping to the next year

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
from scipy.optimize import minimize
from scipy.optimize import Bounds


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

# # flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# # other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)
from mf_utility import get_layer_from_elev, param_load

import Basic_soil_budget_monthly as swb



# %%
# import swb_functions
# reload(swb_functions)
from swb_functions import prep_soil_dict, calc_S, calc_pc, calc_yield, calc_profit, run_swb, mak_irr_con


# %%
def init_h5(h5_fn):
        """ Initiate hdf5 files for the given year before appending data for each crop"""
        with h5py.File(h5_fn, "w") as f:
            grp = f.require_group('array') # makes sure group exists
            grp.attrs['units'] = 'meters/day'
            grp.attrs['description'] = 'Rows represent the soil units and columns represent the days in the season'

# these hdf5 files are written at the start so they can be appended to for each year and crop
def crop_arr_to_h5(arr, crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "a") as f:
        grp = f.require_group('array') # makes sure group exists
        # grp.attrs['units'] = 'meters/day'
        # grp.attrs['description'] = 'Each layer of the array is a day in the water year'
        dset = grp.require_dataset(crop, arr.shape, dtype='f', compression="gzip", compression_opts=4)
        dset[:] = arr

# so if you have a dictionary d and want to access (read) its values with the syntax x.foo instead of the clumsier d['foo'], just do
# convert a dictionary to an object with object style referencing
class cost_variables(object):
  def __init__(self, adict):
    self.__dict__.update(adict)


# %%
# resampled ground surface elevation
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
# testing
year = int(2020)
crop='Grape'
# crop='Corn'
crop='Alfalfa'
# crop='Pasture' # will require extra work due to AUM vs hay
# crop = 'Misc Grain and Hay'

# %%
# # testing
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Regional/'
base_model_ws = loadpth + 'crop_soilbudget'
crop_in = pd.read_csv(join(base_model_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'))
dtw_df = pd.read_csv(join(base_model_ws, 'field_SWB', 'dtw_ft_parcels_'+str(year)+'.csv'), 
                     index_col=0, parse_dates=['dt'])
dtw_df.columns = dtw_df.columns.astype(int)
soil_rep = True

# %%


def load_run_swb(crop, year, crop_in, base_model_ws, dtw_df, soil_rep = False):
    ''' 
    Function to import variables related to soil water budget function
    to then run the function in a profit optimizer before saving the results in hdf5 format
    Input:
    crop: character string of shortened crop name
    year: water year
    crop_in: dataframe with parcel_id, (crop) name, pod (point of diversion)
    base_model_ws: directory to save files to
    dtw_df: dataframe with depth to water for each daily step needed for the year
    soil_rep: boolean to identify whether a representative soil profile should be used instead
                of unique soils for each field
    '''
    # %%
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)

    # %%
    #  get the dates for each yield cycle
    yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
    yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
    # get the total extent of the irrigation season (calculation period)
    strt_date = yield_start.min()
    end_date = yield_end.max()
    dates = pd.date_range(strt_date, end_date, freq='D')
    nper = (end_date-strt_date).days +1
    print('Start', strt_date.date(), ', End', end_date.date(),', No. days', nper)
    # not use, base_model_ws is better since all files have the year attached
    # model_ws = join(base_model_ws, crop+'_'+str(strt_date.date()))

    # %%
    for var in ['percolation','GW_applied_water', 'SW_applied_water']:
        name = join(base_model_ws, 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        init_h5(name)

    # %%
    Kc, Kc_dates = swb.load_Kc(year)
    
    Kc_c = Kc[Kc.Crop==crop]
    Kc_dates_c = Kc_dates[Kc_dates.Crop==crop]
    
    Kc_df = pd.DataFrame()
    # for alfalfa we need to iterate over this for each cutting
    # function failed for Misc. grain and hay with date in next year
    for c in Kc_dates_c.Cycle.unique():
        Kc_df_c = swb.get_Kc_dates(Kc_dates_c[Kc_dates_c.Cycle==c], Kc_c)
        Kc_df = pd.concat((Kc_df, Kc_df_c))


    # %%
    rain, ETo_df = swb.load_hyd(year, dates)
    # get the crop ET
    ETc = ETo_df.values*Kc_df.Kc.values

    # %%
    # create dictionary of crop variables and join multiples into list
    var_crops_dict = {}
    for v in var_crops.index.unique():
        var_crops_dict[v] = var_crops[v].tolist()
    
    # create dictionary of general parameters
    gen_dict = {**var_gen.to_dict(), **var_crops_dict}
    
    gen_dict['y_max'] = var_crops[['y_max']].values # y_max may be multiple so keep array structure
    gen_dict['nper'] = (end_date-strt_date).days +1
    
    # days to index for calculating yield impacts, add 1 to include end day
    yield_ind = np.append([0], (yield_end-strt_date).dt.days.values +1)
    gen_dict['yield_ind'] = np.append([0], (yield_end-strt_date).dt.days.values +1)

    # %%
    # create time series of daily yield response factors
    var_yield['dt'] = swb.ymd2dt(year, var_yield.month, var_yield.day, var_yield.year_adj)
    # adjust for dates in the next year
    K_Y = var_yield.set_index('dt').reindex(dates).ffill()
    K_Y = K_Y['K_Y'].values

    # %%
    gap_irr = var_crops['gap_irr'] # Number of days between irrigations
    n_irr = np.floor(len(dates)/gap_irr).astype(int) + 1 # Calculate number of irrigations
    irr_days = np.arange(0, (n_irr*gap_irr-1), gap_irr).astype(int) # Calculate days on which irrigation takes place



    # %%
    gen_dict['K_Y'] = var_yield.set_index('dt')['K_Y'].reindex(dates).ffill().values
    gen_dict['n_irr'] = n_irr
    gen_dict['irr_days'] = np.arange(0, (n_irr*gap_irr-1), gap_irr).astype(int) # Calculate days on which irrigation takes place

# %% [markdown]
#     # ## Test SWB on one field

    # %%
    avg_irr_eff = pd.read_csv(join(proj_dir, 'model_inputs', 'avg_irr_eff_by_crop.csv'), index_col=0)
    irr_eff_mult = 100/avg_irr_eff.loc[crop_dict[crop]].Avg_eff
    gen_dict['irr_eff_mult'] = 100/avg_irr_eff.loc[crop_dict[crop]].Avg_eff
    print(pred_dict[crop], ':',crop)

    soil_crop = swb.load_soil(pred_dict[crop], crop_in)
    if soil_rep:
        # get the average field properties for a representative field
        # with one for a point of diversion and one without
        soil_crop = soil_crop.groupby('pod').mean(numeric_only=True).reset_index()
        # repeat soil_crop data to iterate over different DTW profiles
        soil_crop = pd.concat([soil_crop]*dtw_df.shape[1])
        
    # another option instead of looping over texture classes would be to sort by texture
    # then the data would be ordered to reference the nearest similar texture
    nfield_crop = len(soil_crop)
    print('Num crops:', nfield_crop)


    # %%
    soil_path = join(uzf_dir,'clean_soil_data')
    # connection of ag fields to all grid cells
    grid_soil = pd.read_csv(join(soil_path, 'parcel_field_to_cell.csv'),index_col=0)
    # subset to fields for crop (not used?)
    # grid_crop = grid_soil.merge(soil_crop[['UniqueID']])


    # %%
    # crop_wells = soil_crop[['UniqueID']].merge(parcel_wells)
    # select parcels in the simulation
    # for rep soil only need one dtw profile
    if soil_rep:
        crop_dtw = dtw_df.copy()
    else:
        crop_dtw = dtw_df.loc[:,soil_crop['UniqueID'].values]
    # select dates being simulated
    crop_dtw = crop_dtw.loc[dates].values

# %% [markdown]
#     # ## Iterate over each unique soil condition

# %%

    # convert dictionary of variables to class for easier referencing, constant over different soil
    gen = cost_variables(gen_dict)
    
    bounds = Bounds(lb = 0)


    # %%
    t0 = time.time()
    # save irrigation, function, time
    irr_all =np.zeros((nfield_crop,2*n_irr))
    p_all = np.ones(nfield_crop) # make initial profit negative
    t_all = np.zeros(nfield_crop)
    
    # tol = 0.01
    irr_lvl = np.zeros(2*n_irr); # Initial irrigation values for optimization
    irr_lvl[:] = (2/12)*0.3048 # irrigate with 2 inches (convert to meters)
    irr_lvl_base = np.copy(irr_lvl)
    
    etc_arr = np.zeros((nper))
    for n in np.arange(nper):
        etc_arr[n] = ETc[n]
    
    for ns in np.arange(0,nfield_crop):
    # for ns in np.arange(0,10):
        soil_ag = soil_crop.iloc[[ns]] #keep as dataframe for consistency 
        nfield = soil_ag.shape[0]
    
        # dtw_arr = dtw_all[:,ns]
        dtw_arr = crop_dtw[:,ns]
    
        # prep_soil(soil_ag, etc_arr, var_crops)
        soil_dict = prep_soil_dict(soil_ag, etc_arr, var_crops)
        soil = cost_variables(soil_dict)
        
        if ns > 1:
            irr_lvl[:] = irr_all[ns-1]
        # if no POD then no SW irrig
        if soil_ag.pod.iloc[0]=='No Point of Diversion on Parcel':
            irr_lvl[:n_irr] = 0
            irr_lvl[n_irr:] *= 2 # put double the irrigation to the GW
    
        linear_constraint = mak_irr_con(soil_ag, n_irr) # will change for POD
        tol = 0.01
        # continue optimizing until profit is positive
        while p_all[ns] >0 :
            # the minimization with 'trust-constr' and no constraints doesn't solve and has increasing WB error
            out = minimize(run_swb, irr_lvl, args = (soil, gen, rain, ETc, dtw_arr), method='trust-constr',
                    constraints = [linear_constraint],
                    bounds=bounds,
            #          options={'verbose':1}
                           tol = tol
                    )
            # decrease tolerance and reset starting irrigation to help solving
            if out.fun >0:
                tol /=10
                irr_lvl[:] = np.copy(irr_lvl_base)
            if tol < 1E-5:
                break # if tolerance gets too small then skip
            irr_all[ns] = out.x
            p_all[ns] = out.fun
            t_all[ns] = out.execution_time
        # print final results
        print('Soil ', str(ns),'%.2f' %(-out.fun),'$ ,in %.2f' %(out.execution_time/60),'min')
    t1 = time.time()
    print('Total time was %.2f min' %((t1-t0)/60), 'for', ns,'parcels')


# %% [markdown]
#     # It wasn't until running grapes which have 3 times the number of irrigations that I realized that each solver takes about 2 min instead of 0.2 min (Corn). Alfalfa had run times of 0.3 min. Misc. grain and hay never found positive profit  (-250 to -300), and took multiple minutes as well.

# %% [markdown]
#     # # Post-processing
#     # 1. run optimized swb with irrigation efficiency
#     # 2. save applied water for pumping
#     # 3. save deep percolation

    # %%
    # runs relatively quickly
    
    # scale by irrigation efficiency of the crop after optimizing
    irr_true = irr_all * irr_eff_mult # one efficiency for each crop type
    # irr_true[0]
    p_true = np.copy(p_all)
    pc_all = np.zeros((nfield_crop, gen.nper))
    for ns in np.arange(0,nfield_crop):
    # for ns in np.arange(0,10):
        # p_true[ns] = run_swb(irr_true[ns], soil, gen, rain, ETc, dtw_arr)
        p_true[ns], pc, K_S  = run_swb(irr_true[ns], soil, gen, rain, ETc, dtw_arr, arrays=True)
        pc_all[ns] = pc[:,0] # original shape meant for multiple fields, but only has one since iteration is over fields

    # %%
    # break down irrigation into groundwater and surface water time series
    irr_sw_out = np.zeros((nfield_crop, gen.nper))
    irr_sw_out[:, irr_days] = irr_true[:, :n_irr] # SW out
    
    irr_gw_out = np.zeros((nfield_crop, gen.nper))
    irr_gw_out[:, irr_days] = irr_true[:, n_irr:] # GW out


# %%
    
    # need separte hdf5 for each year because total is 300MB, group by crop in array
    fn = join(base_model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
    crop_arr_to_h5(pc_all, crop, fn)
    
    # applied water (GW and SW are separate)
    fn = join(base_model_ws, 'field_SWB', "GW_applied_water_WY"+str(year)+".hdf5")
    crop_arr_to_h5(irr_gw_out, crop, fn)
    
    fn = join(base_model_ws, 'field_SWB', "SW_applied_water_WY"+str(year)+".hdf5")
    crop_arr_to_h5(irr_sw_out, crop, fn)
