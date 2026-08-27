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

import pandas as pd
import numpy as np

# standard geospatial python utilities
import shapely
import geopandas as gpd

import h5py


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

import functions.Basic_soil_budget_monthly as swb
# reload(swb)


# %%
import functions.swb_functions
reload(functions.swb_functions)
from functions.swb_functions import prep_soil_dict, calc_S, calc_pc
from functions.swb_functions import calc_yield, calc_profit, choose_water_source
from functions.swb_functions import run_swb, mak_irr_con, mak_irr_con_adj

# %%
from functions.data_functions import init_h5, crop_arr_to_h5, read_crop_arr_h5


# %%

# so if you have a dictionary d and want to access (read) its values with the syntax x.foo instead of the clumsier d['foo'], just do
# convert a dictionary to an object with object style referencing
class cost_variables(object):
  def __init__(self, adict):
    self.__dict__.update(adict)


# %%
# resampled ground surface elevation
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
# # # # # # testing
# year = int(2018)
# crop='Grape'
# # # # # crop='Corn'
# # # crop='Alfalfa'
# # # # crop='Pasture' # will require extra work due to AUM vs hay
# # # # crop = 'Misc Grain and Hay'

# %%
# # # # # # # # # testing
# # # # # loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic/'
# loadpth = 'D://WRDAPP/GWFlowModel/Cosumnes/Economic/'
# # loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Economic/'
# # m_nam = 'input_write_2014_2022'
# m_nam = 'input_write_2014_2025_R203'
# base_model_ws = join(loadpth, m_nam )
# swb_ws = join(base_model_ws, 'crop_soilbudget')
# crop_in = pd.read_csv(join(base_model_ws,'rep_crop_soilbudget', 'field_SWB', 'crop_parcels_'+str(year)+'.csv'), index_col=0)
# # add in the crops from previous years
# for year_previous in np.arange(2015, year):
#     crop_in_previous = pd.read_csv(join(base_model_ws,'rep_crop_soilbudget', 'field_SWB', 'crop_parcels_'+str(year_previous)+'.csv'), index_col=0)
#     # to avoid conflict, add previous years as new columns with convention name_year
#     crop_in_previous = crop_in_previous[['parcel_id','name']].rename(columns={'name':'name_'+str(year_previous)})
#     crop_in = crop_in.merge(crop_in_previous)
        

# soil_rep = True # True is for the complex dtw_df case
# run_opt=True
# soil_rep = False # True is for the complex dtw_df case
# run_opt=False
# field_id = 'parcels'

# sw_con=100
# gw_con=100
# input_name = 'static_model_inputs_no_p_o.xlsx'
# # # # load parcel data for soil_rep=False
# # # # soil_rep=False

# if soil_rep:
#     # dtw_df = pd.read_csv(join(base_model_ws, 'field_SWB', 'dtw','dtw_ft_parcels_'+str(year)+'.csv'), 
#     #                      index_col=0, parse_dates=['dt'])
#     dtw_df = pd.read_csv(join(base_model_ws, 'rep_crop_soilbudget','field_SWB', 'dtw_ft_WY'+str(year)+'.csv'),
#                         index_col=0, parse_dates=['date'])
#     dtw_df.columns = dtw_df.columns.astype(int)
# else:
#     dtw_df = pd.read_csv(join(base_model_ws,'crop_soilbudget','field_dtw', 'dtw_ft_'+crop+'_'+str(year)+'.csv'),index_col=0)
#     dtw_df.index = pd.to_datetime(dtw_df.index)
#     dtw_df.columns = dtw_df.columns.astype(int)
#     from functions.summarize_functions import format_irr_all, adj_irr_rates
#     irr_gw_df_all = pd.read_csv(join(base_model_ws,'rep_crop_soilbudget',  'output', 'irr_gw_all'+str(year)+'.csv'),index_col=0,parse_dates=['date'])
#     irr_sw_df_all = pd.read_csv(join(base_model_ws,'rep_crop_soilbudget', 'output', 'irr_sw_all'+str(year)+'.csv'),index_col=0,parse_dates=['date'])

# base_model_ws = swb_ws


# %%
# # simple representative DTW for linear steps 10 ft to 200 ft
# # with a 5 ft decline from June to December based on observed data
# dtw_avg = pd.DataFrame(pd.date_range(str(year-1)+'-11-1', str(year)+'-12-31'), columns=['date'])
# dtw_avg = dtw_avg.assign(decline = float(0)).set_index('date')
# # dates where a decline date is specified
# decline_dates = dtw_avg.index[dtw_avg.index >=str(year)+'-6-1']
# decline_total = 5
# decline = np.cumsum(np.full(len(decline_dates), decline_total/len(decline_dates)))
# dtw_avg.loc[decline_dates, 'decline'] = decline
# dtw_simple = np.repeat(np.reshape(np.arange(10, 200, 10), (1,-1)), len(dtw_avg), axis=0)
# dtw_simple = dtw_simple + np.reshape(dtw_avg.decline, (-1,1))
# dtw_df =  pd.DataFrame(dtw_simple, dtw_avg.index)

# # plt.plot(dtw_simple[:,0])

# %%


def load_run_swb(crop, year, crop_in, base_model_ws, dtw_df,
                 input_name, 
                 soil_rep = False,
                run_opt=True, irr_all=None, field_id = 'parcels', 
                sw_con=125, gw_con=125,
                ):
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
    run_opt: boolean to identify whether to run the optimization for each field
    irr_all: irrigation for each field to be specified if run_opt=False
    input_name: need to pass on name of spreadsheet with static inputs so can adjust revenue, costs, etc.
    '''
    # %%
    # need to specify year to get specify year values or nearest available
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year=year, input_name = input_name)

    # %%
    # over-ride the default sw_con, gw_con with crop specific values and use more stringent of the pair
    # input irrigation constraints to compare against regional constraints
    irr_con_in = var_irr.set_index('irrigation')['depth_total']
    # print(irr_con_in)
    if irr_con_in.loc['GW'] < gw_con:
        gw_con = irr_con_in.loc['GW']
    if irr_con_in.loc['SW'] < sw_con:
        sw_con = irr_con_in.loc['SW']

    # %%
    #  get the dates for each yield cycle
    yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
    yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
    # get the total extent of the irrigation season (calculation period)
    strt_date = yield_start.min()
    end_date = yield_end.max()
    dates = pd.date_range(strt_date, end_date, freq='D')
    nper = (end_date-strt_date).days +1
    print('\n')
    print('########## Load run SWB ##########')
    print('Start', strt_date.date(), ', End', end_date.date(),', No. days', nper)
    print(base_model_ws)

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
    rain, ETo_df = swb.load_hyd(dates)
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
    gen_dict['nper'] = (end_date-strt_date).days + 1
    
    # days to index for calculating yield impacts, add 1 to include end day
    yield_ind = np.append([0], (yield_end-strt_date).dt.days.values +1)
    gen_dict['yield_ind'] = np.append([0], (yield_end-strt_date).dt.days.values + 1)


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
    print('Irr eff mult', irr_eff_mult)
    print('Eff GW con ', gw_con/irr_eff_mult,' and eff SW con', sw_con/irr_eff_mult)

    soil_crop = swb.load_soil(pred_dict[crop], crop_in)
    if soil_rep:
        # get the average field properties for a representative field
        # with one for a point of diversion and one without
        soil_crop = soil_crop.groupby('pod').mean(numeric_only=True).reset_index()
        # repeat soil_crop data to iterate over different DTW profiles
        # soil_crop = pd.concat([soil_crop]*dtw_df.shape[1]) # old version that fails to sort
        # repeat data no POD first as the output_processing scripts expect this
        soil_crop = pd.concat((pd.concat([soil_crop[soil_crop.pod_bool==0]]*dtw_df.shape[1]),
                 pd.concat([soil_crop[soil_crop.pod_bool==1]]*dtw_df.shape[1])))

        
    # another option instead of looping over texture classes would be to sort by texture
    # then the data would be ordered to reference the nearest similar texture
    nfield_crop = len(soil_crop)
    print('Num crops:', nfield_crop)


    # %%
    # crop_wells = soil_crop[['UniqueID']].merge(parcel_wells)
    # select parcels in the simulation
    # for rep soil only need one dtw profile
    if soil_rep:
        crop_dtw = dtw_df.copy()
        # need to extend the dataset based on number of groups
        ntimes = int(nfield_crop/crop_dtw.shape[1])
        crop_dtw = pd.concat([crop_dtw]*ntimes, axis=1)
    else:
        crop_dtw = dtw_df.loc[:,soil_crop['UniqueID'].values]
        # could also import already ffill data created in 03b_summarize_output
        crop_dtw = crop_dtw.reindex(dates).ffill() 


    # select dates being simulated
    crop_dtw = crop_dtw.loc[dates].values
    # print(crop_dtw.shape)


# %% [markdown]
#     # ## Iterate over each unique soil condition

    # %%
    # convert dictionary of variables to class for easier referencing, constant over different soil
    gen = cost_variables(gen_dict)

    # a scalar indicates it is applied to all variables
    # upper bound here based on crop reasonable upper limit
    bounds = Bounds(lb = 0, ub = (var_irr.depth_max.max()/12)*0.3048)


    # %%
    # # create a dataframe of soil data to save
    soil_keys_keep = ['Ks', 'por', 'eps', 'CN', 'psdi', 'm',
                      'wc_f', 'wc_wp', 'depth', 'taw', 'Smax']
    soil_df_out = pd.DataFrame()

    # %%
    etc_arr = np.zeros((nper))
    for n in np.arange(nper):
        etc_arr[n] = ETc[n]

    # at this point the optimization begins
    # the function could be updated to here to skip the optimization
    # when not needed, in that case we need to specify the irrigation
    # to use in the SWB calculation (irr_all substitute)
    # %%
    # run_opt=True
    min_profit=1000

    # %%
    if run_opt:
        t0 = time.time()
        # save irrigation, function, time
        irr_all =np.zeros((nfield_crop,2*n_irr))
        p_all = np.ones(nfield_crop) # make initial profit negative (positive p_all)
        p_all[:] = min_profit + 1
        t_all = np.zeros(nfield_crop)
        
        for ns in np.arange(0,nfield_crop):
            soil_ag = soil_crop.iloc[[ns]] #keep as dataframe for consistency 
            nfield = soil_ag.shape[0]
        
            dtw_arr = crop_dtw[:,ns]
    
            ## add check for cheaper water source
            # then adjust solver to use only one if significantly cheaper
            water_source = choose_water_source(dtw_arr, gen, mix_fraction=1)
            
            # prep_soil(soil_ag, etc_arr, var_crops)
            soil_dict = prep_soil_dict(soil_ag, etc_arr, var_crops)
            soil = cost_variables(soil_dict)
    
            # each field iteration adds to the soil dataframe
            field_soil_df = pd.DataFrame(np.append([soil_dict[k][0] for k in soil_keys_keep], 
                                                soil_ag.UniqueID.iloc[0])).transpose()
            # add field to dataframe of all field
            soil_df_out = pd.concat((soil_df_out, field_soil_df),axis=0)
            
            # if no POD then no SW irrig
            if soil_ag.pod.iloc[0]=='No Point of Diversion on Parcel':
                sw_con = 0
                n_irr_type=1
                water_source='gw'
            if water_source=='gw':
                sw_con = 0
                n_irr_type=1
            elif water_source=='sw':
                gw_con = 0
                n_irr_type=1
            
            irr_lvl = np.zeros(n_irr_type*n_irr); # Initial irrigation values for optimization
            # review shows usually it should be 4 inches for most crops with 8 for alfalfa
            # starting at 4 seems to run even slower for quick test in jupyter lab
            # starting at 2 also has issues
            irr_lvl[:] = (2/12)*0.3048 # irrigate with 2 inches (convert to meters)
            irr_lvl_base = np.ones(n_irr_type*n_irr)*(2/12)*0.3048
            # the issue with reusing previous irrigation is if a bad decision
            # is made by the solver then it persists so better to start from defaults each time
            if ns > 0:
                if water_source=='gw':
                    irr_lvl[:] = irr_all[ns-1,n_irr:]
                elif water_source=='sw':
                    irr_lvl[:] = irr_all[ns-1,:n_irr]
                else:
                    irr_lvl[:] = irr_all[ns-1]

            print('Irr length:', len(irr_lvl))
            # simple linear keeps both SW/GW
            # the constraint should be divided by the irrigation efficiency multiplier
            # to limit based on the actual water that will be applied
            linear_constraint = mak_irr_con_adj(n_irr, gw_con = gw_con/gen.irr_eff_mult, sw_con = sw_con/gen.irr_eff_mult) 
            # for the linear dtw the start tol (0.01) was too coarse
            # with representative case it makes to use 0.001 tolerance for all and not force positive values
            tol = 0.001  
            # continue optimizing until profit is positive (p is negative and means positive profit)
            # this step may no longer be necessary as several crops are
            # expected to have negative profits (200-400 dollar losses), lower limit to avoid >-1000 dollars
            # and now we are removing p_o so most should be positive
            while p_all[ns] > min_profit:
                # should look at adding irrigation effiency into cost calculation but not irrigation or else the optimizer would simply scale
                # means we should also adjust the linear constraint for each efficiency
                out = minimize(run_swb, irr_lvl, args = (soil, gen, rain, ETc, dtw_arr, water_source),
                               method='trust-constr',
                               # method='SLSQP',
                        constraints = [linear_constraint],
                        bounds=bounds,
                         # options={'verbose':1}
                               tol = tol
                        )
                # decrease tolerance and reset starting irrigation to help solving
                if out.fun > min_profit:
                    tol /=10
                    irr_lvl[:] = np.copy(irr_lvl_base)
                if tol < 1E-4:
                    break # if tolerance gets too small then skip
                # make sure irrigation is saved in the right spot
                if water_source=='gw':
                    irr_all[ns, n_irr:] = out.x
                elif water_source=='sw':
                    irr_all[ns, :n_irr] = out.x
                else:
                    irr_all[ns] = out.x
                p_all[ns] = out.fun
                t_all[ns] = out.execution_time
            # print final results
            print('Soil ', str(ns),'%.2f' %(-out.fun),'$ ,in %.2f' %(out.execution_time/60),'min')
            if (water_source =='sw')&(soil_rep):
                print('Water source is SW so DTW no longer has an impact for representative profiles')
                print('Skipping remaining DTW profiles and holding constant SW irrigation')
                irr_all[ns:, :n_irr] = out.x
                break

        t1 = time.time()
        print('Total time was %.2f min' %((t1-t0)/60), 'for', ns+1,'parcels')


# %% [markdown]
#     #     # # Post-processing
#     #     # 1. run optimized swb with irrigation efficiency
#     #     # 2. save applied water for pumping
#     #     # 3. save deep percolation

# %% [markdown]
# it doesn't make sense to insert a new initial water content to the optimized version because we'd have to average across fields. It does make sense to confirm that initial conditions in spring should be field water content after all the rain.  
#
# To provide more robust estimates of actual cost, the final SWB run should use the actual initial water content, this might result in slightly lower irrigation demand and/or slightly higher percolation.

    # %%
    # if the file doesn't exist then assign value of None which function assumes is field water content
    # load water content from winter SWB module 
    wc_fn = join(base_model_ws, 'crop_soilbudget', 'ag_non_irr',
                                  'swb_water_content_'+str(year)+'.csv')
    if exists(wc_fn):
        wc_df = pd.read_csv(wc_fn,index_col=0, parse_dates=True)
        wc_df.columns = wc_df.columns.astype(int)
        # sample the initial water content for the fields corresponding to the crop
        crop_wc = wc_df.loc[:,soil_crop['UniqueID'].values].copy()
        # sample the water content for the date that irrigation run starts
        crop_wc_init = crop_wc.loc[strt_date].copy()
    else:
        # the first year that uses the fall period from MODFLOW needs to assume initial conditions of field water content
        print('Assuming field water content for initial conditions')
        crop_wc_init = None

# %%
# only for script testing 
# irr_all = format_irr_all(crop, year, crop_in, pred_dict, irr_gw_df_all, irr_sw_df_all)


    # %%
    # runs relatively quickly
    if (not run_opt) & (irr_all is None):
        # we can have irr_all be an optional input
        print('irr_all not specified for soil water budget')
        print('Assuming no irrigation for all fields')
        irr_all = np.zeros((nfield_crop,2*n_irr))
    elif (not run_opt) & (irr_all is not None):
        print('Using input irrigation rates for irr_all')
    print('Calculating true SWB + profit with irrigation efficiency')
    # correct for any irrigation below 0 set by optimization
    irr_all[irr_all<0] = 0
    # scale by irrigation efficiency of the crop after optimizing
    irr_true = irr_all * irr_eff_mult # one efficiency for each crop type
    p_true = np.zeros(nfield_crop) 
    Y_A_true = np.copy(p_true)
    cost_true = np.copy(p_true)
    pc_all = np.zeros((nfield_crop, gen.nper))
    wc_all = np.zeros((nfield_crop, gen.nper))
    ETa_all = np.zeros((nfield_crop, gen.nper))
    rp_all = np.zeros((nfield_crop, gen.nper))

    for ns in np.arange(0,nfield_crop):
        # update variables for each crop
        dtw_arr = crop_dtw[:,ns]
        soil_ag = soil_crop.iloc[[ns]] #keep as dataframe for consistency 
        soil_dict = prep_soil_dict(soil_ag, etc_arr, var_crops)
        soil = cost_variables(soil_dict)
        if crop_wc_init is not None:
            # need to sample to the field that is being run
            crop_wc_init_field = crop_wc_init.loc[soil_ag.UniqueID.values].values 
        else:
            crop_wc_init_field = soil.wc_f
        # run final soil water budget and save output as arrays
        # should now add specification of initial water content after end of winter run
        # it would be good to save the other relevant outputs such as soil moisture
        p_true[ns], pc, K_S, Y_A, wc, ETa, rp, cost  = run_swb(irr_true[ns], soil, gen, rain, ETc, dtw_arr, arrays=True,
                                                        init_wc = crop_wc_init_field
                                                        )
        Y_A_true[ns] = Y_A.sum() # yield comes as an array
        cost_true[ns] = cost # yield comes as an array
        # double check this works for the multiple simple dtw version
        pc_all[ns] = pc[:,0] # original shape meant for multiple fields, but only has one since iteration is over fields
        wc_all[ns] = wc[1:,0] # original shape meant for multiple fields, but only has one since iteration is over fields
        ETa_all[ns] = ETa[:,0] # original shape meant for multiple fields, but only has one since iteration is over fields
        rp_all[ns] = rp[:,0] # original shape meant for multiple fields, but only has one since iteration is over fields

    # %%
    # it might also be helpful to simply save the final time step of water content to save time but can 
    # also post-process this later

    # %%
    # break down irrigation into groundwater and surface water time series
    # here we assume the series goes surface water first then groundwater
    irr_sw_out = np.zeros((nfield_crop, gen.nper))
    irr_sw_out[:, irr_days] = irr_true[:, :n_irr] # SW out
    
    irr_gw_out = np.zeros((nfield_crop, gen.nper))
    irr_gw_out[:, irr_days] = irr_true[:, n_irr:] # GW out


    # %%
    print('Saving output to files')

    # need separte hdf5 for each year because total is 300MB, group by crop in array
    # the current version doesn't get files bigger than 1,300 KB so we could combine years
    # biggest file is percolation because daily basis instead of seasonal or irrigation cycle
    fn = join(base_model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
    crop_arr_to_h5(pc_all, crop, fn)

    # save profit and yield values
    # the profit saved here is negative for minimization
    fn = join(base_model_ws, 'field_SWB', "profit_WY"+str(year)+".hdf5")
    crop_arr_to_h5(p_true, crop, fn)

    fn = join(base_model_ws, 'field_SWB', "cost_WY"+str(year)+".hdf5")
    crop_arr_to_h5(cost_true, crop, fn)
    # does it make sense to simply create another output file to cover the revenue or create a sub file for it under profit?
    # if a new file then it needs to also be initiated at the start and post-processed in f_summarize_output.py

    fn = join(base_model_ws, 'field_SWB', "yield_WY"+str(year)+".hdf5")
    crop_arr_to_h5(Y_A_true, crop, fn)

    # applied water (GW and SW are separate)
    fn = join(base_model_ws, 'field_SWB', "GW_applied_water_WY"+str(year)+".hdf5")
    crop_arr_to_h5(irr_gw_out, crop, fn)
    
    fn = join(base_model_ws, 'field_SWB', "SW_applied_water_WY"+str(year)+".hdf5")
    crop_arr_to_h5(irr_sw_out, crop, fn)

    # %%
    # the best way to save water budget output is likely to put them all in the same hdf5 file
    fn = join(base_model_ws, 'field_SWB', "swb_output_WY"+str(year)+".hdf5")
    crop_arr_to_h5(wc_all, crop, fn, units='meters/day', group='wc')
    crop_arr_to_h5(ETa_all, crop, fn, units='meters/day', group='ETa')        
    crop_arr_to_h5(rp_all, crop, fn, units='meters/day', group='rp')
