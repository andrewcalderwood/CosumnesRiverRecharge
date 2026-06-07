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
# Goal: Convert the final output into a more useable format for post-processing

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import glob
import sys
import re

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
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
# for batch file case we want to ignore the consant h5py deprecation warning
import warnings

warnings.filterwarnings("ignore")
# warnings.filterwarnings("ignore", category=DeprecationWarning) 

# %%

from functions.output_processing import get_wb_by_parcel
from functions.f_gw_dtw_extract import sample_dtw, avg_heads
import functions.Basic_soil_budget_monthly as swb

# initialize HDF5 files for the year
from functions.data_functions import init_h5


# %%
import f_rep_swb_profit_opt
reload(f_rep_swb_profit_opt)
from f_rep_swb_profit_opt import load_run_swb


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 


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
def summarize_output_year(loadpth, m_nam, m_per, parcels, input_name):
    '''
    INPUT:
    m_nam: model name
    m_per: the integer index of the current model period (1 = year 1)
    parcels: geopandas dataframe with UniqueID geometry to allow area calculation
    '''
    parcels['acres'] = parcels.area_m2/(43560*0.3048**2)
    parcels.UniqueID = parcels.UniqueID.astype(int)

# %%

    model_ws = join(loadpth, m_nam)
    out_dir = join(model_ws, 'output_clean')

    # simpler way to get base model workspace is remove R\d{1,3} since all should follow this format
    # made it {1,3} digits to allow for hundreds to represent larger shifts in base conditions (e.g., 200 series have different costs)
    m_nam_base = re.sub(r'_R\d{1,3}', '', m_nam)


    # %%
    # provide representative soil water budget folder
    swb_ws = join(model_ws, 'rep_crop_soilbudget')
    os.makedirs(join(swb_ws, 'output'), exist_ok=True)

    # %%
    # define modflow model WS to reference for modflow input
    m_model_ws = join(dirname(loadpth), 'Regional', m_nam_base)


    # %%
    m_dim = np.loadtxt(join(m_model_ws, 'model_metadata.txt')).astype(int)


    # %%
    # now have this file pre-created by model_connect since it already estimates layer
    parcel_wells = pd.read_csv(join(model_ws, 'crop_modflow', 'parcel_wells_with_layer.csv'))
    parcel_wells.UniqueID = parcel_wells.UniqueID.astype(int)


    # %%
    all_run_dates = pd.read_csv(join(model_ws, 'crop_modflow', 'all_run_dates.csv'), parse_dates=['date'])


    # %%
    # for m_per in np.arange(1, all_run_dates.shape[0]-1):
    # for m_per in [6]:
    m_strt = all_run_dates.iloc[m_per].date
    year = m_strt.year
    print(year)


# %%

    # initialize SWB folder
    os.makedirs(join(model_ws, 'crop_soilbudget', 'field_SWB'), exist_ok=True)
    for var in ['profit', 'yield', 'percolation','GW_applied_water', 'SW_applied_water']:
        name = join(model_ws, 'crop_soilbudget', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        init_h5(name)

    for var in ['swb_output']:
        name = join(model_ws, 'crop_soilbudget', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        init_h5(name, groups=['wc','ETa', 'rp'])

    # %%
    model_ws_last = join(model_ws, 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date()))
    hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')

    # model output dates
    times = hdobj.get_times()
    difftime = np.diff(times)
    m_dates = all_run_dates.loc[m_per-1].date+np.array(times - np.append(difftime,[0])).astype('timedelta64[D]')
    m_dates = pd.DataFrame(m_dates, columns=['dates']).set_index('dates')
    m_dates['kstpkper'] = hdobj.get_kstpkper()
    # determine dates for spring sampling
    spring_dates = m_dates[m_dates.index.month==3]
    # get head value from last 30 days to avoid using extreme single day value
    spring_heads = avg_heads(spring_dates.kstpkper.values, hdobj, m_dim)
    
    # the dtw conversion runs a little slow
    # get the DTW for the wels in the simulation from the last period
    well_dtw = sample_dtw(spring_heads, parcel_wells)
    # need to make integer for join with crop choice
    well_dtw.UniqueID = well_dtw.UniqueID.astype(int)

    # %%
    # get head value from last 30 days to avoid using extreme single day value
    # now the csv and code to load produce the same results it seems
    # well_dtw = pd.read_csv(join(swb_ws,'field_SWB', 'modflow_spring_dtw_ft_WY'+str(year)+'.csv'), index_col=0)
    # need to make integer for join with crop choice
    # well_dtw.UniqueID = well_dtw.UniqueID.astype(int)


    # %%
    crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    # remove unclassified fallow for identifying fields with irrigation
    crop_in_irr = crop_in[crop_in.name!='Unclassified fallow']
    # now after the fact it would be useful to have summary info on dtw for fallow
    # going to see if anything in the code would break by updating this
    # crop_in_irr = crop_in.copy()
    dtw_simple_df = pd.read_csv(join(swb_ws,'field_SWB', 'dtw_ft_WY'+str(year)+'.csv'),index_col='date')

    # subset parcels with wells to those identified in crop choice
    year_wells = parcel_wells[parcel_wells.UniqueID.isin(crop_in_irr.parcel_id)]

    # %%
    fn = join(swb_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
    print('Crops with SWB results')
    with h5py.File(fn) as dset:
        finished_crops = list(dset['array'].keys())
        print(finished_crops)

    # %%
    # load the processed dataframe with all data
    # this looks back at model_ws, 'field_SWB', where model_ws is swb_ws which is join(model_ws, 'rep_crop_soilbudget')
    # this is just repeating what is done in the main code then, better to read in csv which have updated values for everything
    # pc_df_all, irr_gw_df_all, irr_sw_df_all = get_wb_by_parcel(swb_ws, year, 
    #                  crop_in, finished_crops, dtw_simple_df, well_dtw)

    # %%
    # check to see if the get_wb_by_parcel needs to be run again and why it doesn't just reference the csv files
    # i think only the csv files are updated with establishment irrigation so need to look for that
    # it's not clear why I had these commented out
    irr_gw_df_all = pd.read_csv(join(swb_ws, 'output', 'irr_gw_all'+str(year)+'.csv'),index_col=0, parse_dates=['date'])
    irr_sw_df_all = pd.read_csv(join(swb_ws, 'output', 'irr_sw_all'+str(year)+'.csv'),index_col=0, parse_dates=['date'])
    # we don't actually need this percolation data as it is going to be re-calculated in the SWB anyway
    # pc_df_all = pd.read_csv(join(swb_ws, 'output', 'pc_all'+str(year)+'.csv'),index_col=0)

    # info on crop change to inform adjustments to profit and yield
    crop_change_info = pd.read_csv(join(swb_ws, 'output', 'crop_change_info'+str(year)+'.csv'))

# %% [markdown]
# There is an issue with the very first period running 1 day into the next period. After checking the copy_model_modflow it shows that it ends on 3/31/2016 so it doesn't make sense.  
# The model input shows only 182 periods but the get_ts has 183 times, is this because it includes the initial heads?

    # %%
    # well locations to sample for head
    wells_idx = list(zip(year_wells.layer-1, year_wells.row-1, year_wells.column-1))

    # the head data loaded here can be generic across all crops then filtered down
    print('Loading previously simulated heads')
    hd_ts_all = pd.DataFrame()
    # iterate over previous and current year to get complete time series (Nov (-1 year) -Dec)
    for n in [-1,0]:
        model_ws_year = join(model_ws, 'crop_modflow/'+str(all_run_dates.loc[m_per+n].date.date()))
        hdobj_year = flopy.utils.HeadFile(model_ws_year + '/MF.hds')
    
        # model output dates
        times = hdobj_year.get_times()
        difftime = np.diff(np.append([0],times) )
        m_dates_year = all_run_dates.loc[m_per+n].date+np.array(times - difftime).astype('timedelta64[D]')
        m_dates_year = pd.DataFrame(m_dates_year, columns=['dates']).set_index('dates')
        m_dates_year['kstpkper'] = hdobj_year.get_kstpkper()
    
        # get head time series for the wells across the year
        # doesn't take too long
        hd_ts = hdobj_year.get_ts(wells_idx)
        # convert to dataframe for use in load_run_swb
        hd_ts_df = pd.DataFrame(hd_ts)
        # update columns to represent uniqueID
        hd_ts_df.columns = ['date']+year_wells.UniqueID.tolist()
        hd_ts_df['date'] = m_dates_year.index
        hd_ts_all = pd.concat((hd_ts_all, hd_ts_df))


    # %%
    # these all need to be converted to depth to water as they are currently head values (usually in negatives and in meters)
    dtw_df = hd_ts_all.set_index('date').copy()
    # calculate the depth to water as dtw_ft = (DEM (m) - WSE (m))/0.3048
    dtw_df = dtw_df.multiply(-1).add(year_wells.dem.values, axis=1).multiply(1/0.3048)
    # the function (SWB) expects the head values to be provided for all fields then filters by crop_in for the crop

    # %%
    for crop in finished_crops:
    # for crop in ['Alfalfa']:

        print(crop)

        # %%
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop)
        # need to account for when crops aren't predicted and skip them
        # if pred_dict[crop] in pred_crops: 
        print(crop, ':',pred_dict[crop])

        # %%
        yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
        yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
        # get the total extent of the irrigation season (calculation period)
        strt_date = yield_start.min()
        end_date = yield_end.max()
        dates = pd.date_range(strt_date, end_date, freq='D')
        print('Start',strt_date.date(),'end', end_date.date(), 'nper',(end_date-strt_date).days+1)

        # %%
        # irr_all =np.zeros((nfield_crop,2*n_irr))
        # 
        gap_irr = var_crops['gap_irr'] # Number of days between irrigations
        n_irr = np.floor(len(dates)/gap_irr).astype(int) + 1 # Calculate number of irrigations
        irr_days = np.arange(0, (n_irr*gap_irr-1), gap_irr).astype(int) # Calculate days on which irrigation takes place
        # specify dates of irrigation to sample
        irr_dates = pd.Series(strt_date+irr_days.astype('timedelta64[D]'))

        # %%
        print('Preparing irrigation input')
        # get fields for the crop
        crop_df = crop_in[crop_in.name==pred_dict[crop]]
        
        # sample the irrigation rates by crop and irr_date
        irr_gw_crop = irr_gw_df_all.set_index('UniqueID').loc[crop_df.parcel_id].copy()
        irr_gw_crop_dates = irr_gw_crop.reset_index().set_index('date').loc[irr_dates]
        irr_gw_crop_dates = irr_gw_crop_dates.sort_values(['UniqueID','date'])
        irr_gw_crop_dates = irr_gw_crop_dates.reset_index().pivot_table(columns='date', values='rate',index='UniqueID')
        # repeat for Surface water
        irr_sw_crop = irr_sw_df_all.set_index('UniqueID').loc[crop_df.parcel_id].copy()
        irr_sw_crop_dates = irr_sw_crop.reset_index().set_index('date').loc[irr_dates]
        irr_sw_crop_dates = irr_sw_crop_dates.sort_values(['UniqueID','date'])
        # pivot to wide dataframe to match input shape
        irr_sw_crop_dates = irr_sw_crop_dates.reset_index().pivot_table(columns='date', values='rate',index='UniqueID')
    
        # specify irr_all input
        # irr_gw goes first then irr_sw
        irr_all = np.hstack((irr_sw_crop_dates.values,irr_gw_crop_dates.values))

# %%
# the simulation results show that naturally there is little to no irrigation before april so hold off for now
# and I guess set a check to use last years head data or backward fill
    # irr_gw_crop_dates.transpose().plot()

# %% [markdown]
# TODO: alfalfa has two irrigation dates that occur before the new modflow run, in theory we shouldn't need DTW in that time as we were going
# to assume no irrigation or we can sample from the previous year. The actual irrigation estimated is very small in those months so okay for now

        # %%
        # # filter well locations to the crop
        # crop_wells = year_wells[year_wells.UniqueID.isin(crop_df.parcel_id)]
        # save dtw for each crop uniquely for plotting later
        crop_dtw = dtw_df.loc[:,crop_df['parcel_id'].values]
        # dtw_arr = crop_dtw.loc[dates]# save output to reference
        # re-index and forward fill to account for the last year of simulation which will end on 9/30 
        # but crops like alfalfa need data until 10/4
        dtw_arr = crop_dtw.reindex(dates).ffill() # save output to reference
        dtw_arr.to_csv(join(model_ws,'crop_soilbudget','field_dtw', 'dtw_ft_'+crop+'_'+str(year)+'.csv'))

        # %%
        print('Running soil water budget with irrigation and updated DTW to re-calculate yield, profit, and percolation')
        # in theory the best function to use if it works
        load_run_swb(crop, year, crop_in, join(model_ws,'crop_soilbudget'), dtw_df, input_name = input_name, soil_rep = False,
                        run_opt=False, irr_all=irr_all, field_id = 'parcels')
        sys.stdout.flush()

    # %%
    name = join(model_ws, 'crop_soilbudget', 'field_SWB', 'profit' + '_WY'+str(year)+'.hdf5')
    with h5py.File(name) as dset:
        crop_list = list(dset['array'].keys())


    # %%
    print('Updating profit and yield to account for crop change')
    # need to overwrite yield and profit as 0, be careful in this section as it could repeatedly subtract
    fn = join(model_ws, 'crop_soilbudget', 'field_SWB', "profit_WY"+str(year)+".hdf5")
    for crop in crop_list:
    # for crop in ['Grape']:
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year = year, input_name=input_name)
        # identify which parcels had crop changes
        crop_change = crop_in[crop_in.name==pred_dict[crop]].copy()
        crop_change = crop_change.merge(crop_change_info.rename(columns={'UniqueID':'parcel_id'}), how='left')
        # identify integer indexes for referencing with hdf5 files
        change_inds = np.arange(0, len(crop_change))[crop_change.changed==True]
        # only apply cost to change in the first year
        change_cost_inds = np.arange(0, len(crop_change))[crop_change.year_offset==1]    
        if crop=='Grape':
            # only need to update yield for Grape
            with h5py.File(join(model_ws, 'crop_soilbudget', 'field_SWB', "yield_WY"+str(year)+".hdf5"), "r+") as f:
                grp = f['array'] 
                dset = grp[crop]
                # extract yield to subtract revenue below
                change_yield = dset[change_inds][:]
                # where the crop was changed make yield 0
                # in year three it could be up to 5 tons but ignore for simlicity
                dset[change_inds] = 0
        # open file in read/write model
        with h5py.File(fn, "r+") as f:
            grp = f['array'] 
            dset = grp[crop]
            # grape has 0 profit for the first three years
            # to account for irrigation cost still, we could simply subtract the revenue (yield * $/ton)
            if crop=='Grape':
                # where the crop was changed make profit 0
                # dset[change_inds] = 0
                # better to instead remove the revenue from yield (profit is negative so need to add to remove the revenue)
                dset[change_inds] = dset[change_inds][:]  + change_yield*var_crops['p_c']
                # temp = dset[change_inds][:] + change_yield*var_crops['p_c']
            # apply the cost to switch cost to the profit (profit is reverse sign)
            dset[change_cost_inds] += var_crops['p_d']



# %%
# after running the updated swb for a year then it would make sense to load in the hdf5 for profit to calculate the average value but this could also be done
# in a secondary script called by model_connect
# code comes from plot_final_output.py
    print('Loading output to calculate crop average profit and yield')
    # %%
    # load SWB folder
    crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    print('\n', year, end=' - ')
    df_all = pd.DataFrame()

    # for var in ['profit', 'yield', 'percolation','GW_applied_water', 'SW_applied_water']:
    for var in ['profit', 'yield']:
        print(var, end=',')
        name = join(model_ws, 'crop_soilbudget', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        with h5py.File(name) as dset:
            finished_crops = list(dset['array'].keys())
            print(finished_crops, end='.')
        for crop in finished_crops:
            # need dates for time series water budget output
            var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year = year, input_name=input_name)

            # extract output and convert to dataframe with ID columns
            arr = read_crop_arr_h5(crop, name)
            df = pd.DataFrame(arr, columns=['value']).assign(crop=crop, year=year, var=var)
            # add parcel information back
            df = pd.concat((df,crop_in[crop_in.name==pred_dict[crop]].reset_index()),axis=1)
            df_all = pd.concat((df_all, df))

    # correct profit from negative to positive
    df_all.loc[df_all['var']=='profit','value'] *= -1
    # fix name before ID join
    df_all = df_all.rename(columns={'parcel_id':'UniqueID'})
    # rename as econ for plotting reference
    df_econ = df_all.merge(parcels[['UniqueID','acres']])

    # %%
    # scale value rates (1/acre) into totals 
    df_econ['total_value'] = df_econ['value']*df_econ.acres
    # we want to aggregate yield and profit by the profit/acre and yield/acre to the total
    # look at average rate, and summed total (scaled by acreage)
    df_econ_agg = df_econ.groupby(['crop','name','var','year'])[['total_value','value']].agg({'total_value':'sum', 'value':'mean'})
    # df_econ.groupby(['crop','name','var','year'])['total_value'].agg(['sum', 'mean'])

    df_econ_agg = df_econ_agg.reset_index()
    df_econ_agg.year = df_econ_agg.year.astype(str)

    # %%
    print('Saving average profit and yield')
    # save data for read in by model_connect to inform next year's crop choice
    df_econ_agg.to_csv(join(out_dir, 'profit_yield_long_'+str(year)+'.csv'))
