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
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

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
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic'

# update to different modflow models here, next step is using the 20 year model
# base_model_ws = loadpth + 'crop_soilbudget'
m_nam = 'historical_simple_geology_reconnection'
m_nam = 'input_write_2014_2020'
model_ws = join(loadpth, m_nam)


# %%
# provide representative soil water budget folder
swb_ws = join(model_ws, 'rep_crop_soilbudget')
os.makedirs(join(swb_ws, 'output'), exist_ok=True)

# %%
# define modflow model WS to reference for modflow input
m_model_ws = join(dirname(loadpth), 'Regional', m_nam)

load_only=['DIS', 'BAS6']

m = flopy.modflow.Modflow.load('MF.nam', model_ws= m_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only)

# %%
# bottom array is needed for referencing well layer
botm = m.dis.botm.array

# %%
from mf_utility import get_layer_from_elev
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')
nlay,nrow,ncol = botm.shape
# also need shapefile of pumping well locations for each parcel
parcel_wells = gpd.read_file(join(gwfm_dir, 'WEL_data', 'parcels_to_wells', 'parcels_to_wells.shp'))
frow = parcel_wells.row-1
fcol = parcel_wells.column-1
# # parcel_wells layers (make 1-based
parcel_wells['layer'] = get_layer_from_elev(dem_data[frow,fcol] - parcel_wells.depth_m*0.9, botm[:, frow,fcol], nlay) + 1
# get elevation
parcel_wells['dem'] = dem_data[parcel_wells.row-1, parcel_wells.column-1]
parcel_wells = parcel_wells[['UniqueID','dem','layer', 'row','column']]
parcel_wells.UniqueID = parcel_wells.UniqueID.astype(int)


# %%
all_run_dates = pd.read_csv(join(model_ws, 'crop_modflow', 'all_run_dates.csv'), parse_dates=['date'])


# %% [markdown]
# Iteration starts here by year and crop lower down

# %%

for m_per in np.arange(1, all_run_dates.shape[0]-1):
# for m_per in [1]:
    m_strt = all_run_dates.iloc[m_per].date
    year = m_strt.year
    print(year)


# %%

    # initialize SWB folder
    os.makedirs(join(model_ws, 'crop_soilbudget', 'field_SWB'), exist_ok=True)
    for var in ['profit', 'yield', 'percolation','GW_applied_water', 'SW_applied_water']:
        name = join(model_ws, 'crop_soilbudget', 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        init_h5(name)

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
    spring_heads = avg_heads(spring_dates.kstpkper.values, hdobj, m)
    
    # the dtw conversion runs a little slow
    # get the DTW for the wels in the simulation from the last period
    well_dtw = sample_dtw(spring_heads, botm)
    # need to make integer for join with crop choice
    well_dtw.UniqueID = well_dtw.UniqueID.astype(int)

    # %%
    # well_dtw = pd.read_csv(join(swb_ws,'field_SWB', 'modflow_spring_dtw_ft_WY'+str(year)+'.csv'))


    # %%
    crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'),index_col=0)
    # remove unclassified fallow for identifying fields with irrigation
    crop_in_irr = crop_in[crop_in.name!='Unclassified fallow']
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
    # load the processed dataframe with all datas
    pc_df_all, irr_gw_df_all, irr_sw_df_all = get_wb_by_parcel(swb_ws, year, 
                     crop_in, finished_crops, dtw_simple_df, well_dtw)

    # %%
    # # this output with the parcel data needs to be saved as well
    # pc_df_all = pd.read_csv(join(swb_ws, 'output', 'pc_all'+str(year)+'.csv'))
    # irr_gw_df_all = pd.read_csv(join(swb_ws, 'output', 'irr_gw_all'+str(year)+'.csv'))
    # irr_sw_df_all = pd.read_csv(join(swb_ws, 'output', 'irr_sw_all'+str(year)+'.csv'))

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
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
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
        irr_all = np.hstack((irr_gw_crop_dates.values,irr_sw_crop_dates.values))

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
        dtw_arr = crop_dtw.loc[dates]# save output to reference
        # re-index and forward fill to account for the last year of simulation which will end on 9/30 
        # but crops like alfalfa need data until 10/4
        dtw_arr = crop_dtw.reindex(dates).ffill() # save output to reference
        dtw_arr.to_csv(join(model_ws,'crop_soilbudget','field_dtw', 'dtw_ft_'+crop+'_'+str(year)+'.csv'))

        # %%
        print('Running soil water budget with irrigation and updated DTW to re-calculate yield, profit, and percolation')
        # in theory the best function to use if it works
        load_run_swb(crop, year, crop_in, join(model_ws,'crop_soilbudget'), dtw_df, soil_rep = False,
                        run_opt=False, irr_all=irr_all, field_id = 'parcels')

# %% [markdown]
# # Comparison of output

# %%
# out_ws = join(model_ws,'crop_soilbudget')

# %%
# # # applied water (GW and SW are separate)
# fn = join(out_ws, 'field_SWB', "yield_WY"+str(year)+".hdf5")
# Y_A = read_crop_arr_h5(crop, fn)

# fn = join(out_ws, 'field_SWB', "profit_WY"+str(year)+".hdf5")
# # the saved profit is mutliplied by negative for minimization so need to make into real profit
# pi = -read_crop_arr_h5(crop, fn)
