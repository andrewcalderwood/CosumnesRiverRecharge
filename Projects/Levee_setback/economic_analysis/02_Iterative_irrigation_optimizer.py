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
# from importlib import reload
# reload(swb)
# from Basic_soil_budget_monthly import calc_yield, calc_profit, prep_soil
#, calc_S, calc_pc


# %%
    
import h5py
# using write deletes other datasets, may need to use append instead
def crop_arr_to_h5(arr, crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "a") as f:
        grp = f.require_group('array') # makes sure group exists
        grp.attrs['units'] = 'meters/day'
        grp.attrs['description'] = 'Rows represent the soil units and columns represent the days in the season'
        dset = grp.require_dataset(crop, arr.shape, dtype='f', compression="gzip", compression_opts=4)
        dset[:] = arr
    
def crop_yr_to_h5(arr, crop, year, dates, soil_crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "a") as f:
        grp = f.require_group(crop) # makes sure group exists
        grp.attrs['units'] = 'meters/day'
        grp.attrs['description'] = 'Rows represent the soil units and columns represent the days in the season'
        dset = grp.require_dataset(year, arr.shape, dtype='f', compression="gzip", compression_opts=4)
        dset[:] = arr
        dset = grp.require_dataset(year, arr.shape, dtype='f', compression="gzip", compression_opts=4)
        dset[:] = arr


# %%
# resampled ground surface elevation
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Regional/'
model_ws = loadpth+'historical_simple_geology_reconnection'
load_only = ['DIS','BAS6']
m = flopy.modflow.Modflow.load('MF.nam', model_ws=model_ws, 
                                exe_name='mf-owhm.exe', version='mfnwt', load_only=load_only)
botm = np.copy(m.dis.botm.array)
os.makedirs(join(model_ws, 'field_SWB'), exist_ok=True)

# %% [markdown]
# A simple way to operate the script is to have it read in a file with the year that is specified by the main script.

# %%
model_ws

# %%
year = int(2015)
# crop='Grape'
crop='Corn'
# crop='Alfalfa'
# crop='Pasture' # will require extra work due to AUM vs hay
# crop = 'Misc Grain and Hay'
crop_list = ['Corn','Alfalfa','Pasture','Misc Grain and Hay', 'Grape']



# %%

def init_h5(h5_fn):
    """ Initiate hdf5 files for the given year before appending data for each crop"""
    with h5py.File(h5_fn, "w") as f:
        grp = f.require_group('array') # makes sure group exists
        grp.attrs['units'] = 'meters/day'
        grp.attrs['description'] = 'Rows represent the soil units and columns represent the days in the season'

for var in ['percolation','GW_applied_water', 'SW_applied_water']:
    name = join(model_ws, 'field_SWB', var + '_WY'+str(year)+'.hdf5')
    init_h5(name)

# %%
for crop in crop_list:
    print(crop)

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


    # %%
    # the model will run the irrigation optimizer on specified dates (multiple crops can be done at once or in sequence)
    # the modflow model will be run for the periods between the specified irrigation optimizer dates
    loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Regional/'
    
    # model_ws = loadpth + 'crop_soilbudget/'+crop+'_'+str(strt_date.date())


    # %%
    # fn = join(data_dir,'static_model_inputs.xlsx')
    
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
    ## Load crop type predictions
    
    # double check with Yusuke whether crop2016 or crop_cat2016 is the original or predicted
    crop_wide = pd.read_excel(join(proj_dir, 'crop_prediction', 'wide_dataset.xls'))
    # if a dataset has a diversion then enable decision between SW and GW
    # very few parcels are considered to have a POD (0.5%)
    frac_no_pod = (crop_wide.pod=='No Point of Diversion on Parcel').sum()/len(crop_wide)
    print('%.2f parcels have no POD' %(frac_no_pod*100))
    # crop_wide.waterbody_type.unique()
    # filter to relevant columns for SWB
    keep_cols = crop_wide.columns.str.extract(r'(crop\d{4})').dropna().iloc[:,0].tolist()
    crop_wide = crop_wide.set_index(['parcel_id','pod'])[keep_cols]
    # crop_wide = crop_wide.set_index(['parcel_id','pod'])[['crop2016']]

    # %%
    # long format then filter for year
    crop_long = crop_wide.melt(ignore_index=False,value_name='name', var_name='crop_yr')
    crop_long['year'] = crop_long.crop_yr.str.extract(r'(\d+)').astype(int)


    # %%
    # create dictionary of crop variables and join multiples into list
    var_crops_dict = {}
    for v in var_crops.index.unique():
        var_crops_dict[v] = np.asarray(var_crops[v].tolist())
    
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
    # Misc. grain and hay causes NAs


# %%
    
    gap_irr = var_crops['gap_irr'] # Number of days between irrigations
    n_irr = np.floor(len(dates)/gap_irr).astype(int) + 1 # Calculate number of irrigations
    irr_days = np.arange(0, (n_irr*gap_irr-1), gap_irr).astype(int) # Calculate days on which irrigation takes place
    # irr_days


    # %%
    gen_dict['K_Y'] = var_yield.set_index('dt')['K_Y'].reindex(dates).ffill().values
    gen_dict['n_irr'] = n_irr
    gen_dict['irr_days'] = np.arange(0, (n_irr*gap_irr-1), gap_irr).astype(int) # Calculate days on which irrigation takes place

# %% [markdown]
#     # ## Test SWB on one field

    # %%
    # decide later best way to identify but in future the predicted data will be available for each year
    crop_in = crop_long[crop_long.year==2016].reset_index()

    # %%
    avg_irr_eff = pd.read_csv(join(proj_dir, 'model_inputs', 'avg_irr_eff_by_crop.csv'),index_col=0)
    irr_eff_mult = 100/avg_irr_eff.loc[crop_dict[crop]].Avg_eff
    gen_dict['irr_eff_mult'] = 100/avg_irr_eff.loc[crop_dict[crop]].Avg_eff
    print(pred_dict[crop], crop)
    soil_crop = swb.load_soil(pred_dict[crop], crop_in)
    # another option instead of looping over texture classes would be to sort by texture
    # then the data would be ordered to reference the nearest similar texture
    
    nfield_crop = len(soil_crop)


    # %%
    soil_path = join(uzf_dir,'clean_soil_data')
    # connection of ag fields to all grid cells
    grid_soil = pd.read_csv(join(soil_path, 'parcel_field_to_cell.csv'),index_col=0)
    # subset to fields for crop
    grid_crop = grid_soil.merge(soil_crop[['UniqueID']])


    # %%
    # also need shapefile of pumping well locations for each parcel
    parcel_wells = gpd.read_file(join(gwfm_dir, 'WEL_data', 'parcels_to_wells', 'parcels_to_wells.shp'))
    frow = parcel_wells.row-1
    fcol = parcel_wells.column-1
    # # parcel_wells layers (make 1-based
    parcel_wells['layer'] = get_layer_from_elev(dem_data[frow,fcol] - parcel_wells.depth_m*0.9, botm[:, frow,fcol], m.dis.nlay) + 1
    # subset to crops
    parcel_wells.UniqueID = parcel_wells.UniqueID.astype(int)
    crop_wells = soil_crop[['UniqueID']].merge(parcel_wells)
    crop_dem = dem_data[crop_wells.row-1, crop_wells.column-1]

    # %%
    # sample heads
    # model_ws_last = loadpth + 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date())
    # hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')
    hdobj = flopy.utils.HeadFile(model_ws + '/MF.hds')
    # sp_last = hdobj.get_kstpkper()[-1]


    # %%
    m_strt = pd.to_datetime(str(year-1)+'-10-1')


    # %%
    # it's going to be very slow to sample the time series for all wells for all stress periods and time steps
    # it does seem more efficient to sample on a monthly scale like a farmer would (irrigation dates) and interpolate
    well_kij = list(zip(parcel_wells.layer-1, parcel_wells.row-1, parcel_wells.column-1))
    # hdobj.get_ts(well_kij)
    head = hdobj.get_data(kstpkper = (0,0))
    head[crop_wells.layer-1, crop_wells.row-1, crop_wells.column-1].shape

    # %%
    # for spd in np.arange(1,
    crop_dtw = np.zeros((nper,nfield_crop))
    spd_strt = (strt_date-m_strt).days
    dtw_days = np.hstack((irr_days,[nper-1])) # irr_days alternate
    for n, kper in enumerate(dtw_days):
    # for kper in np.arange(0,nper): # daily sampling
        head = hdobj.get_data(kstpkper = (0,kper+spd_strt))
        crop_dtw[kper,:] = crop_dem - head[crop_wells.layer-1, crop_wells.row-1, crop_wells.column-1]
        # if n >=1: # not necessary except for plotting
        #     # fill in between as farmer might expect linear change
        #     crop_dtw[irr_days[n-1]:kper+1,:] = np.linspace(crop_dtw[irr_days[n-1]], crop_dtw[kper], kper-irr_days[n-1]+1)
    # convert from meters to feet
    crop_dtw /= 0.3048

    # %%
    import f_gw_dtw_extract
    reload(f_gw_dtw_extract)
    from f_gw_dtw_extract import get_dtw


    # %%
    dtw = np.loadtxt(join(model_ws, 'field_SWB', 'dtw_ft_parcels_'+str(year)+'.txt'))
    dtw_df = pd.read_csv(join(model_ws, 'field_SWB', 'dtw_ft_parcels_'+str(year)+'.csv'), index_col=0)


# %% [markdown]
#     # ## Iterate over each unique soil condition

    # %%
    # import swb_functions
    # reload(swb_functions)
    from swb_functions import prep_soil_dict, calc_S, calc_pc, calc_yield, calc_profit, run_swb, mak_irr_con


    # %%
    # so if you have a dictionary d and want to access (read) its values with the syntax x.foo instead of the clumsier d['foo'], just do
    # convert a dictionary to an object with object style referencing
    class cost_variables(object):
      def __init__(self, adict):
        self.__dict__.update(adict)
    
    # convert dictionary of variables to class for easier referencing, constant over different soil
    gen = cost_variables(gen_dict)
    
    from scipy.optimize import Bounds
    bounds = Bounds(lb = 0)


    # %%
    from scipy.optimize import minimize

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
            # now I know some crops like Misc grain and hay have losses then this doesn't work
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
    print('Total time was %.2f min' %((t1-t0)/60),'for', crop)

# %%
# troubleshooting code
# pi, pc, K_S = run_swb(irr_lvl, soil, gen, rain, ETc, dtw_arr,arrays=True)
# Y_A = calc_yield(ETc, K_S, gen)
# calc_profit(Y_A, dtw_arr, irr_gw, irr_sw, gen)

    # %%
    # import matplotlib.pyplot as plt
    # plt.bar(x = np.arange(0, irr_all.shape[1]), height=irr_all[0])

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
    
    p_true = np.copy(p_all)
    pc_all = np.zeros((nfield_crop, gen.nper))
    # for ns in np.arange(0,nfield_crop):
    for ns in np.arange(0,10):
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
# def arr_2_df(arr,name):
#     """
#     Clean array with rows of fields and columns of dates into a long form dataframe
#     """
#     df = pd.DataFrame(arr, columns=dates, index=soil_crop.UniqueID)
#     df_long = df.melt(ignore_index=False, var_name='date',value_name=name)
#     return(df_long)

# %%
# pd.merge(pd.merge(arr_2_df(pc_all, 'pc'), arr_2_df(irr_sw_out, 'irr_sw')), arr_2_df(irr_gw_out, 'irr_gw'))
# arr_2_df(pc_all, 'pc').merge(arr_2_df(irr_sw_out, 'irr_sw'), on=['UniqueID','date'])
# arr_2_df(pc_all, 'pc')
# arr_2_df(irr_sw_out, 'irr_sw')

# %%
    
    # need separte hdf5 for each year because total is 300MB, group by crop in array
    fn = join(model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
    crop_arr_to_h5(pc_all, crop, fn)
    
    # applied water (GW and SW are separate)
    fn = join(model_ws, 'field_SWB', "GW_applied_water_WY"+str(year)+".hdf5")
    crop_arr_to_h5(irr_gw_out, crop, fn)
    
    fn = join(model_ws, 'field_SWB', "SW_applied_water_WY"+str(year)+".hdf5")
    crop_arr_to_h5(irr_sw_out, crop, fn)

# %%
