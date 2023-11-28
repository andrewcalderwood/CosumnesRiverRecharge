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

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename
import sys
import glob
from importlib import reload

import pandas as pd
import numpy as np
from scipy.stats import hmean, gmean

# import calendar
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import pyproj # for converting proj4string
import shapely
import geopandas as gpd
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm



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
# run_dir = 'C://WRDAPP/GWFlowModel'
run_dir = 'F://WRDAPP/GWFlowModel'

# loadpth = run_dir +'/Cosumnes/levee_setback/streamflow/'
# # model_nam = 'setback_streamflow'
# model_nam = 'historical_streamflow'

loadpth = run_dir +'/Cosumnes/Regional/'
# model_nam = 'setback_streamflow'
model_nam = 'historical_simple_geology'

base_model_ws = loadpth+model_nam


# %%
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt')

nrow,ncol,nlay,delr,delc = (m.dis.nrow, m.dis.ncol, m.dis.nlay, m.dis.delr, m.dis.delc)
m.model_ws = loadpth + 'crop_modflow'

if 'LPF' in m.get_package_list():
    gel_nam = 'LPF'
else:
    gel_nam = 'UPW'
gel = m.__getattr__(gel_nam)

# doesn't change between realizations
gel.write_file()
# test to see if model will run with longer itemp, owhm might auto correct
m.chd.write_file()
m.ghb.write_file()

m.nwt.write_file()

# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')


# %%
dem_data = np.copy(m.dis.top.array)
botm = np.copy(m.dis.botm.array)

# %%
import h5py
uzf_dir = join(gwfm_dir, 'UZF_data')
nrow_p, ncol_p = 100,230


# %%
wel_dir = join(gwfm_dir, 'WEL_data')
uzf_dir = join(gwfm_dir, 'UZF_data')

# %% [markdown]
# Write static modflow files into the main directory including LPF, GHB, CHD. LPF (21 MB) will not need to be written as there is no dependence on stress periods. GHB, CHD and SFR (20 MB) will need to be overwritten or saved multiple times as they have a change due to stress periods with ITMP. Pre-processing and writing output for each of these will save runtime later, but take up about 1.5 GB of storage.
#
#
# The RCH package is 492 MB and well package is 1.59 GB but these file sizes will be subdivided for each period so won't take up much more storage than before.

# %%
all_strt_date = pd.to_datetime(m.dis.start_datetime)
all_dates = all_strt_date + (m.dis.perlen.array.cumsum()-1).astype('timedelta64[D]')
all_end_date = all_dates[-1]
print(all_strt_date, all_end_date)
months = pd.date_range(all_strt_date, all_end_date, freq='MS')
years = pd.date_range(all_strt_date, all_end_date, freq='YS').year.values


# %%
def ymd2dt(year, month, day):
    date = pd.to_datetime(str(year)+'-'+ month.astype(str)+'-'+day.astype(str))
    return(date)



# %%
# choose crops on first day of year
month_crop = pd.Series(1)
day_crop = pd.Series(1)

# load summary excel sheet on irrigation optimization
# this will specify the date ranges to run and pause
fn = join(data_dir,'static_model_inputs.xlsx')
season = pd.read_excel(fn, sheet_name='Seasons', comment='#')


# %%
## specify dates where modflow will start 
all_run_dates = pd.DataFrame()
# yn = 0
# y = years[yn]
for y in years:
    run_dates = ymd2dt(y, season.month_run, season.day_run)
    run_dates = run_dates.drop_duplicates().sort_values()
    run_dates = pd.DataFrame(run_dates).assign(use='irrigation')
    crop_date = ymd2dt(y, month_crop, day_crop)
    crop_date = pd.DataFrame(crop_date).assign(use='crop')
    all_run_dates = pd.concat((all_run_dates, crop_date, run_dates))
    
all_run_dates = pd.concat((pd.DataFrame([all_strt_date]).assign(use='start'), all_run_dates))
all_run_dates = pd.concat((pd.DataFrame([all_end_date]).assign(use='end'), all_run_dates))
all_run_dates=all_run_dates.sort_values(0).reset_index(drop=True).rename(columns={0:'date'})

# %%
all_run_dates

# %%


# load the existing deep percolation and ETc (i.e., AW) datasets (hdf5)

# Load the optimized irrigation rates and corresponding deep percolation values
# translate from fields to grid cells (may be already done in other script)

# where the irrigation optimizer ran overwrite the default DP and AW data

# write the recharge and WEL packages


# using the existing full model rewrite the GHB/CHD package for the given dates

# rewrite BAS6 with start heads from previous output

# rewrite geology

# rewrite OC, NWT

# %%
## Potential ETo spatial interpolation from CIMIS
fn = glob.glob(join(uzf_dir,'Cosumnes_dailyET_precip*.csv'))
daily_data = pd.DataFrame()
for file in fn:
    new_data = pd.read_csv(file, index_col = ['Date'], parse_dates = True)
    daily_data = pd.concat((daily_data, new_data))
# units of mm
data_in = daily_data[daily_data['Stn Name']=='Fair Oaks']
# clean up data so columns are by location, units of Precip are in mm
rain_in = data_in.pivot_table(index = 'Date', columns = 'Stn Name', values = 'Precip (mm)')
rain_m = rain_in/1000
# clean up data so columns are by location, units of Precip are in mm
ETo_in = data_in.pivot_table(index = 'Date', columns = 'Stn Name', values = 'ETo (mm)')
ETo_m = ETo_in/1000


# %%
def dwr_etc(strt_date, end_date):
    nper_tr = (end_date-strt_date).days+1
    natETc = np.zeros((nper_tr,nrow_p,ncol_p))
    agETc = np.zeros((nper_tr,nrow_p,ncol_p))

    per_n = 0 
    for y in np.arange(strt_date.year, end_date.year+1):
        # set start and end date for range for the year to be iterated over
        yr_strt = pd.to_datetime(str(y)+'-01-01')
        yr_end = pd.to_datetime(str(y)+'-12-31')
        # get the length of the date range needed for that year
        yearlen = len(pd.date_range(yr_strt, yr_end))
        if yr_strt < strt_date:
            yr_strt = strt_date
        if yr_end > end_date:
            yr_end = end_date
        yr_len = len(pd.date_range(yr_strt, yr_end))
        # load hdf5 files
        f_irr = h5py.File(join(uzf_dir, "dwr_ETc/irrigated_"+str(y)+".hdf5"), "r")
        agETc[per_n:per_n+yr_len,:,:] = f_irr['array'][str(y)][:][yr_strt.dayofyear-1:yr_end.dayofyear,:,:]
        f_irr.close()
        f_nat = h5py.File(join(uzf_dir, "dwr_ETc/native_"+str(y)+".hdf5"), "r")
        natETc[per_n:per_n+yr_len,:,:] = f_nat['array'][str(y)][:][yr_strt.dayofyear-1:yr_end.dayofyear,:,:]
        f_nat.close()
        per_n += yr_len
    # make sure the return value is separate from the loop
    return(agETc, natETc)


# %%
def load_perc(strt_date, end_date):
    nper_tr = (end_date-strt_date).days+1
    # years and array index 
    years = pd.date_range(strt_date,end_date,freq='AS-Oct')
    yr_ind = (years-strt_date).days
    perc = np.zeros((nper_tr, nrow_p,ncol_p))
    # need separte hdf5 for each year because total is 300MB
    for n in np.arange(0,len(yr_ind)-1):
    #     arr = pc[yr_ind[n]:yr_ind[n+1]]
        fn = join(uzf_dir, 'basic_soil_budget',"percolation_WY"+str(years[n].year+1)+".hdf5")
        f = h5py.File(fn, "r")
        arr = f['array']['WY'][:]
        perc[yr_ind[n]:yr_ind[n+1]] = arr
    #     arr_to_h5(arr, fn)
        f.close()
    return(perc)


# %%
agETc, natETc = dwr_etc(m_strt, m_end)
# net ETc should be ETc from ag and native plants joined
ETc = agETc + natETc

# %%
ag_well_depth_arr = np.loadtxt(wel_dir+'/ag_well_depth_arr.tsv', delimiter='\t')

# %%
# simplified ag well layer with just one layer per well
ag_row, ag_col = np.where(ET_ag.sum(axis=0)>0)
ag_well_lay = get_layer_from_elev((dem_data-ag_well_depth_arr*0.9)[ag_row, ag_col], 
                                  botm[:, ag_row, ag_col], m.dis.nlay)
ag_well_lay.shape, ag_row.shape
ag_well_lay = pd.DataFrame(np.transpose((ag_row,ag_col, ag_well_lay)), columns=['row','column','layer'])

# %%
# load prepared daily domestic use data
dom_use = pd.read_csv(join(wel_dir, 'domestic_water_use.csv'), index_col=0, parse_dates=True)
dom_use = dom_use[all_strt_date:all_end_date]

# load data of locations of domestic wells
dom_loc = pd.read_csv(join(wel_dir, 'ag_res_parcel_domestic_wells.csv'), index_col=0)
# make row,column 0 based
dom_loc.row = (dom_loc.row-1).astype(int)
dom_loc.column = (dom_loc.column -1).astype(int)
# aggregate to the cell level, summing area will keep water usage scaling correct
dom_loc = dom_loc.groupby(['node','row','column', 'CITY']).sum(numeric_only=True).reset_index()
# get domestic well layers
dom_wel_bot = (dem_data[dom_loc.row, dom_loc.column]- dom_loc.fill_depth_m).values
dom_loc['layer'] = get_layer_from_elev(dom_wel_bot, botm[:,dom_loc.row, dom_loc.column], m.dis.nlay)

# use either the total area or expected fraction of irrigated area
# dom_loc['pump_scale'] = dom_loc.used_area_acres
dom_loc['pump_scale'] = dom_loc.area_acres

# %%
# the coefficient for open water is 1.2 at all times
ET_water = ETo_m[all_strt_date:all_end_date]*1.2

water_surf = gpd.read_file(join(uzf_dir,'county_landuse','ag_lu_locally_defined.shp'))
water_surf = gpd.overlay(water_surf, grid_p)
water_surf['area_m2'] = water_surf.geometry.area
# make row,column 0 based
water_surf.row = (water_surf.row-1).astype(int)
water_surf.column = (water_surf.column -1).astype(int)
# determine layer
water_surf['depth_m'] = ag_well_depth_arr[water_surf.row, water_surf.column]
wel_bot_elev = dem_data[water_surf.row, water_surf.column] - water_surf.depth_m
water_surf['layer'] = get_layer_from_elev(wel_bot_elev, botm[:,water_surf.row, water_surf.column], m.dis.nlay)

# %%
# create empty dictionary to fill with stress period data
wel_ETc_dict = {}
# end date is not included as a stress period, starting at 1st TR spd (2)
for t in np.arange(0,nper):
    wel_i, wel_j = np.where(ET_ag[t, :, :]>0)
    new_xyz = ag_well_lay.loc[list(zip(wel_i,wel_j))] 
#     wel_ETc = -ET_ag[t-1,wel_i,wel_j]*delr*delr
# use new row,cols because there are more layers to use
#     wel_ETc = -ET_ag_layered[t, new_xyz.rowi, new_xyz.colj]*delr*delr
    wel_ETc = -ET_ag[t, new_xyz.rowi, new_xyz.colj]*delr*delr
    # ['layer','row','column', 'flux'] are necessary for WEL package
    spd_ag = np.stack((new_xyz.layer, new_xyz.rowi, new_xyz.colj,wel_ETc),axis=1)
    # correct by dropping any rows or cols without pumping as some may be added
    spd_ag = spd_ag[spd_ag[:,-1]!=0,:]
    spd_all = np.copy(spd_ag)
    wel_ETc_dict[t] = spd_all

# %%
wel_dict = {}
for t in np.arange(0, nper):
    # for each stress period specify the flux of parcels from the expected domestic well flux time series
    dom_loc['flux'] = - dom_use.loc[dates[t-time_tr0],'flux_m3d']*dom_loc.pump_scale
    wells_dom = dom_loc[['layer','row','column','flux']].values
    # for each stress period specify the flux of water surfaces 
    water_surf['flux'] = -ET_water.loc[dates[t-time_tr0],'Fair Oaks']*water_surf.area_m2
    wells_ws = water_surf[['layer','row','column','flux']].values
    spd_noag = np.vstack((wells_dom, wells_ws))
    spd_all = np.vstack((wel_ETc_dict[t],spd_noag)) 
    wel_dict[t] = spd_all

# %% [markdown]
# Write out modflow files that are not impacted by irrigation optimization (GHB, CHD, UPW, OC, NWT, DIS)

# %%
##############################################################################################
## write out the irrigation independent inputs (GHB, CHD, UPW, OC, NWT, DIS)

# m_per = 1
for m_per in np.arange(1,5): # runs first year to next crop choice
# for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date
    m_end = all_run_dates.iloc[m_per+1].date

    dates = pd.date_range(m_strt, m_end)

    # The number of periods is the number of dates 
    nper = len(dates) 
    # Each period has a length of one because the timestep is one day, have the 1st stress period be out of the date range
    # need to have the transient packages start on the second stress period
    perlen = np.ones(nper)
    # Steady or transient periods
    steady = np.zeros(nper).astype('bool').tolist()
    # Reduce the number of timesteps to decrease run time
    nstp = np.ones(nper)
    # Identify periods to pull data from
    all_dates[(all_dates>=m_strt)&(all_dates<m_end)]
    spd = np.where((all_dates>=m_strt)&(all_dates<m_end))[0]

    ##############################################################################################
    model_ws = loadpth + 'crop_modflow/'+str(m_strt.date())

    # switch to modflow nwt to enable option bloack for use in owhm
    m_month = flopy.modflow.Modflow(modelname = 'MF', exe_name = 'mf-owhm.exe', 
                              version = 'mfnwt', model_ws= model_ws)

    #lenuni = 1 is in ft, lenuni = 2 is in meters
    # itmuni is time unit 5 = years, 4=days, 3 =hours, 2=minutes, 1=seconds
    dis = flopy.modflow.ModflowDis(nrow=nrow, ncol=ncol, 
                                   nlay=nlay, delr=delr, delc=delc,
                                   model=m_month, lenuni = 2, itmuni = 4,
    #                                xul = xul, yul = yul,rotation=rotation, proj4_str=proj4_str,
                                  nper = nper, perlen=perlen, nstp=nstp, steady = steady,
                                  start_datetime = m_strt)
    #
    m_month.dis.botm = np.copy(m.dis.botm.array)
    m_month.dis.top = np.copy(m.dis.top.array)

    # overwrite files that change
    ghb_spd = dict()
    chd_spd = dict()
    for n, t in enumerate(spd):
        ghb_spd[n] = m.ghb.stress_period_data[t]
        chd_spd[n] = m.chd.stress_period_data[t]

    ghb_month = flopy.modflow.ModflowGhb(model=m_month, stress_period_data = ghb_spd, ipakcb=55)
    chd_month = flopy.modflow.ModflowChd(model=m_month,stress_period_data =  chd_spd)

    # For later model runs when all the data is needed to be saved
    oc_spd = {}
    oc_spd = { (j,0): ['save head', 'save budget'] for j in np.arange(0,nper,1)}
    oc_spd[0,0] = ['save head', 'save budget','print budget']
    oc = flopy.modflow.ModflowOc(model = m_month, stress_period_data = oc_spd, compact = True)

    # add LPF connection and re-write name file
    upw_month = flopy.modflow.ModflowUpw.load(loadpth + 'crop_modflow/MF.upw', model=m_month)

    # load NWT file
    nwt_month = flopy.modflow.ModflowNwt.load(loadpth + 'crop_modflow/MF.nwt', model=m_month)

    # m_month.write_name_file()
#     m_month.write_input()
    dis.write_input()
    ghb_month.write_input()
    chd_month.write_input()
    upw_month.write_input()
    oc.write_input()
    nwt_month.write_input()

# %% [markdown]
# Write out files that receive input from the irrigation input.

# %%
##############################################################################################

for m_per in np.arange(1,5): # runs first year to next crop choice
# for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date
    m_end = all_run_dates.iloc[m_per+1].date

    dates = pd.date_range(m_strt, m_end)
    # The number of periods is the number of dates 
    nper = len(dates) 
    # Identify periods to pull data from
    all_dates[(all_dates>=m_strt)&(all_dates<m_end)]
#     spd = np.where((all_dates>=m_strt)&(all_dates<m_end))[0]

    ##############################################################################################
    model_ws = loadpth + 'crop_modflow/'+str(m_strt.date())

    # switch to modflow nwt to enable option bloack for use in owhm
    load_only=['UPW','DIS','OC','NWT']
    m_month = flopy.modflow.Modflow.load('MF.nam',  model_ws= model_ws,
                                        load_only=load_only,
                                        )  
    # vka needed for other packages
    vka = np.copy(m.upw.vka.array)
    ## update deep percolation 
    perc = load_perc(m_strt, m_end)
    # percolation can't exceed vertical conductivity (secondary runoff)
    perc = np.where(perc >vka[0,:,:], vka[0,:,:], perc)

    # have transient recharge start after the 1st spd
    rech_spd = {}
    for j in np.arange(0,nper):
        rech_spd[j] = perc[j,:,:] 

    # nrchop = 3, to highest active cell
    rch = flopy.modflow.ModflowRch(model = m_month, nrchop=3, rech = rech_spd, ipakcb=55)
    
    ## update pumping
    agETc, natETc = dwr_etc(m_strt, m_end)
    # net ETc should be ETc from ag and native plants joined
#     ETc = agETc + natETc
    # already filtering by land type above
    ET_ag = np.copy(agETc)

    wel_spd = dict()
    for n, t in enumerate(spd):
        wel_spd[n] = m.wel.stress_period_data[t]
    # Create well flopy object
    wel_month = flopy.modflow.ModflowWel(m_month, stress_period_data=wel_spd,ipakcb=55)

    strt = np.ones((nlay, nrow, ncol), dtype = np.float32)
    # The model should start in hydraulic connection
    if spd[0]==0:
        strt[:,:,:] = m.dis.top[:,:] #maybe the mean of starting heads i causing issues?
    else:
        model_ws_last = loadpth + 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date())
        hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')
        sp_last = hdobj.get_kstpkper()[-1]
        strt[:,:,:] = hdobj.get_data(sp_last)
    ibound = np.ones([nlay, nrow,ncol])

    # for the first period we use dem as starting point
    # if solver criteria are not met, the model will continue if model percent error is less than stoperror
    bas_month = flopy.modflow.ModflowBas(model = m_month, ibound=ibound, strt = strt)

    # m_month.write_name_file()
#     m_month.write_input()
    bas_month.write_input()
    rch_month.write_input()
    wel_month.write_input()

#     success, buff = m_month.run_model()

# %%
