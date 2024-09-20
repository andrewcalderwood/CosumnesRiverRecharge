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
# Version to use for running the connected model code with modflow split up into years as needed
#
# - This script copy and pastes input files that do not change into each year's run folder
# - The script updates inputs that have minor changes such as SFR/GHB/CHD
# - The script pre-processes additional inputs for WEL/RCh which will be read by the connection script when determining the recharge and pumping from each year

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename
import sys
import glob
from importlib import reload
import shutil

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
run_dir = 'C:/WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
# run_dir = 'D://WRDAPP/GWFlowModel'

# loadpth = run_dir +'/Cosumnes/levee_setback/streamflow/'
# # model_nam = 'setback_streamflow'
# model_nam = 'historical_streamflow'

loadpth = run_dir +'/Cosumnes/Regional/'

model_nam = 'historical_simple_geology_reconnection'
model_nam = 'input_write_2014_2020'

base_model_ws = loadpth+model_nam


# %%
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt')

nrow,ncol,nlay,delr,delc = (m.dis.nrow, m.dis.ncol, m.dis.nlay, m.dis.delr, m.dis.delc)

if 'LPF' in m.get_package_list():
    gel_nam = 'LPF'
else:
    gel_nam = 'UPW'
gel = m.__getattr__(gel_nam)

# %%
loadpth = run_dir +'/Cosumnes/Economic/'

m.model_ws = join(loadpth, model_nam, 'crop_modflow')
# # drop HOB since we don't want to update it
# m.remove_package('HOB')
# # re-write name file before copying
# # causes issue with file path for tab file
# m.write_name_file()

# %%


# doesn't change between realizations
gel.write_file()
# test to see if model will run with longer itemp, owhm might auto correct
m.chd.write_file()
m.ghb.write_file()
m.nwt.write_file()

m.sfr.write_file()
m.lak.write_file()

# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')


# %%
dem_data = np.copy(m.dis.top.array)
botm = np.copy(m.dis.botm.array)

# %%
# For the tab files the left column is time (in model units) and the right column is flow (model units)
# Time is days, flow is cubic meters per day
# USGS presents flow in cfs (cubic feet per second)
inflow = pd.read_csv(join(gwfm_dir, 'SFR_data', 'MB_daily_flow_cfs.csv'), index_col = 'datetime', parse_dates = True)
# covnert flow from cubic feet per second to cubic meters per day
inflow['flow_cmd'] = inflow.flow_cfs * (86400/(3.28**3))

# save tabfiles dict to variable for reference
tabfiles_dict = m.sfr.tabfiles_dict

# %%
# deer creek doesn't flow in dry-season
# have minimum flow of 100/200 cfs to start deer creek flowing
# deer creek is approximated as about 10% of Cosumnes flow during wet season
dc_min_flow = 100*(86400/(3.28**3))
inflow_dc = inflow.copy()
# anything above the flow threshold is 10% of Cosumnes
inflow_dc.flow_cmd *= 0.1
# anything below the flow threshold is 0
inflow_dc.loc[inflow_dc.flow_cmd<dc_min_flow,'flow_cmd'] = 0

# %%
import h5py
uzf_dir = join(gwfm_dir, 'UZF_data')
nrow_p, ncol_p = 100,230


# %%
wel_dir = join(gwfm_dir, 'WEL_data')
uzf_dir = join(gwfm_dir, 'UZF_data')

# %% [markdown]
# Write static modflow files into the main directory including LPF, GHB, CHD. LPF (21 MB) with LAK (10 MB) and SFR (40 KB) will not need to be written as there is no dependence on stress periods. GHB, CHD, LAK, and SFR (20 MB) will need to be overwritten or saved multiple times as they have a change due to stress periods with ITMP. Pre-processing and writing output for each of these will save runtime later, but take up about 1.5 GB of storage.
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
# ## specify dates where modflow will start 
# all_run_dates = pd.DataFrame()
# # yn = 0
# # y = years[yn]
# for y in years:
#     run_dates = ymd2dt(y, season.month_run, season.day_run)
#     run_dates = run_dates.drop_duplicates().sort_values()
#     run_dates = pd.DataFrame(run_dates).assign(use='irrigation')
#     crop_date = ymd2dt(y, month_crop, day_crop)
#     crop_date = pd.DataFrame(crop_date).assign(use='crop')
#     all_run_dates = pd.concat((all_run_dates, crop_date, run_dates))

# all_run_dates
# simple code to set dates for april 1
all_run_dates = pd.date_range(all_strt_date, all_end_date,freq='AS-Apr')
all_run_dates = pd.DataFrame(all_run_dates).assign(use='irrigation')
# add total start and end dates
all_run_dates = pd.concat((pd.DataFrame([all_strt_date]).assign(use='start'), all_run_dates))
all_run_dates = pd.concat((pd.DataFrame([all_end_date]).assign(use='end'), all_run_dates))
all_run_dates=all_run_dates.sort_values(0).reset_index(drop=True).rename(columns={0:'date'})

# %%
all_run_dates.to_csv(join(m.model_ws, 'all_run_dates.csv'), index=False)
# all_run_dates

# %% [markdown]
# Write out modflow files that are not impacted by irrigation optimization (GHB, CHD, UPW, OC, NWT, DIS)
# - I may want to see about turning this into a function.
# - We can remove HOB since we don't need to update it and it causes issues

# %% [markdown]
# Write out files that receive input from the irrigation input.
#

# %%
# lak_month = flopy.modflow.ModflowLak.load(loadpth + 'crop_modflow/MF.lak', 
#                                           model=m_month, 
                                          # ext_unit_dict={1: {'inuit': 55}}
                                          # tabdata= True, tab_files='MF.bath', tab_units=[57],
                                         # )


# %%
##############################################################################################
## write out the irrigation independent inputs (GHB, CHD, UPW, OC, NWT, DIS)
print('Writing static input files')
# for m_per in [0]:
for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date
    m_end = all_run_dates.iloc[m_per+1].date

    dates = pd.date_range(m_strt, m_end)[:-1]

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
    # all_dates[(all_dates>=m_strt)&(all_dates<m_end)]
    spd = np.where((all_dates>=m_strt)&(all_dates<m_end))[0]

    ##############################################################################################
    model_ws = join(m.model_ws,str(m_strt.date()))

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
    # this code would break if the ghb/chd weren't explicit
    # for a stress period
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
    
    # load SFR file
    sfr_month = flopy.modflow.ModflowSfr2.load(join(m.model_ws,'MF.sfr'), model=m_month)
    # need to update the options for SFR to represent the number of tab file rows
    # Add option block at the top of the sfr input file for tabfiles
    for k in tabfiles_dict.keys():
        tabfiles_dict[k]['numval'] = nper
    sfr_month.tabfiles_dict = tabfiles_dict
    num_tab = str(len(tabfiles_dict.keys()))
    options_line = ' reachinput transroute tabfiles '+num_tab+' ' + str(nper) + ' no_reach_layer_change'
    tab_option = flopy.utils.OptionBlock(options_line = options_line, package = sfr_month, block = True)
    sfr_month.options = tab_option
    # since it doesn't seem to work to copy sfr then we may have to just 
    # redefine the input for all pieces
    # sfr_month = flopy.modflow.ModflowSfr2(model = m_month, nstrm = sfr.nstrm, nss = sfr.nss, nparseg = sfr.nparseg, 
    #                                 const = sfr.const, dleak = sfr.dleak, ipakcb = sfr.ipakcb, istcb2 = sfr.istcb2, 
    #                           isfropt = sfr.isfropt, nstrail = sfr.nstrail, isuzn = sfr.isuzn, irtflg = sfr.irtflg, 
    #                           numtim = sfr.numtim, weight = sfr.weight, flwtol = sfr.flwtol,
    #                                 reachinput=True, transroute=True, tabfiles=True,
    #                                 tabfiles_dict=sfr.tabfiles_dict)
    
    dis.write_file()
    ghb_month.write_file()
    chd_month.write_file()
    sfr_month.write_file()
    oc.write_file()
    
    print('Done with:', str(m_strt.date()))

# %%
sfr_month.tabfiles_dict

# %%
for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date
    m_end = all_run_dates.iloc[m_per+1].date

    # update SFR tab file for Michigan Bar
    # filter out data between the stress period dates
    mon_inflow = inflow.loc[m_strt:m_end].copy()
    mon_inflow['time'] = (mon_inflow.index-m_strt).days
    # convert to array for writing output
    time_flow = np.transpose((mon_inflow.time, mon_inflow.flow_cmd))
    # save file to each model ws for reference
    np.savetxt(join(model_ws,'MF.tab'), time_flow, delimiter = '\t')
    
    # update SFR tab file for Michigan Bar
    mon_inflow = inflow_dc.loc[m_strt:m_end].copy()
    mon_inflow['time'] = (mon_inflow.index-m_strt).days
    # convert to array for writing output
    time_flow = np.transpose((mon_inflow.time, mon_inflow.flow_cmd))
    # save file to each model ws for reference
    np.savetxt(join(model_ws,'MF.dc.tab'), time_flow, delimiter = '\t')

# %%
m_strt, m_end, nper, spd.shape
m_strt+pd.DateOffset(days=nper)


# %%
# for files that don't need updates
# they can be diretcly copied to save file formatting write time
# # copy over basic package to enable flopy to read in the model, will update start heads later

files_copy = ['MF.lak','MF.bath','MF.evt', 'MF.gage', 'MF.upw','MF.nwt', 'MF.bas']

for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date
    model_ws = join(m.model_ws, str(m_strt.date()))
    for f in files_copy:
        shutil.copy(join(base_model_ws,f), join(model_ws,f))

# %%
model_ws

# %%
# name file had HOB manually removed
files_copy = ['MF.nam']

for m_per in np.arange(0, all_run_dates.shape[0]-1):
    m_strt = all_run_dates.iloc[m_per].date
    model_ws = join(m.model_ws, str(m_strt.date()))
    for f in files_copy:
        shutil.copy(join(m.model_ws,f), join(model_ws,f))

# %% [markdown]
# This code should also run the first period from October to April so that the next script can reference the model output from March.
# - we can't just re-use the standard gridded output I made becaues this assumes some irrigation for farm land. The soil budget must be done separate for all fields.
#

# %%
# for now the lake package only needs flux_data updated
# future versions might use a tab file for flux_data intead
# actually the lake file simply re-uses itself
# m.lak.flux_data[0]

# %%
m_per = 0
m_strt = all_run_dates.iloc[m_per].date
m_end = all_run_dates.iloc[m_per+1].date
spd = np.where((all_dates>=m_strt)&(all_dates<m_end))[0]

model_ws = join(m.model_ws, str(m_strt.date()))

load_only = ['DIS']
# switch to modflow nwt to enable option bloack for use in owhm
m_month = flopy.modflow.Modflow.load('MF.nam',  model_ws= model_ws,
                                    load_only = load_only)

# overwrite files that change
wel_spd = dict()
rech_spd = dict()
# the first model period is October to April where we will assume pumping is 0
# for Ag
for n, t in enumerate(spd):
    wel_spd[n] = np.array([[0,0,0,0]]) # m.wel.stress_period_data[t]
    # wel_spd[n] = np.array([[int(nlay/2), int(nrow/2),int(ncol/2),0]] )
    rech_spd[n] = m.rch.rech.array[t,0]

wel_month = flopy.modflow.ModflowWel(model=m_month, stress_period_data = wel_spd, ipakcb=55)
rch_month = flopy.modflow.ModflowRch(model=m_month, nrchop = 3, rech =  rech_spd, ipakcb=55)

wel_month.write_file()
rch_month.write_file()


success, buff = m_month.run_model()

# %%

# %% [markdown]
# The previous version that simply rewrote the well and recharge from old inputs will change to use the new soil water budget output for all fields, I may want to write another intermediate script to transfer the soil water budget from by crop to one dataset.

# %%
