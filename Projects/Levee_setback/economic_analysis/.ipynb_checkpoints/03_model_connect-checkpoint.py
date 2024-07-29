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
# Script to run the code for the crop choice model, irrigation optimizer and modflow in a cycle for the years of interest.
#
# 1. Run the crop choice model
# 2. Run the soil water budget optimization
# 3. Update WEL/RCH packages then run MF
# 4. Start the next year

# %%
import h5py
import numpy as np
import pandas as pd
import geopandas as gpd
import os
from os.path import join,exists, dirname, basename, expanduser

# %%
from importlib import reload

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')

# %%
import sys
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)


# %%
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# %%
git_dir = join(doc_dir,'GitHub','CosumnesRiverRecharge')
add_path(join(git_dir,'python_utilities'))

# %%
from report_cln import base_round
from mf_utility import get_layer_from_elev

# %%
from parcelchoicemodelupdate.f_predict_landuse import predict_crops

# %%
import functions.Basic_soil_budget_monthly as swb
reload(swb)
# import f_swb_profit_opt
# reload(f_swb_profit_opt)
# from f_swb_profit_opt import load_run_swb

import f_rep_swb_profit_opt
reload(f_rep_swb_profit_opt)
from f_rep_swb_profit_opt import load_run_swb

# %%
import functions.f_gw_dtw_extract
reload(functions.f_gw_dtw_extract)
from functions.f_gw_dtw_extract import sample_dtw, avg_heads
# from functions.f_gw_dtw_extract import get_dtw
from functions.f_gw_dtw_extract import calc_simple_dtw

# %%
from functions import data_functions
reload(data_functions)
from functions.data_functions import read_crop_arr_h5

# %%
import functions.output_processing
reload(functions.output_processing)
from functions.output_processing import get_local_data, out_arr_to_long_df
from functions.output_processing import get_wb_by_parcel

# %% [markdown]
# If we are going to run the crop/swb for a certain period then there should already have been a MODFLOW model run for that same period with estimates of streamflow and precipitation to drive the other boundary conditions.

# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
# load parcel data for reference as needed
parcels = gpd.read_file(join(proj_dir,'Parcels shapefile/parcels.shp'))
parcels['area_m2'] = parcels.geometry.area

# load dataframe that assigns each field to cells for recharge
field_id = 'parcel'
field_df = pd.read_csv(join(uzf_dir, 'clean_soil_data', field_id+'_field_to_cell.csv'),index_col=0)

# load dataframe with location (row,column) of wells
well_loc_df = gpd.read_file(join(gwfm_dir, 'WEL_data','parcels_to_wells','parcels_to_wells.shp'))
well_loc_merge = well_loc_df[['UniqueID','row','column','depth_m']].copy()
well_loc_merge = well_loc_merge.merge(parcels[['UniqueID', 'area_m2']])
well_loc_merge.UniqueID = well_loc_merge.UniqueID.astype(int)

# %%
# should we may the loadpth economic instead of Regional?
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic/'

# base_model_ws = loadpth + 'crop_soilbudget'
m_model_ws = join(dirname(loadpth),'Regional', 'historical_simple_geology_reconnection')

# %%
load_only=['DIS', 'BAS6']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= m_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only)

# %%
# save ibound array for when updating bas6 for each run
ibound = m.bas6.ibound.array
# bottom array is needed for referencing well layer
botm = m.dis.botm.array

# %%
## add domestic pumping
wel_dir = join(gwfm_dir, 'WEL_data')
# load prepared daily domestic use data
dom_use = pd.read_csv(join(wel_dir, 'domestic_water_use.csv'), index_col=0, parse_dates=True)
# load data of locations of domestic wells
dom_loc = pd.read_csv(join(wel_dir, 'ag_res_parcel_domestic_wells.csv'), index_col=0)
# make row,column 0 based
dom_loc[['row','column']] = (dom_loc[['row','column']]-1).astype(int)
# aggregate to the cell level, summing area will keep water usage scaling correct
dom_loc = dom_loc.groupby(['node','row','column', 'CITY']).sum(numeric_only=True).reset_index()
# get domestic well layers
dom_wel_bot = (m.dis.top[dom_loc.row, dom_loc.column]- dom_loc.fill_depth_m).values
dom_loc['layer'] = get_layer_from_elev(dom_wel_bot, botm[:,dom_loc.row, dom_loc.column], m.dis.nlay)

# use either the total area or expected fraction of irrigated area
dom_loc['pump_scale'] = dom_loc.used_area_acres # avg of 25% less than area_acres (sometimes more though)
# alternate to avoid over-estimation just stick with 2AF/year
dom_loc['pump_scale'] = 1 # avg of 25% less than area_acres (sometimes more though)


# %%
# columns are the ntaive land use fields
pc_df = pd.read_csv( join(proj_dir,'native_parcel_zonalstats', 
                   'swb_percolation.csv'),index_col=0)
# long format to join with row,column
pc_df_long = pc_df.melt(ignore_index=False, 
                        var_name='UniqueID',
                       value_name='pc')
pc_df_long.UniqueID=pc_df_long.UniqueID.astype(int)


# %%
soil_path = join(uzf_dir,'clean_soil_data')
field_ids = 'native'
# reference for field to grid cell
grid_soil = pd.read_csv(join(soil_path, field_ids+'_field_to_cell.csv'), index_col=0)

# %%
# join field to row,column
native_pc = pc_df_long.reset_index().merge(grid_soil)
native_pc = native_pc.rename(columns={'index':'date'})
# calculate recharge rate by scaling for fraction of cell covered 
native_pc['rch_rate'] = native_pc.pc*native_pc.field_area/(200*200)
# update index to datetime
native_pc = native_pc.set_index('date')
native_pc.index=pd.to_datetime(native_pc.index)

# %%
# drop dates and cells with zero recharge as the default recharge is 0
native_pc = native_pc[native_pc.rch_rate!=0]

# %%
# quick budget checking
# native_pc['year'] = native_pc.index.year
# yr_sum = native_pc.groupby(['year','row','column'])['rch_rate'].sum()
# yr_sum = yr_sum.reset_index()
# yr_sum.groupby('year').mean()*3.28*12 # m/yr, 0.67 inches/yr for max year
# we are underestimating recharge

# %%
# well_loc_merge.depth_m
## sample well layer
# get top elevations to calculate well screen elevation
well_loc_merge['top_elev'] = m.dis.top.array[well_loc_merge.row-1, well_loc_merge.column-1]
well_loc_merge['well_elev'] = well_loc_merge.top_elev - well_loc_merge.depth_m
# make layer 1 based
well_loc_merge['layer'] = get_layer_from_elev(well_loc_merge.well_elev.values, 
                    botm[:, well_loc_merge.row-1, well_loc_merge.column-1], m.dis.nlay)+1
# make 0-based columns
well_loc_merge[['k','i','j']] = well_loc_merge[['layer','row','column']]-1

# %% [markdown]
# Specify dates for the entire model period and seasonal dates to start/stop the model.

# %%
all_strt_date = pd.to_datetime(m.dis.start_datetime)
all_dates = all_strt_date + (m.dis.perlen.array.cumsum()-1).astype('timedelta64[D]')
all_end_date = all_dates[-1]
print(all_strt_date, all_end_date)
months = pd.date_range(all_strt_date, all_end_date, freq='MS')
years = pd.date_range(all_strt_date, all_end_date, freq='YS').year.values

# %%

# load summary excel sheet on irrigation optimization
# this will specify the date ranges to run and pause
fn = join(data_dir,'static_model_inputs.xlsx')
season = pd.read_excel(fn, sheet_name='Seasons', comment='#')

# %%
# choose crops on first day of year
# month_crop = pd.Series(1)
# day_crop = pd.Series(1)
# month_irr = pd.Series(4)
# day_irr = pd.Series(1)

# ## specify dates where modflow will start 
# all_run_dates = pd.DataFrame()
# # yn = 0
# # y = years[yn]
# for y in years:
#     run_dates = swb.ymd2dt(y, season.month_run, season.day_run, season.start_adj)
#     run_dates = run_dates.drop_duplicates().sort_values()
#     run_dates = pd.DataFrame(run_dates).assign(use='irrigation')
#     crop_date = swb.ymd2dt(y, month_crop, day_crop, season.start_adj)
#     crop_date = pd.DataFrame(crop_date).assign(use='crop').dropna()
#     irr_date = swb.ymd2dt(y, month_irr, day_irr, season.start_adj)
#     irr_date = pd.DataFrame(irr_date).assign(use='irrigation').dropna()
#     # all_run_dates = pd.concat((all_run_dates, crop_date, run_dates))
#     all_run_dates = pd.concat((all_run_dates, crop_date, irr_date))

# all_run_dates
# simple code to set dates for april 1
all_run_dates = pd.date_range(all_strt_date, all_end_date,freq='AS-Apr')
all_run_dates = pd.DataFrame(all_run_dates).assign(use='irrigation')
# add total start and end dates
all_run_dates = pd.concat((pd.DataFrame([all_strt_date]).assign(use='start'), all_run_dates))
all_run_dates = pd.concat((pd.DataFrame([all_end_date]).assign(use='end'), all_run_dates))
all_run_dates=all_run_dates.sort_values(0).reset_index(drop=True).rename(columns={0:'date'})

# %% [markdown]
# Review season dates, the plan was to change the start of the irrigation date for misc. grain and hay since irrigators don't typically start until summer even though it grows in winter. Simplify to crop choice Jan 1 and Irrig. run Apr 1
# - also at some point we may have discussed just doing them at the same time? in that case we just need to iterate over years

# %%
print(all_run_dates)

# %% [markdown]
# Cameron asked about cross-validation and hind-sight in the model. Maybe look at best way to re-run SWB or modflow after seeing the model run.

# %% [markdown]
# # Begin yearly iteration here 

# %% [markdown]
# For the first period we want to run things simply to get started. Thus we won't have output from the irrigation model.
# - check which soil budget data is available from the pre-existing dataset or if we should set up something simple to run on the hydrology of the current year assuming full irrigation on rain-ETc deficit.
#     - this should be a question for Yusuke, in the off-season of October to March we are going to assume no irrigation so there will no groundwater pumping. Just dispersed recharge.

# %%
# # this loop was set to run for the years of interest
for m_per in np.arange(3, all_run_dates.shape[0]-1):
# for m_per in [2]:
    m_strt = all_run_dates.iloc[m_per].date
    m_end = all_run_dates.iloc[m_per+1].date-pd.DateOffset(days=1)
    use = all_run_dates.iloc[m_per].use
    dates = pd.date_range(m_strt, m_end)

    # %%
    print(m_strt.date(),'to', m_end.date())

    # %%
    # year = 2015
    year = m_strt.year
    # year
    # crop='Corn'

# %% [markdown]
# # Define groundwater elevation sampling

    # %%
    model_ws_last = loadpth + 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date())
    hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')

    # %%
    # model output dates
    m_dates = all_run_dates.loc[m_per-1].date+np.array(hdobj.get_times()).astype('timedelta64[D]')
    m_dates = pd.DataFrame(m_dates, columns=['dates']).set_index('dates')
    m_dates['kstpkper'] = hdobj.get_kstpkper()
    # subset to 1 month of output
    # determine dates for fall sampling
    fall_dates = m_dates[m_dates.index.month==10]
    # determine dates for spring sampling
    spring_dates = m_dates[m_dates.index.month==3]


# %% [markdown]
# Irrigation decisions are made with spring levels but the crop choice model was developed with fall levels

    # %%
    # get head value from last 30 days to avoid using extreme single day value
    fall_heads = avg_heads(fall_dates.kstpkper.values, hdobj, m)
    
    # the dtw conversion runs a little slow
    # get the DTW for the wels in the simulation from the last period
    well_dtw = sample_dtw(fall_heads, botm)
    # need to make integer for join with crop choice
    well_dtw.UniqueID = well_dtw.UniqueID.astype(int)

# %% [markdown]
# # Crop choice model
# Evetnaully this should use the updated DTW from each previous year.

    # %%
    # load Sac Valley WYT
    wyt = pd.read_csv(join(proj_dir, 'model_inputs', 'sacramento_WY_types.txt'))
    # whether it is a critical or dry year or not
    wyt['wy_dc'] = 0
    # for wet or above normal years the boolean will be 0 
    wyt.loc[wyt['Yr-type'].isin(['C','D']),'wy_dc'] = 1

    # %%
    crop_choice_dir = 'parcelchoicemodelupdate'
    # Read logit coefficients
    logit_coefs = pd.read_csv(join(crop_choice_dir, 'logit_coefs.csv'))
    
    # the parcel data needs the dtwfa (avg fall dtw in feet for the parcel) and wy_dc (pulled from Sac wy type dataset and switched to dry boolean)
    # missing WY type prediction? 
    # Read parcel data
    data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
    # still needs to be updated to auto update DTW and WY type
    # data['wy_dc'] = np.where(data['year'] == 2020, 1, 0) # should be pulled from Sac WY type data instead
    data['wy_dc'] = wyt.loc[wyt.WY==year, 'wy_dc'].values[0]
    # update DTW to use simulated value instead of contoured
    well_dtw_merge = well_dtw[['UniqueID','dtw_ft']].rename(columns={'dtw_ft':'dtwfa'})
    data = data.drop(columns=['dtwfa', 'dtwsp'])
    data = data.merge(well_dtw_merge, left_on='parcel_id',right_on='UniqueID')
    
    # Import prior year revenue data by crop
    rev_prior_yr_df = pd.read_csv(join(crop_choice_dir, "data_model/rev_prior_yr.csv"))

# %% [markdown]
# Of note is that the simulated DTW for fall is 4x time greater than the observed so we may see a sharp drop in crops because of this.
# - it leads to an increase in unclassified fallow

    # %%
    # expect updated column name for pod
    # add column to POD that was available in previous dataset
    pod = data[['parcel_id','pod']].copy().rename(columns={'pod':'pod_bool'})
    pod['pod'] = 'No Point of Diversion on Parcel'
    pod.loc[pod.pod_bool==1, 'pod'] = 'Point of Diversion on Parcel'

    # %%
    # crop choice model uses "_" instead of " "
    # the irrigation model set up uses " "
    data_out = predict_crops(data.copy(), rev_prior_yr_df, logit_coefs)
    data_out['Crop_Choice'] = data_out.Crop_Choice.str.replace('_',' ')
    # update naming of Corn
    data_out.Crop_Choice = data_out.Crop_Choice.str.replace('Corn  ','Corn, ')
    
    data_out.to_csv(join(crop_choice_dir, 'parcel_crop_choice_'+str(year)+'.csv'), index=False)

# %% [markdown]
# # Load MF output
# The depth to water function should sample from the previous model run which may be a year or less long.
#
# The original get_dtw function was set up assuming a continuous model run which won't be the case.

    # %%
    # get head value from last 30 days to avoid using extreme single day value
    spring_heads = avg_heads(spring_dates.kstpkper.values, hdobj, m)
    
    # the dtw conversion runs a little slow
    # get the DTW for the wels in the simulation from the last period
    well_dtw = sample_dtw(spring_heads, botm)
    # need to make integer for join with crop choice
    well_dtw.UniqueID = well_dtw.UniqueID.astype(int)

# %% [markdown]
# # Irrigation submodel

# %% [markdown]
# The crop choice model didn't predict any corn.
# - this submodel may need to be in a separate script to run multiprocessing (parallel) as this re-loads the entire active script each time.

    # %%
    crop_list = ['Corn','Alfalfa','Pasture','Misc Grain and Hay', 'Grape']


    # %%
    #
    crop_in = data_out.merge(pod)
    crop_in = crop_in.rename(columns={'Crop_Choice':'name'})
    crop_in.to_csv(join(base_model_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'))
    # crop_in[crop_in.name==crop]
    pred_crops = crop_in.name.unique()
    print(pred_crops)

    # %%
    print(crop_in.groupby('name')[['parcel_id']].count())
    # paper by Snyder, statewide values for pasture or range

# %% [markdown]
# ## Simplified representation of DTW range from min to max.
# Instead of using the full record of the DTW, the mdoel should just sample the average DTW for the month of interest

    # %%
    ## simple representative DTW for linear steps from dtw_min to dtw_max
    ## with a 5 ft decline from June to December based on observed data
    dtw_simple_df = calc_simple_dtw(well_dtw, year, dtw_step = 10)
    dtw_simple_df.to_csv(join(loadpth, 'rep_crop_soilbudget','field_SWB', 'dtw_ft_WY'+str(year)+'.csv'))

# %%
# pd.DataFrame(dtw_simple[:,0], dtw_avg.index).loc['2020-02-12':'2020-10-04'].plot()
# plt.plot(dtw_simple.mean(axis=0))
# pd.concat([dtw_simple_df]*2, axis=1)
# dtw_simple_df.iloc[:, ::2]

# %% [markdown]
# ## No irrigation SWB
# We will need to use a script similar to reference_swb.py since we need the recharge rates outside of the typical irrigation season and for fallow fields. For simplicity we can run this for the full model period (Apr to Apr) for all the selected crops.

# %%
# ## simple soil water budget without optimization to calculate
# ## recharge rates when there is no irrigation

# # create a simple array of zero depth to water to inform the soil water budget
# # without optimization
dtw_zero = pd.DataFrame(np.zeros((len(dtw_simple_df), len(data_out))))
# hopefully the parcel ids match as datatype (int vs str)
dtw_zero.columns = data_out.parcel_id.values
dtw_zero.index = dtw_simple_df.index

# for crop in crop_list:
#     var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
#     # need to account for when crops aren't predicted and skip them
#     if pred_dict[crop] in pred_crops: 
#         # to equalize the situation we might use a simple DTW profile

#         load_run_swb(crop, year, crop_in, join(loadpth, 'crop_soilbudget_no_irr'),
#                      dtw_zero, 
#                      soil_rep=False,
#                     run_opt=False) 

# %% [markdown]
# ## Iterate over all crops to save the representative profiles
# The way the irrigation models are set up, they still run for the entire season but they are simulated in spring which means that they can't decide winter water pumping for the first year. In each year following the simulation runs april to april so they can influence pumping in the following year.
# - just need to note that the first period of 6 months is spin-up in a way.

# %% [markdown]
# After updating function need to check everything still works.
#
# - fixed issue hdf5 overwrite
# - found that Corn, sorghum name match was off

    # %%
    print(pred_crops)
    for crop in crop_list:
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
        # need to account for when crops aren't predicted and skip them
        # if pred_dict[crop] in pred_crops: 
        print(crop, ':',pred_dict[crop])

    # %%
    # initialize HDF5 files for the year
    from functions.data_functions import init_h5
    base_model_ws = join(loadpth, 'rep_crop_soilbudget')
    # initialize SWB folder
    os.makedirs(join(base_model_ws, 'field_SWB'), exist_ok=True)
    for var in ['profit', 'yield', 'percolation','GW_applied_water', 'SW_applied_water']:
        name = join(base_model_ws, 'field_SWB', var + '_WY'+str(year)+'.hdf5')
        init_h5(name)
    
    # check profit/yield if annual or daily

# %%
    
    # for crop in ['Corn']:
    for crop in crop_list:
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
        # need to account for when crops aren't predicted and skip them
        if pred_dict[crop] in pred_crops: 
            # to equalize the situation we might use a simple DTW profile
            load_run_swb(crop, year, crop_in, join(loadpth, 'rep_crop_soilbudget'),
                         dtw_simple_df, 
                         soil_rep=True) 

    # %%
    fn = join(model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
    print('Crops with SWB results')
    with h5py.File(fn) as dset:
        finished_crops = list(dset['array'].keys())
        print(finished_crops)
    # only grape was completed?

# %% [markdown]
# 1. Load the representative results and sample for each field by crop type to back calculate the irrigation requirements. use the estimated irrigation as an input to the modflow model for pumping and percolation for recharge.
#     -   use the DTW id to reference to the irrigation in the full array, need to group by SW, GW or mixed.
#     -   if we wanted we could re-run the SWB one time with the specified irrigation rates to get the exact recharge rates with field specific values
# 2. RUn the modflow model to get the resultant DTW profile
# 3. re-calculate the profit using the irrigation and actual DTW profile on a soil by soil basis (non-optimization) after running the next modflow chunk. Actually the re-run for the true profit could be done if profits aren't needed mid-simulation

    # %%
    model_ws = join(loadpth, 'rep_crop_soilbudget')
    # load the processed dataframe with all datas
    pc_df_all, irr_gw_df_all, irr_sw_df_all = get_wb_by_parcel(model_ws, crop, year, 
                     crop_in, finished_crops, dtw_simple_df, well_dtw)

# %% [markdown]
# ## RCH input

    # %%
    # join water budget data to field area for scaling to cell for recharge
    rch_df = pc_df_all.merge(field_df)
    # calculate the fraction of a cell covered by the field
    rch_df['field_scale'] = rch_df.field_area/(200*200)
    # calculate the recharge rate for each field
    rch_df['rch_rate'] =  rch_df.rate*rch_df.field_scale
    # set date as index for sampling
    rch_df = rch_df.set_index('date')
    # we only need to assign rch when non-zero
    rch_df = rch_df[rch_df.rch_rate!=0]

# %% [markdown]
# Take the combined dataframe to inform the recharge array with recharge from the irrigation season. The next step is to insert recharge from non-irrigated times and native lands.

    # %%
    ## add native lands recharge
    rch_all = pd.concat((rch_df, native_pc))
    # aggregate to the row, column level
    rch_all = rch_all.groupby(['date','row','column'])[['rch_rate']].sum()
    rch_all = rch_all.reset_index(['row','column'])
    # is groupby faster or just adding individual datasets to the array??

# %%
## add off-season farmlands recharge



    # %%
    # allocate dictionary for recharge
    rch_dict = dict()
    for t, d in enumerate(dates):
        # get data for the stress period
        rch_arr = np.zeros((m.dis.nrow, m.dis.ncol))
        if d in rch_df.index:
            # spd_df = rch_df.loc[d]
            spd_df = rch_all.loc[d]
            rch_arr[spd_df.row-1, spd_df.column-1] = spd_df.rch_rate
        # assign values to recharge dict
        rch_dict[t] = rch_arr

# %% [markdown]
# ## WEL input

# %%
    
    # join water budget data to field area for application to cell for well
    wel_df = irr_gw_df_all.merge(well_loc_merge)
    # calculate the volumetric flux of the well 
    wel_df['flux'] = -wel_df.rate*wel_df.area_m2
    # # set date as index for sampling
    wel_df = wel_df.set_index('date')
    # wel_df

# %%
# allocate dictionary for recharge
wel_dict = dict()
for t, d in enumerate(dates):
    # get data for the stress period
    if d in wel_df.index:
        spd_df = wel_df.loc[d]
        wel_arr = spd_df[['k','i','j','flux']].values
    else:
        wel_arr = [[0,0,0,0]]
    # append domestic well pumping to pumping dataset
    dom_loc['flux'] = - dom_use.loc[d,'flux_m3d']*dom_loc.pump_scale
    wells_dom = dom_loc[['layer','row','column','flux']].values
    
    # assign input to well dictionary
    wel_dict[t] = np.append(wel_arr, wells_dom, axis=0)


# %%
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# plt.imshow(rch_dict[100])
# plt.colorbar(norm = mpl.colors.LogNorm)

# %% [markdown]
# # Run update modflow

    # %%
    # load the existing modflow model for the next year
    mf_ws = loadpth + 'crop_modflow/'+str(m_strt.date())
    m_month = flopy.modflow.Modflow.load(join(mf_ws,'MF.nam'), model_ws = mf_ws, load_only=['DIS'])
    # m.model_ws = mf_ws
    # The model should start in hydraulic connection
    if m_per > 0:
        mf_ws_last = loadpth + 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date())
        hdobj = flopy.utils.HeadFile(mf_ws_last + '/MF.hds')
        sp_last = hdobj.get_kstpkper()[-1]
        strt = hdobj.get_data(sp_last)
    
    # for the first period we use dem as starting point
    # if solver criteria are not met, the model will continue if model percent error is less than stoperror
    bas_month = flopy.modflow.ModflowBas(model = m_month, ibound=ibound, strt = strt)
    # write bas input with updated starting heads
    bas_month.write_file()

# %%

    # update flopy packages with new data
    wel_month = flopy.modflow.ModflowWel(model=m_month, stress_period_data = wel_dict, ipakcb=55)
    rch_month = flopy.modflow.ModflowRch(model=m_month, nrchop = 3, rech =  rch_dict, ipakcb=55)
    
    wel_month.write_file()
    rch_month.write_file()
    


    # %%
    success, buff = m_month.run_model()

# %% [markdown]
# The MODFLOW run-time for one year was 49 min. This was for the version with 19 layers which we can cut down to 10 or so by removing TPROGs and this could be sped up by reducing the amont of output saved since we may not need to review the CBC output.

# %% [markdown]
# # Re-run soil water budget to calculate actual profit
# The soil water budget needs to be re-run at the end of each year to calculate the profit based on the actual water availability for surface water and groundwater. This could be done at the very end of the simulation if the farmer doesn't need previous year profits between years.
