# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Script to run the code for the soil water budget model simulated for representative soil profiles and DTW values for cases with and without points of diversion. The soil profiles should ideally be for each map unit but for simplicity we have been using the average for fields selected, now that it is independent of crop choice we will need to average across the entire domain of parcels or by soil map unit. The DTW will be across the practical range of DTW values and apply a uniform assumption of 10 ft of decline in the summer although for simplicity we may want to keep this static to reduce the number of assumptions we make?

# %%
import sys
import os
from os.path import join,exists, dirname, basename, expanduser
import time 
import re

import h5py
import numpy as np
import pandas as pd
import geopandas as gpd

from importlib import reload

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')

# %%
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

add_path('../')


# %%
import functions.Basic_soil_budget_monthly as swb
# reload(swb)

# import f_rep_swb_profit_opt
# reload(f_rep_swb_profit_opt)
from f_rep_swb_profit_opt import load_run_swb


# %%
# import functions.f_gw_dtw_extract
# reload(functions.f_gw_dtw_extract)
from functions.f_gw_dtw_extract import sample_dtw, avg_heads
# from functions.f_gw_dtw_extract import get_dtw
from functions.f_gw_dtw_extract import calc_simple_dtw


# %%
# from functions import data_functions
# reload(functions.data_functions)
from functions.data_functions import read_crop_arr_h5

from functions.data_functions import init_h5


# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
in_data = sys.argv

if 'ipykernel' in in_data[0]:
    # input_name = 'static_model_inputs.xlsx'
    input_name = 'static_model_inputs_no_p_o.xlsx'
    year_load_var_in = "False"

else:
    # added option to specify different static_model_inputs to allow easier adjustment of operating costs, revenue
    input_name = in_data[1]
    year_load_var_in = in_data[2]

print('sys.argv[1] (input_name) is...')
print(input_name)

print('sys.argv[2] (year_load_var_in) is...')
print(year_load_var_in)
print('\n')

swb_version = input_name.replace('static_model_inputs','').replace('.xlsx','')

year_load_var_in = "False"

if year_load_var_in != "False":
    year_load_var_in = int(year_load_var_in)
    swb_version += '_'+year_load_var_in

print('swb version:',swb_version)


# %%
swb_ws = join(proj_dir, 'model_inputs', 'swb_rep', 'version'+swb_version)
os.makedirs(swb_ws, exist_ok=True)

# %%
# save log by date so we can see old versions
os.makedirs(join(swb_ws, 'log'), exist_ok=True)
if 'ipykernel' in in_data[0]:
    print('print to jupyter')
else:
    print('printing remaining output to log file')
    sys.stdout = open(join(swb_ws, 'log', 'rep_SWB_log_'+str(pd.to_datetime('today').date())+'.txt'), 'w')

print('Run started at:', pd.to_datetime('today'))
print('swb version:',swb_version)


# %%
# loda model scenario reference sheet
scenario_summary = pd.read_excel(join(data_dir, 'scenario_summary.xlsx'))
# all scenarios use constraints of 100 inches except for the gw pumping constraint test of 36 for one
sw_con = 100
gw_con = 100
# # print out additional summary info
print('SW Constraint %.2f' %sw_con, 'inches')
print('GW Constraint %.2f' %gw_con, 'inches')


# %%
t_start = time.time()


# %%
# load parcel data for reference as needed
parcels = gpd.read_file(join(proj_dir,'Parcels shapefile/parcels.shp'))
parcels['area_m2'] = parcels.geometry.area
parcels.UniqueID = parcels.UniqueID.astype(int)

field_ids = 'parcel'
soil_path = join(uzf_dir,'clean_soil_data')
grid_soil_ag = pd.read_csv(join(soil_path, field_ids+'_field_to_cell.csv'), index_col=0)


# %%
crop_choice_dir = '../parcelchoicemodelupdate'
data = pd.read_csv(join(crop_choice_dir, "data_model/parcel_data_test.csv"))
# add column to POD that was available in previous dataset
pod = data[['parcel_id','pod']].copy().rename(columns={'pod':'pod_bool'})
pod['pod'] = 'No Point of Diversion on Parcel'
pod.loc[pod.pod_bool==1, 'pod'] = 'Point of Diversion on Parcel'


# %%
# crop_in isn't creatd from prediction so for each step we average properties across all parcels
crop_in_parcels = parcels[['UniqueID']].rename(columns={'UniqueID':'parcel_id'})
crop_in_parcels = crop_in_parcels.sort_values('parcel_id').reset_index(drop=True)
crop_in_parcels = crop_in_parcels.merge(pod)

# %%
# load cleaned soil data for ag fields
soil_path = join(uzf_dir,'clean_soil_data')
# soil data for each ag field
soil_ag_all = pd.read_csv(join(soil_path, 'soil_for_'+field_ids+'_fields.csv'), index_col=0)
soil_ag_gdf = parcels[['UniqueID','geometry']].merge(soil_ag_all,on='UniqueID')
# soil_ag_gdf.plot('Texture',legend=True,  legend_kwds = {'loc':(0.01,0.7)})
# the parcels are predominantly clay loam or loam
# soil_ag_gdf.plot('Texture',legend=True,  legend_kwds = {'loc':(0.01,0.7)})
# soil_ag_gdf.hist();
# there is not a huge spread for results
# field content is from 20-30 mostly
# ksat is prety much uniformly below 1
# porosity is between 35-50
# likely minimal impacts if we assume average conditions
# could acknowledge that while some farmers might account for this, we are focusing
# on the impact of DTW and

#the main adjustment that could be made is subsetting to historical fields for a crop
# but this would add some potential bias so safer to use average across all cultivated fields

# %%
# load summary excel sheet on irrigation optimization
# this will specify the date ranges to run and pause
fn = join(data_dir ,input_name) # new version allows alternates

season = pd.read_excel(fn, sheet_name='Seasons', comment='#')
irr_est = pd.read_excel(fn, sheet_name='Irr_establish', comment='#')

# read in crop choice parameter csv for list of possible crops to run
# at this point we may want to try to standardize names a bit more
# if trying to clean up the files
crop_dict = pd.read_excel(fn, sheet_name='Name_dict', comment='#')
# names in the static_model_inputs sheet are the short form
# spreadsheet Name_dict does the translation from parameter to DWR to crop prediction name
# the others are Other (no parameters) and fallow
crop_list = crop_dict.Crop.values


# %%
all_strt_date = pd.to_datetime('2014-10-1')
all_end_date = pd.to_datetime('2025-9-30')

# simple code to set dates for april 1
all_run_dates = pd.date_range(all_strt_date, all_end_date,freq='YS-APR')
all_run_dates = pd.DataFrame(all_run_dates).assign(use='irrigation')
# add total start and end dates
all_run_dates = pd.concat((pd.DataFrame([all_strt_date]).assign(use='start'), all_run_dates))
all_run_dates = pd.concat((pd.DataFrame([all_end_date]).assign(use='end'), all_run_dates))
all_run_dates=all_run_dates.sort_values(0).reset_index(drop=True).rename(columns={0:'date'})
all_run_dates.to_csv(join(swb_ws, 'all_run_dates.csv'))

# %%
for m_per in np.arange(1, all_run_dates.shape[0]-1):
# for m_per in [1]:
# for m_per in [4]:

    m_strt = all_run_dates.iloc[m_per].date
    m_end = all_run_dates.iloc[m_per+1].date-pd.DateOffset(days=1)
    use = all_run_dates.iloc[m_per].use
    dates = pd.date_range(m_strt, m_end)

    # %%
    print(m_strt.date(),'to', m_end.date())

    # %%
    year = m_strt.year
    year_load_var = year

    # %%
    ## simple representative DTW for linear steps from dtw_min to dtw_max
    ## with a 5 ft decline from June to December based on observed data
    # this part then evaluates the min and max range of DTW
    # a 300 foot decline should be the upper limit as there would be significant well impacts
    # beyond this point as quite a few wells are 200-300 ft
    well_dtw = pd.DataFrame([0,300], columns=['dtw_ft'])
    # and applies a simple 5 ft decline from spring outward
    dtw_simple_df = calc_simple_dtw(well_dtw, year, dtw_step = 10)
    # for simplicity it may be better to simply use constant DTW values for each year

    # %%
    print(crop_list)
    for crop in crop_list:
    # for crop in ['Alfalfa']:
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year = year_load_var, input_name=input_name)
        # need to account for when crops aren't predicted and skip them
        print(crop, ':',pred_dict[crop])

    # %%
    print('SWB workspace:',swb_ws)
    print('Initializing hdf5 files')
    # maybe for these we could move toward a different array for each year?
    # get started first then adapt

    # %%
    # initialize HDF5 files for the year
    # TODO static crop/AW: if skipping SWB optimization then don't need to intialize these
    # initialize SWB folder
    for var in ['profit', 'cost', 'yield', 'percolation','GW_applied_water', 'SW_applied_water']:
        name = join(swb_ws,  var + '_WY'+str(year)+'.hdf5')
        init_h5(name)

    for var in ['swb_output']:
        name = join(swb_ws, var + '_WY'+str(year)+'.hdf5')
        init_h5(name, groups=['wc','ETa', 'rp'])

    # %%
    print('Start irrigation optimization run\n\n')
    sys.stdout.flush()

# %%
# missing variable crop_in (parcels defined with each crop)
# we shouldn't need these so may need to specify as None or dummy?

    # soil_crop = swb.load_soil(pred_dict[crop], crop_in)
    # if soil_rep:
    #     # get the average field properties for a representative field
    #     # with one for a point of diversion and one without
    #     soil_crop = soil_crop.groupby('pod').mean(numeric_only=True).reset_index()
    #     # repeat soil_crop data to iterate over different DTW profiles
    #     soil_crop = pd.concat([soil_crop]*dtw_df.shape[1])
# crop_in identifies parcels crop
# then load_soil samples them for each parcel that has the crop of interest

# %%
    
    # for crop in ['Alfalfa']:
    # TODO static crop/AW: if skipping SWB optimize then skip this
    for crop in crop_list:
        print(crop)
        # assign name for all parcels to crop of interest for that step in the loop
        crop_in = crop_in_parcels.copy()
        # will need to add year to swb.load_var(crop, year) if we want to use year specific profit and cost
        # variables, this would be useful for comparing against baseline while future should use average
        # within load_run_swb, these variables are called again specified by year
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop, year = year_load_var, input_name=input_name)
        crop_in['name'] = pred_dict[crop]
        # to equalize the situation we might use a simple DTW profile
        load_run_swb(crop, year, crop_in, swb_ws,
                     dtw_simple_df, 
                     soil_rep=True,
                     sw_con=sw_con, gw_con=gw_con,
                     input_name = input_name
                     ) 
        sys.stdout.flush()

# %%
t_final = time.time()
print('Total time was %.2f hours' %((t_final-t_start)/3600))

