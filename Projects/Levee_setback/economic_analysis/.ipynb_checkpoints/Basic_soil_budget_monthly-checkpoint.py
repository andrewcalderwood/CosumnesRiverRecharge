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
# Code to build a basic soil water budget model based on the mathematical equations used by IDC with simplification as here the key is identifying percolation rather than applied water.

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import glob
import sys
import time
import h5py

import pandas as pd
import numpy as np

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import shapely
import geopandas as gpd

# specify variables needed by simulation
global var_gen, var_crops_all, var_yields_all, crop_dict, pred_dict
global season_all, rain_m, ETo_m, Kc, Kc_dates, soil_ag_all

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%

loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Regional/'

model_ws = loadpth + 'crop_modflow'

# %%
def ymd2dt(year, month, day, year_adj):
    """ Only year is expected to come as int the others as series/1D arrays"""
    year = (year + year_adj).astype(str)
    date = pd.to_datetime(year+'-'+ month.astype(str)+'-'+day.astype(str))
    return(date)



# %% [markdown]
# ## Hydrologic data
# Could be shared globally

# %%
# pd.read_excel(fn, sheet_name=None, comment='#').keys()


# %%
# var_crops = var_crops_all[var_crops_all.crop==crop].copy()
# year=2017
def get_var_year(df, year):
    """
    Given a crop variable with different values by year
    select the closest year or use an average, in a tie it uses the lower year
    """
    # df = var_crops_all[var_crops_all.crop==crop].copy()

    # find variables that change by year
    df_year = df[~df.year.isna()]
    if year is not None:
        # find the year closest to the simulated year
        df_year = df_year.loc[[(df_year.year-year).abs().idxmin()]]
    else:
        # if no year is supplied then use the average of the data
        df_year = df_year.dropna(axis=1, how='all') # can't group by columns of all NAs
        grp_cols = df_year.columns[~df_year.columns.isin(['year','value'])].tolist()
        df_year = df_year.groupby(grp_cols).mean(numeric_only=True).reset_index()
    # add year data back to the static data
    df = pd.concat((df[df.year.isna()], df_year))
    # since the season exists for Alfalfa/Pasture sort to make sure the variables are applied correctly
    df = df.sort_values(['variable', 'season'])
    return df

# year = None
# get_var_year(var_crops, year)


# %%
def load_var(crop, year=None):
    # crop = 'Alfalfa'
    fn = join(data_dir,'static_model_inputs.xlsx')
    var_gen = pd.read_excel(fn, sheet_name='General', comment='#')
    var_gen = var_gen.set_index('variable')['value'] # adjust for quick pulling of variables

    # new version has prices by year
    # var_crops_all = pd.read_excel(fn, sheet_name='Crops', comment='#').dropna(axis=1, how='all')
    var_crops_all = pd.read_excel(fn, sheet_name='Crops_byyear', comment='#').dropna(axis=1, how='all')
    var_yields_all = pd.read_excel(fn, sheet_name='Yield', comment='#').dropna(axis=1, how='all')
    
    # cross-reference between simple and DWR full crop name
    crop_dict = pd.read_excel(fn, sheet_name='Name_dict', comment='#')
    pred_dict = crop_dict.set_index('Crop')['pred_name'].to_dict()
    crop_dict = crop_dict.set_index('Crop')['DWR_name'].to_dict()
    
    # season dates - primarily for calculating yield adjustments
    season_all = pd.read_excel(fn, sheet_name='Seasons', comment='#')
    
    # subset to crop for current simulation if looping, could also do ID
    var_crops = var_crops_all[var_crops_all.crop==crop]
    var_crops = get_var_year(var_crops, year)
    var_crops = var_crops.set_index('variable')['value'] # adjust for quick pulling of variables
    var_yield = var_yields_all[var_yields_all.crop==crop]
    
    season = season_all[season_all.crop==crop].sort_values(['crop','season'])
    return(var_gen, var_crops, var_yield, season, pred_dict, crop_dict)

# %%
# get_var_year(var_crops, 2020)


# %%
def load_hyd(year, dates):
    ## Potential ETo spatial interpolation from CIMIS
    fn = glob.glob(join(uzf_dir,'CIMIS','Cosumnes_dailyET_precip*.csv'))
    daily_data = pd.DataFrame()
    for file in fn:
        new_data = pd.read_csv(file, index_col = ['Date'], parse_dates = True)
        daily_data = pd.concat((daily_data, new_data))
    # units of mm
    data_in = daily_data[daily_data['Stn Name']=='Fair Oaks']
    
    # clean up data so columns are by location, units of Precip are in mm
    rain_in = data_in.pivot_table(index = 'Date', columns = 'Stn Name', values = 'Precip (mm)')
    rain_m = rain_in/1000
    
    # clean up ETo data
    ETo_in = data_in.pivot_table(index = 'Date', columns = 'Stn Name', values = 'ETo (mm)')
    ETo_m = ETo_in/1000

    # fill in missing days of rain with zero values
    rain_df = rain_m.reindex(dates)
    rain_df[rain_df.isna()]=0
    # convert to array for calculations
    rain = rain_df['Fair Oaks'].values
    
    # linearly interpolate missing days for ET
    ETo_df = ETo_m[dates.min():dates.max()].resample('D').interpolate('linear')['Fair Oaks']
    # a zero value of ETo caused a divide by zero fault for pasture in the yield calculation
    # smallest value in 2015 is 4E-5, assume min value of 1E-9 to avoid divide by 0
    ETo_df[ETo_df==0] = 1E-9
    
    return(rain, ETo_df)


# %%
def load_Kc(year):
    fn = join(data_dir,'static_model_inputs.xlsx')
    Kc = pd.read_excel(fn, sheet_name='Kc', comment='#')
    Kc = Kc.set_index('Period')
    Kc_dates = pd.read_excel(fn, sheet_name='Kc_dates', comment='#')
    Kc_dates = Kc_dates[Kc_dates.Growth_stage!='Growing season'] # don't need overall dates for Kc calculation
    Kc_dates['Start_date'] = pd.to_datetime({'year': year+Kc_dates['start_adj'], 'month':Kc_dates['Start_month'], 'day': Kc_dates['Start_day']})
    Kc_dates['End_date'] = pd.to_datetime({'year': year+Kc_dates['end_adj'], 'month':Kc_dates['End_month'], 'day': Kc_dates['End_day']})
    Kc_dates.Start_date = pd.to_datetime(Kc_dates.Start_date)
    Kc_dates.End_date = pd.to_datetime(Kc_dates.End_date)
    Kc_dates = Kc_dates.set_index('Growth_stage')
    return Kc, Kc_dates

# %%
def get_Kc_dates(Kc_dates_c, Kc_c):
    """ Given a dataframe with the dates for Kc periods and a dataframe with those Kc values create a timeseries
    of Kc values"""
    Kc_df_all = pd.DataFrame()
    # initial period and mid season are constant, crop dev and late season are linear
    # initial period
    dates_p = pd.date_range(Kc_dates_c.loc['Initial period', 'Start_date'], Kc_dates_c.loc['Initial period', 'End_date'])
    Kc_df = pd.DataFrame(dates_p, columns=['date']).assign(Kc=Kc_c.loc['ini', 'Value'])
    Kc_df_all = pd.concat((Kc_df_all, Kc_df))
    # crop development
    dates_p = pd.date_range(Kc_dates_c.loc['Crop development', 'Start_date'], Kc_dates_c.loc['Crop development', 'End_date'])
    Kc_df = pd.DataFrame(dates_p, columns=['date'])
    Kc_df['Kc'] = np.linspace(Kc_c.loc['ini', 'Value'], Kc_c.loc['mid', 'Value'], len(Kc_df))
    Kc_df_all = pd.concat((Kc_df_all, Kc_df))
    # mid-season
    dates_p = pd.date_range(Kc_dates_c.loc['Mid-season', 'Start_date'], Kc_dates_c.loc['Mid-season', 'End_date'])
    Kc_df = pd.DataFrame(dates_p, columns=['date']).assign(Kc=Kc_c.loc['mid', 'Value'])
    Kc_df_all = pd.concat((Kc_df_all, Kc_df))
    # late season
    dates_p = pd.date_range(Kc_dates_c.loc['Late season', 'Start_date'], Kc_dates_c.loc['Late season', 'End_date'])
    Kc_df = pd.DataFrame(dates_p, columns=['date'])
    Kc_df['Kc'] = np.linspace(Kc_c.loc['mid', 'Value'], Kc_c.loc['end', 'Value'], len(Kc_df))
    Kc_df_all = pd.concat((Kc_df_all, Kc_df))
    return(Kc_df_all)



# %% [markdown]
# # Spatial data
# Yusuke's updated data set has UniqueID for each parcel. This must then be joined with the predicted crop type.
#
#

# %%
def load_soil(crop, crop_in):
    ''' Imports soil data for all fields and subsets to those for the parcels
    with the current crop
    crop: string specifying the crop
    crop_in: dataframe with parcel_id, name (crop), pod'''
    field_ids = 'parcel' # 'ag'
    # load cleaned soil data for ag fields
    soil_path = join(uzf_dir,'clean_soil_data')
    # soil data for each ag field
    soil_ag_all = pd.read_csv(join(soil_path, 'soil_for_'+field_ids+'_fields.csv'), index_col=0)

    # curve numbers
    CN = pd.read_csv(join(soil_path, field_ids+'_field_CN.csv'),index_col = 0)
    
    soil_ag_all = pd.merge(soil_ag_all, CN)
    
    # add crop data to fields, NAs are introduced for when a crop isn't specified
    soil_ag_all = soil_ag_all.merge(crop_in, left_on='UniqueID', right_on='parcel_id')
    # for now drop fields without crop specified 
    soil_ag_all = soil_ag_all[~soil_ag_all.name.isna()]
    
    # subset for fields with the crop (slight differences in naming occurred switching around datasets (DWR to LandIQ)
    # soil_crop = soil_ag_all[soil_ag_all.name==crop]
    soil_crop = soil_ag_all[soil_ag_all.name==crop]
    return(soil_crop)


# %%
def prep_soil(soil_ag, etc_arr, var_crops):
    global soil_Ks, soil_por, soil_eps, soil_CN
    global soildepth, soil_m, wc_f, wc_wp, taw, Smax, raw, P
    
    # # when soil_K_low is missing using a substitute of Ksat/10
    soil_Ks = np.where(soil_ag.Ksat_Low==0, soil_ag.Ksat/10, soil_ag.Ksat_Low)
    soil_por = soil_ag.Porosity.values/100
    soil_eps = soil_ag.EPS.values
    soil_CN = soil_ag.CN.values
    
    soildepth = soil_ag.SoilDepth.values
    psdi =  soil_ag.PSDI.values
    # parameter for Mualem, van Genuchten
    soil_m = psdi/(psdi+1)
    wc_f =  soil_ag.w3rdbar.values/100 #field capacity
    wc_wp =  soil_ag.w15bar.values/100 #wilting point 

    # Calculate total available water in the root zone
    taw = (wc_f - wc_wp)*soildepth 

    # for runoff, convert CN from fraction to CN
    Smax = (1000/soil_CN) - 10
    
    p_table22 = var_crops['p_table22'] # Soil water depletion fraction for no stress
    # requires ETc in mm/day
    P = p_table22 + 0.04*((5-(etc_arr*1000))) # Calculate adjusted daily soil water depletion fraction for no stress
    raw = taw*P # Calculate readily available water in the root zone

# %% [markdown]
# Total available water (TAW) is the water available between field content and the wilting point times the depth of the root zone (this could be assumed to be the soil depth at times).  
# Per the FAO report, rooting depths are:
# - Corn:1-1.7 m which is about where soil depth would put it  
# - Alfalfa:1-2 m  
# - Barley, Wheat: 1-1.5m  
#
# So all of these fall within the range where the soil depth should determine rooting depth.

# %% [markdown]
# If we need to remove ET from each step then we should reduce the number of crops used to simplify calculations unless the code runs quickly for all cells.
# - Need to assume initial water content of 0. 
#  
# * water content can't exceed porosity. It is assumed that θ can be computed to be greater than θT, and the difference between the θ and θT represents the depth of the pond.

# %%

# %% [markdown]
# # Run the model
# The ETc, and rain will be preset inputs for the irrigation simulation (the ETc will be calculated when the crop is chosen at the start of the year).  
# The irrigation optimizer will be a function of crop - ETc, irrigation (variable calibrated)
#
#

# %% [markdown]
# After the irrigation season ends, the model needs to be run from the end of irrigation to the start of the next season for calculating recharge from rainfall and maintaining the soil water budget.  
#
# The FAO model assumes "that water can be stored in the root zone until field capacity is reached. Although following heavy rain or irrigation the water content might temporally exceed field capacity, the total amount of water above field capacity is assumed to be lost the same day by deep percolation, following any ET for that day. By assuming that the root zone is at field capacity following heavy rain or irrigation, the minimum value for the depletion Dr, i is zero."  
#
# -> this contrasts with IDC where water content is allowed to exceed porosity with the assumption that it has become ponded water.

