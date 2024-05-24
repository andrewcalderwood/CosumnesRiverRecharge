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
# Goal: Run scripts needed to extend the model into the future.
#
# Data downloads:
# - CIMIS potential evapotranspiration and rainfall (soil water budget)
# - USGS Streamflow at Michigan Bar (stream inflow)
# - DWR groundwater elevation data base (lateral boundaries and observation matching)
# - DWR stream stage at McConnell (observation matching)
#
# Scripts to run:
# - SFR/ add script to download USGS data, I have an example of this somewhere
# - UZF/Basic soil budget - fields, check if gridded is still needed for natives?
#     - need to download another year of CIMIS precip/rainfall
# - GHB/ need to run scripts 1 to 3 but these need updates to move testing code to separate notebooks

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import sys
import glob
from importlib import reload

import pandas as pd
import numpy as np


# standard geospatial python utilities
# import pyproj # for converting proj4string
# import shapely
import geopandas as gpd
import rasterio



# %%
from dataretrieval import nwis


# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
print(gwfm_dir)
sfr_dir = gwfm_dir+'/SFR_data/'


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')

import flopy 

from importlib import reload
# importlib.reload
# reload(flopy)

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

from mf_utility import get_layer_from_elev, param_load


# %% [markdown]
# # Load data files
#
# streamflow and precipitation

# %% [markdown]
# ## streamflow

# %%
strt_date = '2000-1-1'
# strt_date = '2019-1-1'
end_date = '2023-12-31'

# %%
## specify the USGS site code for which we want data. Michigan Bar
site = '11335000'

# get instantaneous values (iv)
# df = nwis.get_record(sites=site, service='iv', start='2017-12-31', end='2018-01-01')
mb = nwis.get_record(sites=site, service='dv', start=strt_date, end=end_date, parameterCd='00060')

# get basic info about the site
mb_meta = nwis.get_record(sites=site, service='site')

# %%
mb.index = mb.index.tz_localize(None)


# %%
mb_out = mb.rename(columns={'00060_Mean':'flow_cfs'})[['flow_cfs']]
mb_out.to_csv(join(gwfm_dir, 'SFR_data', 'MB_daily_flow_cfs.csv'))

# %% [markdown]
# ## climate
# the cimis API has a max of 1750 data points per request, 3 years is 1094 steps and 4 years is 1460 steps

# %%
strt_date = '2021-1-1'
end_date = '2023-12-31'

# %%
# %%
# should pull from text file
cimis_key = '43cb8818-53b1-431d-b73b-f2d94e98ea6c'


# %%
# example URL
# http://et.water.ca.gov/api/data?appKey=YOUR-APP-KEY&targets=2,8,127&startDate=2010-01-01&endDate=2010-01-05


# %%
cimis_stns = pd.read_csv(join(gwfm_dir, 'UZF_data','CIMIS','CIMIS_station_lat_long.csv'))


# %%
from requests.auth import HTTPBasicAuth
import requests

def get_cimis_data(sites, strt_date, end_date, cimis_key, 
                   data = ['day-precip', 'day-eto'], unit='M'):
    """
    pull from the CIMIS API for ETo and Precip
    INPUT:
    sites: name of stations of interest [str]
    strt_date, end_date: start and end of period to sample [str]
    cimis_key: API key for the user [str]
    data: list of dataItems to sample from the API (e.g., day-precip)
    unit: M is metric, E is english
    OUTPUT:
    cimis_json: json format of the cimis data requested
    """
    # url = 'http://et.water.ca.gov/api/station'
    url = 'http://et.water.ca.gov/api/data'
    targets = ','.join(sites)
    dataItems = ','.join(data)
    ndays = (pd.to_datetime(end_date) - pd.to_datetime(strt_date) ).days
    if ndays >1750:
        print('No. days exceeds max CIMIS days')
        
    params = {'appKey': cimis_key,'unitOfMeasure':unit,
              'targets': targets, 'startDate': strt_date, 'endDate': end_date,
              # defaults to daily data for all parameters
             'dataItems':dataItems}
    # sample the url
    req = requests.get(url, 
                       params = params,
                       verify=False
                      )
    # pull json format data
    cimis_json = req.json()
    return(cimis_json)



# %%

# %%
sites = np.array([131]).astype(str)
cimis_json = get_cimis_data(sites, strt_date, end_date, cimis_key)

# %%
# manual way to create link for download
# url_param = ''
# for k in params.keys():
#     url_param +=  k+'='+params[k]+'&'

# url_param = url_param[:-1]
# url+'?'+url_param


# %%
# cimis_json['Data']['Providers'][0]['Records'][0].keys()

# %%
def clean_cimis_json(cimis_json):
    """ convert cimis json to long dataframe format """
    # standard cleaning for CIMIS data
    # default is ['Data']['Providers'][0]['Records']
    cimis_df = pd.DataFrame(cimis_json['Data']['Providers'][0]['Records'])
    # identify columns with data values
    val_cols = cimis_df.columns[cimis_df.columns.str.contains('day|hly', case=False)].values
    # convert data to long format to more easily clean up
    cimis_long = cimis_df.melt(id_vars=['Date','Station'], value_vars=val_cols)
    # convert columns with dictionary of values and units to columns
    cimis_long_values = pd.json_normalize(cimis_long['value'])
    # join cleaned value data to date
    cimis_clean = pd.concat((cimis_long[['Date','Station','variable']], cimis_long_values),axis=1)
    # convert values from object to numeric, coerce to allow nans
    cimis_clean['Value'] = pd.to_numeric(cimis_clean['Value'], errors='coerce')
    # clean up date
    cimis_clean.Date = pd.to_datetime(cimis_clean.Date)
    cimis_clean = cimis_clean.set_index('Date')
    return(cimis_clean)


# %%
def prep_cimis(cimis_clean):
    """ given a cimis dataframe convert to wide format to add to existing manual data downloads"""
    # wide format to save a csv to append to the other results
    cimis_out = cimis_clean.pivot_table(columns='variable', values='Value', index=['Date','Station']).reset_index('Station')
    cimis_out = cimis_out.rename(columns={'DayEto':'ETo (mm)','DayPrecip':'Precip (mm)'})
    cimis_out.Station = cimis_out.Station.astype(int)
    return cimis_out


# %%
cimis_clean = clean_cimis_json(cimis_json)
cimis_out = prep_cimis(cimis_clean)

# %%
cimis_all = pd.DataFrame()
# iterate over all cimis stations and append date
for s in cimis_stns['Stn Id'].values:
    sites = np.array([s]).astype(str)
    cimis_json = get_cimis_data(sites, strt_date, end_date, cimis_key)
    cimis_clean = clean_cimis_json(cimis_json)
    cimis_all = pd.concat((cimis_all, cimis_out))

# %%
cimis_all_out = cimis_all.rename(columns={'Station':'Stn Id'}).reset_index().merge(cimis_stns[['Stn Id', 'Stn Name']])
# cimis_all_out
cimis_all_out.to_csv(join(gwfm_dir, 'UZF_data',(strt_date+'_to_'+end_date).replace('-','_')), index=False)

# %% [markdown]
# # UZF update
# Currently new csv of the most recent CIMIS precipitation and ETo data have to be manually downloaded, then the code auto appends them together.   
# - `00a_DWR_LU_preparation.py` can be re-run as new LU maps are available, less urgent
# - `00b_Crop_ET.py` must be run to extend model into present, creates long_ETc_all_lu.hdf5
# - `00c_ETc_arrays.py` takes ETc and land use to make arrays of ETc for each day, creates irrigated_year.hdf5 and native_year.hdf5
# - `02_Basic soil budget - fields.py` runs the field by field soil water budget
#     - the script uses a start and end date to filter the input
# ## WEL
# Need to run `01_domestic_pumping_timeseries.py` to update the domestic pumping estimate using the updated ETc

# %% [markdown]
# # GHB Update
#
# All the scripts in the GHB folder should be run in numeric order.  
#
# - `00_Preliminary DWR GW data analysis.py` loads the DWR periodic measurements and crops to the buffered extent of the model domain. Then saves the points as shapefiles by season and year. The DWR periodic measurements file has to be manually downloaded currently, should update to scripted.  
# - `01_Kriging_for_GHB_boundary.py` takes the point shapefiles and applies kriging to each year/season combination to create rasters of groundwater elevation. It takes specific years overwhich to run for the contours since we don't want to have to always re-create past years.
#
# - `02_Distance GHB.py` and `03_GHB head correction.py` don't depend on the years as they run on the entire datasets available  
