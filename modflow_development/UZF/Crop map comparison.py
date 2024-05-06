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
# Compare the *LandIQ* and *USDA CDL* data sets for refinement

# %%
# standard python utilities
import os
import sys
import glob
import pandas as pd
import numpy as np

# standard geospatial python utilities
# import pyproj # for converting proj4string
import shapely
import geopandas as gpd
import fiona
import rasterio

import flopy

# %%
from os.path import basename, dirname, exists, join
## Set up directory referencing
# Package data
usr_dir = os.getcwd()
while basename(usr_dir) !='Users':
    temp = basename(usr_dir)
    usr_dir = dirname(usr_dir)
usr_dir += '/'+temp
box_dir = usr_dir+'/Box/'
gwfm_dir = usr_dir+'/Box/research_cosumnes/GWFlowModel/'


dis_dir = gwfm_dir+'DIS_data/'
uzf_dir = gwfm_dir+'UZF_data/'

print(gwfm_dir)

# %%
# New model domain 52.9 deg
m_domain = gpd.read_file(dis_dir+'NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')

xll, yll = list(m_domain.geometry.values[0].exterior.coords)[0]
#Maribeth's model parameters, had to switch nrow and ncol due to her issue in xul, yul
nrow=100
ncol=230
delr=np.repeat(200,ncol)
delc=np.repeat(200,nrow)
rotation=52.9
modelgrid = flopy.discretization.StructuredGrid(xoff=xll, yoff=yll, proj4='EPSG:32610', angrot=rotation,
                                   delr=delr, delc=delc, nrow=nrow,ncol=ncol)

# %%
proj_dir = box_dir+'/SESYNC_Paper1'
# parcel data for 2022 from Yusuke (Sacramento, San Joaquin), already clipped to domain
parcels = gpd.read_file(proj_dir+'/Parcels - All counties/parcels_all_counties_model.shp')
# change crs to utm zone 10
parcels = parcels.to_crs('epsg:32610')

# %%
# burn river shapefile into the 10 meter dem and then read it out to find the cells where it is
# Full size dem of northern sac valley
# raster_name = gwfm_dir+"/DEM_data/USGS_ten_meter_dem/modeldomain_10m_transformed.tif"

# dem = rasterio.open(raster_name)
# dem_10m = dem.read((1,))[0,:,:]

# # affine = dem.affine # didn't work
# affine = dem.meta['transform']

# %%
# runs quickly
# from rasterstats import gen_zonal_stats
# zs_gen = gen_zonal_stats(parcels, raster=raster_name, stats=['min', 'max', 'mean', 'median', 'majority','std'])

# %%
# takes several minutes
# zs_parcels = zonal_stats(parcels, raster=raster_name, stats=['min', 'max', 'mean', 'median', 'majority','std'])
# # zonal_stats(parcels, dem_10m, affine=affine, stats=['min', 'max', 'mean', 'median', 'majority', 'std'], nodata=-999)
# # convert to dataframe
# zs_df = pd.DataFrame(zs_parcels)
# # join zone stats of DEM to parcel data
# zs_df = parcels.join(zs_df)
# # save to shapefile
# zs_df.to_file(proj_dir+'/parcel_zonalstats/parcel_elevation_m_statistics.shp')

# %%


# crop_path = uzf_dir+'Newmodeldomain_cropdata_2007_2019'
crop_path = uzf_dir+'Modeldomain_cropdata_2007_2021'

crop_raster_names = glob.glob(crop_path+'/*.tif')
crop_dbf_names = glob.glob(crop_path+'/*.dbf')

import pathlib

crop_raster_list = list(pathlib.Path(crop_path).glob('*.tif'))
crop_dbf_list = list(pathlib.Path(crop_path).glob('*.dbf'))


# %%
# landiq data
landIQ_fn = glob.glob(join(uzf_dir, 'LandIQ/*/*.shp'))
landIQ_dbf = glob.glob(join(uzf_dir, 'LandIQ/*/*.dbf'))

n=-1
# takes a while to load
land_iq = gpd.read_file(landIQ_fn[n])
land_dbf = gpd.read_file(landIQ_dbf[n])

# %%
# crop to model domain, very slow as well
land_iq_m = gpd.overlay(land_iq.to_crs(m_domain.crs), m_domain)

# %%
land_iq_m.CLASS1.unique()
land_iq_m.SUBCLASS1.unique()

# %%
land_dbf.columns

# %%
# useful code but ended up manually summarizing the classes and subclasses
# landIQ_xml = glob.glob(join(uzf_dir, 'LandIQ/*/lyr_converted.xml'))
# pd.read_xml? #namespaces, iterparse
# pd.read_xml(landIQ_xml[0], namespaces={'symbols':'symbol'}) 
# xpath seems to specify the level of interest
# xml_names = pd.read_xml(landIQ_xml[0], xpath=".//symbol")
# xml_names = xml_names[['name']]
# xml_names[['class','class_name']] = xml_names.name.str.split('|',expand=True)
# xml_names

# %%
# A-B, A-C, A-D are the percent of the season for each Kc
# The dates are the dates of the growing season
Kc = pd.read_csv(uzf_dir+'Kc/Kc_Current.csv',skiprows = 1)
Kc = Kc.rename(columns={'Unnamed: 0' : 'Match_kc'})


# %%
def crop_raster2array(raster_file_name, dbf_filename, Kc):
    # read in crop raster and database of crop values
    src = rasterio.open(raster_file_name)
    crop_dbf = pd.DataFrame(gpd.read_file(dbf_filename))
        
    band1 = src.read(1)
    band1 = band1.astype('int16')
    band1 = band1.astype(float)
    band3 = band1.reshape(1, band1.shape[0], band1.shape[1])[:,:,:]
    # turn the rasterio format into a flopy format
    # for some reason flopy will not directly read in the raster .tif
    croprio = Raster(band3, bands = (1,), crs = src.crs, transform = src.transform, 
                 nodataval = 255)
    # no longer need to crop as data from USDA was downloaded to fit domain
    # crop the raster to the model domain
#     croprio.crop(vertices, invert=False)
    
    # The original crop raster has a cell size of 56 by 56 m so if there are less than 4 cells of one crop
    # then for certain they will not fill one cell and most likely have minimal impact considering there are 
    # 6300 model cells in one layer
    crop_hist = np.histogram(band3, bins = np.arange(0,257))
    # only need to filter out those that have no cells in the domain
    crops_in_domain = crop_hist[1][:-1][crop_hist[0]>1]
#     crops_in_domain = crop_hist[1][:-1]
    
    domain_dbf = crop_dbf.iloc[crops_in_domain]

    domain_dbf['CLASS_NAME'] = domain_dbf.CLASS_NAME.str.replace('Dbl Crop ','')
    domain_dbf['crop_hist'] = crop_hist[0][crops_in_domain]
    # remove the 0/background value from the domain_dbf because it messes with the crop histogram
    domain_dbf = domain_dbf.drop(0)
    # then create a column with the percent of the crop
    domain_dbf['crop_percent'] = 100*domain_dbf.crop_hist/domain_dbf.crop_hist.sum()
    Kcmatch = pd.read_csv(uzf_dir+'Kc/Cosumnes_crops.csv', index_col = 0)

    # domain_dbf['crop1'] = domain_dbf.CLASS_NAME.str.split('/', expand = True)[0]
    # domain_dbf['crop2'] = domain_dbf.CLASS_NAME.str.split('/', expand = True)[1]

    domain_dbf = domain_dbf.merge(Kcmatch, on = 'CLASS_NAME')
    domain_dbf = domain_dbf.merge(Kc, left_on = 'Match_kc', right_on = 'Match_kc', how = 'left')
    return(croprio, domain_dbf)



# %%
def calc_kc_dates(year, domain_dbf):
    # The year for each crop for each set of dates needs to change iteratively for each crop individually because
    # some crops have dates that extend into the next year that must not change until the final date of the 
    # season is reached (e.g. 2018-11-01 to 2019-09-17 must stay 2018 and 2019 until 2019-09-17 is reached)
#     i = 2018
    dates = domain_dbf.loc[:,['Beg Month','Beg Day', 'End Month', 'End Day', 'A-B', 'A-C', 'A-D']]

    # Set the pandas datetime from the start and end dates of crops
    # need to just takes .values or indexing will be wrong and mismatch dates to rows
    dates['A'] = pd.to_datetime({'year': year, 'month':dates['Beg Month'].values, 'day': dates['Beg Day'].values}).values
    dates['E'] = pd.to_datetime({'year': year, 'month':dates['End Month'].values, 'day': dates['End Day'].values}).values
    # Make correction for any end dates that are in the next year
    dates.E.loc[dates.E < dates.A] = dates.E.loc[dates.E < dates.A] + pd.offsets.DateOffset(years=1)

    # Get the length of the growing periods
    dates['num_days'] = dates.E-dates.A

    # set the end date of growing period A/ start of period B
    dates['B'] = dates.A + dates.num_days*(dates['A-B']/100)
    # Round the dates, as we will be one a daily time step
    dates.B = pd.to_datetime(dates.B.dt.date)

    # # set the end date of growing period B/ start of period C
    dates['C'] = dates.B + dates.num_days*((dates['A-C']-dates['A-B'])/100)

    # # set the end date of growing period C/ start of period D
    dates['D'] = dates.C + dates.num_days*((dates['A-D']-dates['A-C'])/100)
    return(dates)


# %%
def ETc_calc(ET_final, dates, domain_dbf):
    ETc = np.zeros((yearlen,nrow,ncol))
    Kc_arr = np.zeros((yearlen,nrow,ncol))

    time = 0
    for dt in pd.date_range(yr_strt, yr_end):
        # First step is to get the current Kc for each crop for the time step
        domain_dbf.Kc.loc[dt > dates.A] = domain_dbf.loc[dt > dates.A, 'Kc1']
        domain_dbf.Kc.loc[dt > dates.B] = domain_dbf.loc[dt > dates.B, 'Kc2']
        domain_dbf.Kc.loc[dt > dates.C] = domain_dbf.loc[dt > dates.C, 'Kc3']
        domain_dbf.Kc.loc[dt > dates.D] = domain_dbf.loc[dt > dates.D, 'Kc4']
        domain_dbf.Kc.loc[dt > dates.E] = domain_dbf.loc[dt > dates.E, 'Kc4']
        for i,j in zip(domain_dbf.index.values, domain_dbf.Kc.values):
            Kc_arr[time,crop_data==i] = j
        ETc[time,:,:] = Kc_arr[time,:,:]*ET_final[time,:,:]
        time += 1
    return(ETc, Kc_arr)


# %%
from flopy.utils import Raster


# %%

# %%
# iterate over each year and calculate Kc with crop coefficients
# for y in np.arange(pd.to_datetime(strt_date).year, pd.to_datetime(end_date).year+1):
for y in [2019]:
    fn_r = glob.glob(crop_path+'/CDL_'+str(y)+'*.tif')[0]
    fn_d = glob.glob(crop_path+'/CDL_'+str(y)+'*.dbf')[0]

#     out_fn = uzf_dir+'ETa_all_txt_arrays/ETa_array_'+str(y)+'.tsv'
#     # only create new arrays if they don't exist
#     if exists(out_fn):
#         print('Already exists')
#     elif not exists(out_fn):
#         # set start and end date for range for the year to be iterated over
#         yr_strt = pd.to_datetime(str(y)+'-01-01')
#         yr_end = pd.to_datetime(str(y)+'-12-31')
#         if yr_strt < pd.to_datetime(strt_date):
#             yr_strt = pd.to_datetime(strt_date)
#         if yr_end > pd.to_datetime(end_date):
#             yr_end = pd.to_datetime(end_date)
    if y>0:

        # for each year, import the new crop raster and resample to the model grid
        # and filter out the database of crops to match those in the domain
        croprio, domain_dbf = crop_raster2array(fn_r, fn_d, Kc)
#         file_num +=1

#         crop_data = croprio.resample_to_grid(modelgrid,
#                                     band=croprio.bands[0], method="nearest")
        crop_data = croprio.resample_to_grid(modelgrid,
                                    band=croprio.bands[0], method="nearest")
        # adjust domain_dbf to account for resampling
        resampled_crops = np.append(np.unique(crop_data).astype(int),np.unique(crop_data)[-1])
        resampled_hist =  np.histogram(crop_data, bins = resampled_crops)
        # convert histogram to dataframe to join with domain database info for crops
        resampled_hist = pd.DataFrame(np.transpose(np.vstack((resampled_hist[0], resampled_hist[1][:-1]))))
        resampled_hist.columns =  ['crop_hist','VALUE']
        resampled_hist.VALUE = resampled_hist.VALUE.astype(int)

        resampled_df = resampled_hist.set_index('VALUE').join(domain_dbf.drop('crop_hist', axis=1).set_index('VALUE') , on = 'VALUE', how = 'inner')
        resampled_df.crop_percent = 100*resampled_df.crop_hist/resampled_df.crop_hist.sum()
        resampled_df['Kc'] = 0


# %%
resampled_df

# %%
