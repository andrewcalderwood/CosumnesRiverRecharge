# ---
# jupyter:
#   jupytext:
#     formats: py:percent
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
# # Identify the spatial extent of cooarse pathways
# The Cosumnes model is based only on geologic data in the Cosumnes River area and doesn't include the incised valley fill (IVF) from the American River fan. The AEM data shows there is a strong appearnce of the IVF below the middle lower Cosumnes River. This means it should be included in the geologic model, but rather than produce a TProGS model of the American and crop out the channel deposits we can use the AEM data combined with the cross-sections and maps from Meirovitz (2010) to identify the true channel extents.

# %%
# standard python utilities
import os
import sys
from os.path import basename, dirname, join, exists, expanduser
import glob
import time

import pandas as pd
import numpy as np
import numpy.ma as ma
from scipy.stats import gmean

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
# import pyproj # for converting proj4string
import geopandas as gpd
from osgeo import gdal
import rasterio
from rasterstats import zonal_stats

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study


sfr_dir = gwfm_dir+'/SFR_data/'
upw_dir = gwfm_dir+'/UPW_data/'

proj_dir = join(gwfm_dir, 'Projects','EcoFIP')
gis_dir = join(proj_dir,'GIS')


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

from importlib import reload

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

import map_cln
reload(map_cln)

# from mf_utility import get_layer_from_elev
from report_cln import base_round
from map_cln import plt_cln

# %% [markdown]
# # Load data

# %%
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')

# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_elev = gpd.read_file(join(gwfm_dir,'DIS_data','grid_elevation_m_statistics.shp'))

# %%
grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')

# %%
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
rivers = gpd.read_file(join(sfr_dir,"Sac_valley_rivers/Sac_valley_rivers.shp"))
rivers = rivers.to_crs('EPSG:32610')

rivers_clip = gpd.clip(rivers, m_domain)
rivers_clip = rivers_clip[['GNIS_ID','GNIS_Name','LengthKM', 'geometry']]
cr = rivers_clip[rivers_clip.GNIS_Name=='Cosumnes River']
mr = rivers_clip[rivers_clip.GNIS_Name=='Mokelumne River']


# %% [markdown]
# ## ecofip reference
# Shapefiles of EcoFIP parcels for analysis on the 10-acre grid scale and parcel scale.
# - simplify to ID columns and sample elevation statistics for each parcel
# - save geodataframes as shapefiles for future reference

# %%
# ecofip parcel scale analysis
ecofip_grid_name = 'accumAvg_Hist_regionRanks-g0'
sq_grid = gpd.read_file(join(gis_dir, ecofip_grid_name+'_10AcGrid',ecofip_grid_name+'.shp'))
parcel_grid = gpd.read_file(join(gis_dir, ecofip_grid_name+'_Parcels',ecofip_grid_name+'.shp'))
# FullArea and AnlysArea are in acres, 
# FullArea is likely parcel full area and AnlysArea is area after cropping for river/flooding
grid_cols = ['Region','FullArea','AnlysArea','RchIds','geometry']
# simplify to parcel ids only
sq_id = sq_grid[grid_cols]
parcel_id = parcel_grid[grid_cols]


# %%
def elev_stats(raster_name, gdf):
    """ apply raster stats and join to geodataframe """
    # takes several minutes
    zs_parcels = zonal_stats(gdf, raster=raster_name, stats=['min', 'max', 'mean', 'std'])
    # convert to dataframe
    zs_df = pd.DataFrame(zs_parcels)
    # join zone stats of DEM to parcel data
    zs_df = gdf.join(zs_df)
    return zs_df
    


# %%
dem_raster_name = gwfm_dir+"/DEM_data/USGS_ten_meter_dem/modeldomain_10m_transformed.tif"
parcel_elev = elev_stats(dem_raster_name, parcel_id.to_crs(raster_crs))
parcel_elev.to_file(join(gis_dir, 'analysis_unit_reference', 'parcels_elevation.shp'))

sq_elev = elev_stats(dem_raster_name, sq_id.to_crs(raster_crs))
sq_elev.to_file(join(gis_dir, 'analysis_unit_reference', 'sq_10ac_elevation.shp'))


# %% [markdown]
# ## load interpolated groundwater data


# %%
ghb_dir = join(gwfm_dir, 'GHB_data')
# elevations are in feet
ghb_fp = join(ghb_dir, 'interpolated_data')

ghb_files =  os.listdir(ghb_fp)
y=2014
raster_name = join(ghb_fp, 'fall'+str(y)+'_kriged.tif')
with rasterio.open(raster_name) as r:
    # print(r.meta)
    raster_crs = r.crs


# %%
def gdf_sample_gwe(gdf_elev, gdf_id, grid_cols, season='fall'):
    """ Given a geodataframe of the EcoFIP analysis grid return the groundwater elevation
    statistics for each year specified"""
    gdf_gw = gdf_elev.copy()
    gdf_gw = gdf_gw[grid_cols+['mean']].rename(columns={'mean':'gse_m'})
    # sampling is pretty quick for an individual year
    # sample the groundwater elevation to the parcel level 
    gwe_yrs = np.arange(2000,2021)
    for y in gwe_yrs:
        raster_name = join(ghb_fp, season+str(y)+'_kriged.tif')
        with rasterio.open(raster_name) as r:
            raster_crs = r.crs
        gw_gdf = elev_stats(raster_name, gdf_id.to_crs(raster_crs))
        gdf_gw['gwe'+str(y)] =  gw_gdf['mean'] * 0.3048
    # return sampled groundwater elevations for saving
    return gdf_gw


# %%
def clean_gdf_gw(gdf_gw, grid_cols, grid_id, season):
    """ Clean up sampled GWE data and ave to shapefile"""
    # convert data to long format for easier data use
    gdf_gw_long = gdf_gw.melt(id_vars = grid_cols+['gse_m'], value_name='gwe_m')
    # identify groundwater elevation by year and season sampled
    gdf_gw_long['year'] = gdf_gw_long.variable.str.extract(r'(\d{4})').astype(int)
    gdf_gw_long['season'] = 'fall'
    # calculate depth to groundwater
    gdf_gw_long['dtw_m'] = gdf_gw_long.gse_m - gdf_gw_long.gwe_m
    # save long format (110 MB)
    gdf_gw_long.to_file(join(gis_dir, 'analysis_unit_reference', 
                             'GW_elevations_long_'+grid_id+'_'+season+'.shp'))
    return None


# %% [markdown]
# ## fall groundwater levels

# %%
season='fall'
parcel_gw = gdf_sample_gwe(parcel_elev, parcel_id, grid_cols, season)
clean_gdf_gw(parcel_gw, grid_cols, 'parcel', season)


# %%
sq_gw = gdf_sample_gwe(sq_elev, sq_id, grid_cols, season)
clean_gdf_gw(sq_gw, grid_cols, 'sq', season)


# %% [markdown]
# ## spring groundwater level

# %%
season='spring'
parcel_gw = gdf_sample_gwe(parcel_elev, parcel_id, grid_cols, season)
clean_gdf_gw(parcel_gw, grid_cols, 'parcel', season)


# %%
sq_gw = gdf_sample_gwe(sq_elev, sq_id, grid_cols, season)
clean_gdf_gw(sq_gw, grid_cols, 'sq', season)


# %%
