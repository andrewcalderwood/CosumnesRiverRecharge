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
#
# Newest version of script to identify the minimum reach elevations in a more fluid way. This works for both the Cosumnes River and Deer Creek.
# This script supercedes SFR_reach_input_setup.py which used a simpler elevation sampling directly along the NHD lines.

# Originally developed under CosumnesRiverRecharge/modflow_development/SFR/SFR_XS_to_reach.py


# %%
# standard python utilities
import os
from os.path import join, basename, dirname, exists
import sys
import pandas as pd
import numpy as np
import time


# standard geospatial python utilities
import geopandas as gpd
import rasterio

# %%
import fiona
from shapely.geometry import shape, mapping
from shapely.ops import linemerge

from shapely.geometry import LineString

# %% [markdown]
# # Cross-sections to find minimum
# Better approach than directly sampling the NHD lines because using a cross-section perpendicular to a cell should identify the true channel in a cell.

# %%
def make_transects(geom, dline, xs_width):
    """ Create regular transects along a line
    INPUT:
    geom: shapely geometry object, usually created from union_all of a geodataframe geometry
    dline: how often to interpolate a point along a line, typically should be no less than the raster
        resolution that will be used for sampling (e.g., 10m for 10m USGS DEM)
    xs_width: cross-section width (units of crs) that should be used in total (1/2 each river side)
        should extend beyond actual desired width to allow for some uncertainty in where the thalweg will be
    OUTPUT:
    transg: geodataframe with transects interpolated at distance dline with width of xs_width (no crs)
    """
    # dline: how often to interpolate a point along a line
    # xs_width: 300 meter width so that 100 meter can be cutoff to have 200m around true thalweg
    long_dline = np.copy(dline)
    # # length of the LineString
    length = int(geom.length)
    
    transects = pd.DataFrame(np.zeros((int(np.ceil(length/dline)),1)), columns = ['line'])
    transects['geometry'] = LineString([(0,0),(0,1)]) #initiate LineString geometry column
    
    for i, distance in enumerate(np.arange(0, int(length), dline)):
        short_line = LineString([geom.interpolate(distance),geom.interpolate(distance+dline)])
        geom_left = short_line.parallel_offset(xs_width/2,'left', resolution = 32, join_style = 2)
        geom_right = short_line.parallel_offset(xs_width/2,'right', resolution = 32, join_style = 2)
        # old method was geom_left.boundary but broke with update
        perp_line = LineString([list(geom_left.coords)[0], list(geom_right.coords)[0]])
        transects.loc[i,'geometry'] = perp_line
    # convert into geodataframe
    transg = gpd.GeoDataFrame(transects)
    transg['line'] = np.arange(0,len(transg))
    return(transg)


def create_points(transg, xs_width, dline=10):
    """
    Given transect lines, split them up into data points for sampling a raster
    INPUT:
    transg: geodataframe with transects interpolated at a regular interval with width of xs_width (no crs)
    OUTPUT:
    xs_all: geodataframe with the transects converted into points with spacing of dline along each transect
    """
    xs_all = gpd.GeoDataFrame(pd.DataFrame(columns=['xs_num','dist_from_right_m','geometry']))
    xs = gpd.GeoDataFrame(pd.DataFrame(np.zeros((int(xs_width/dline),2)), columns=['xs_num','dist_from_right_m']))
    
    for j in np.arange(0,len(transg)):
        xs['geometry'] = Point([(0,0)])
        xs['xs_num'] = j
        
        # pick one geometry at a time
        geom = transg.iloc[j].geometry
    
        # # # length of the LineString
        length = int(geom.length)
        # create discrete points for each lien
        for i, distance in enumerate(range(0, int(length), dline)):
            point = geom.interpolate(distance)
            xs.loc[i,'geometry'] = point
            # xs.loc[i,'dist_from_right_m'] = i
            xs.loc[i,'dist_from_right_m'] = distance
        # append individual cross section to all dataframe
        xs_all = pd.concat((xs_all, xs))
    return(xs_all)



# %%
def sample_points(xs_all, crs, raster_name):
    """
    given a geodatafarme with points sample the elevation from a raster
    """
    xs_all.crs=crs
    xs_all['Easting'] = xs_all.geometry.x
    xs_all['Northing'] = xs_all.geometry.y
    point = xs_all.loc[:,['Easting','Northing']].values
    
    with rasterio.open(raster_name) as src:
        xs_all['z_m'] = [sample[0] for sample in src.sample(point)]
    
    with rasterio.open(raster_name) as dem_f:
        dem_nodata = dem_f.nodata
    
    # remove any NA values picked up from DEM raster
    xs_all.loc[xs_all['z_m']==dem_nodata,['z_m']] = np.nan
    xs_all.index = np.arange(0,len(xs_all))
    return(xs_all)





# %%
# trying to incorporate a more fluid version
# will need to do quite a bit more work to make functional
# the hope was to use this for deer creek.
def clean_profile(pnts, z_col, dline=10, min_slope=1E-6, max_slope=1E-3, window=10):
    # find minimum value in XS related to thalweg
    pnts['z_m_min'] = pnts[z_col]

    #roling mean of 6 window centered removes any negative slope
    pnts['z_m_min_cln'] = pnts.z_m_min.rolling(6,center=False).mean()
    pnts['z_m_min_mean'] = pnts.z_m_min_cln.copy()
    
    # calculate slope and fill NAs, fill slope with nearby
    z_cln_diff = pnts.z_m_min_cln.diff().bfill()
    pnts['slope'] = z_cln_diff.abs()/dline
    # correct slope less than 1E-4 (flopy suggested minimum)
    pnts.loc[pnts.slope<min_slope,'slope'] = min_slope
    pnts.loc[pnts.slope>max_slope,'slope'] = max_slope
    # rolling mean of slope to clean up slope for manning's
    pnts['slope_raw'] = pnts.slope.copy()
    # issue with the rolling fill is it overweights the high slope values
    # need to use a geometric mean
    pnts.slope = pnts.slope.rolling(window, center=True, min_periods=1).apply(hmean).bfill().ffill()
    # pnts.slope = pnts.slope.rolling(window, center=True, min_periods=1).apply(np.mean).bfill().ffill()
    
    # fix str bot so all is downward sloping
    for i in pnts.index[-2::-1]:
    # fill NAs due to rolling mean, with backward filling
        if np.isnan(pnts.loc[i,'z_m_min_cln']):
            pnts.loc[i,'z_m_min_cln'] = pnts.loc[i+1,'z_m_min_cln'] + pnts.loc[i,'slope']*dline
    
    for i in pnts.index[:-1]:
        # if the elevations diverge well below existing channel then decrease slope to minimum
        if pnts.loc[i,'z_m_min_cln'] < pnts.loc[i,'z_m_min'] - 0.5:
            pnts.loc[i, 'slope'] = np.min((min_slope, pnts.loc[i, 'slope']))
        # correct elevation to ensure always downsloping
        if pnts.loc[i+1,'z_m_min_cln'] >= pnts.loc[i,'z_m_min_cln']:
            pnts.loc[i+1,'z_m_min_cln'] = pnts.loc[i,'z_m_min_cln'] - pnts.loc[i,'slope']*dline
        # correct down slope to avoid extreme drop-offs
        if (pnts.loc[i,'z_m_min_cln'] - pnts.loc[i+1,'z_m_min_cln'])/dline > max_slope:
            pnts.loc[i+1,'z_m_min_cln'] = pnts.loc[i,'z_m_min_cln'] - max_slope*dline
    
    # calculate the elevation if we use slope only
    pnts['z_m_slope'] = pnts[z_col].max() - (dline*pnts.slope).cumsum()
    
    avg_slope = (pnts.z_m_min_cln.max() - pnts.z_m_min_cln.min())/(dline*len(pnts))
    # new column for easier modflow consistency
    pnts['z_m_final'] = pnts.z_m_slope.copy()
    return(pnts)


# %%
# sample elevations from 10 m raster
def sample_raster(raster_name, river_pts):
    with rasterio.open(raster_name) as rio:
        # sample missing elevation values
        rio_sample = river_pts.to_crs(rio.crs)
        # sample elevations from the raster
        z = rio.sample(list(zip(rio_sample.geometry.x, rio_sample.geometry.y)))
        # extract elevations from generator object
        river_pts['z'] = [sample[0] for sample in z]
    # calculate average elevations by cell_id
    avg_reach_elev = river_pts.groupby('cell_id')['z'].mean().reset_index()
    return(avg_reach_elev)


