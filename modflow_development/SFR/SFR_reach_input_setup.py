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

# %% [raw]
# Cosumnes Model 
# @author: Andrew

# %%
# standard python utilities
import os
from os.path import join, basename, dirname, exists
import sys
import glob
import pandas as pd
import numpy as np
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import pyproj # for converting proj4string
import shapely
import geopandas as gpd
from osgeo import gdal
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
flopy_dir = doc_dir+'/GitHub/flopy/'
if flopy_dir not in sys.path:
    sys.path.append(flopy_dir)
# sys.path
import flopy 

# from importlib import reload
# # importlib.reload
# reload(flopy)


# %%
from flopy.utils.geometry import Polygon, LineString, Point
# New model domain 52.9 deg
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')

# Need to check this when changing model domains
xul, yul = list(m_domain.geometry.values[0].exterior.coords)[1]
list(m_domain.geometry.values[0].exterior.coords)
# m_domain.geometry.values[0].exterior

# %%
sensors = gpd.read_file(gwfm_dir+"/DEM_data/allsensor_latlong.csv")
sensors.Latitude = sensors.Latitude.astype(np.float64)
sensors.Longitude = sensors.Longitude.astype(np.float64)
sensors.geometry = gpd.points_from_xy(sensors.Longitude, sensors.Latitude)
sensors.crs = 'epsg:4326'
sensors = sensors.to_crs('epsg:32610')

# %%
# import contextily as ctx
# fig, ax = plt.subplots(figsize = (10,10))
# mg.plot(ax=ax, alpha = 0.3)
# sensors.plot(ax =ax)

# ctx.add_basemap(ax, source = ctx.providers.Esri.WorldImagery, crs='epsg:26910', alpha = 0.6)
# # ctx.add_basemap(ax, source = ctx.providers.Esri.WorldStreetMap, crs=tg.crs.to_string())

# plt.xlabel('Easting (m)')
# plt.ylabel('Northing (m)')
# plt.savefig('Plots/Model_SFR_UZF_progress/Model_grid_overlay_satellite.png', dpi = 300)

# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
# grid_p = gpd.read_file(gwfm_dir+'/DIS_data/44_7_grid/44_7_grid.shp')


# Find Michigan Bar location
mb_gpd = sensors[sensors.Sensor_id == "MI_Bar"]
mb_grid = gpd.sjoin(mb_gpd, grid_p, how = 'left', predicate = 'intersects')

# %% [markdown]
# ### Raster files can be loaded using the `Raster.load` method

# %%
raster_name = gwfm_dir+"/DEM_data/USGS_ten_meter_dem/modeldomain_10m_transformed.tif"


# %%

# Based on Maribeth's grid aligned with Alisha's TPROGS model
# dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_52_9_200m_nearest.tsv', delimiter = '\t')
dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_52_9_200m_mean.tsv')

# import seaborn as sns
# sns.heatmap(dem_data, cmap = 'viridis', vmin = 0)

# %% [markdown]
# ## Import soil data for Lake Package, UZF Package, SFR Package hydraulic parameters

# %%
uzf_dir = gwfm_dir+'\\UZF_data'


# %%
grid_uzf = gpd.read_file(uzf_dir+'/final_grid_uzf/griduzf.shp')

# %%
soil_dir = join(uzf_dir, 'final_soil_arrays')
soilKs_array = np.loadtxt(soil_dir+'/final_soilKs.tsv', delimiter = '\t')
soiln_array = np.loadtxt(soil_dir+'/final_soiln.tsv', delimiter = '\t')
soileps_array = np.loadtxt(soil_dir+'/final_soileps.tsv', delimiter = '\t')


# %% [markdown]
# # SFR

# %%
sfr_dir = gwfm_dir+'/SFR_data/'

# %%

# Rivers and creeks in the larger area encompassing Cosumnes River in both South American and Cosumnes Subbasins
rivers = gpd.read_file(join(sfr_dir, "Sac_valley_rivers/Sac_valley_rivers.shp")

rivers = rivers.to_crs('EPSG:32610')
rivers_clip = gpd.clip(rivers, m_domain)
rivers_clip.plot()
# rivers_clip.GNIS_Name.unique()

# %%
print(rivers_clip.columns)
# Split into individual streams/creeks
cr_ind = rivers_clip[rivers_clip.GNIS_Name == 'Cosumnes River']
dc_ind = rivers_clip[rivers_clip.GNIS_Name == 'Deer Creek']
cc_ind = rivers_clip[rivers_clip.GNIS_Name ==  'Coyote Creek']
# Pull out data for each river/creek
cr = rivers_clip.loc[cr_ind.index,]
dc = rivers_clip.loc[dc_ind.index,]
cc = rivers_clip.loc[cc_ind.index,]


# %%
import fiona
from shapely.geometry import shape, mapping
from shapely.ops import linemerge

def line_2_points(cr, dline=10):
    # dline: how often to interpolate a point
    geom = linemerge(cr.geometry.unary_union)
    # # length of the LineString
    length = int(geom.length)
    point = np.zeros((int(length/dline)+1,3))
    for i, distance in enumerate(range(0, int(length), dline)):
             point[i,:] = geom.interpolate(distance).coords[:][0]
    point = point[:,[0,1]]
    return(point, geom)
point, geom = line_2_points(cr)
plt.plot(point[:,0],point[:,1])


# %%

# dem10 = rasterio.open(raster_name)
def sample_elev(raster_name, point):
    """ Sample a raster for an array of points
    raster_name: name of raster to be opened with rasterio
    point: n by 2 array of x,y coordinates
    """
    # create blank dataframe
    pnts = pd.DataFrame()
    # sample the raster and add a column to the dataframe
    with rasterio.open(raster_name) as src:
        pnts['z'] = [sample[0] for sample in src.sample(point)]
    # add the x and y coordinates to the dataframe
    pnts['easting'] = point[:,0]
    pnts['northing'] = point[:,1]
    # remove points that are null values (need to udpate)
    pnts = pnts[pnts.z > -1E4]
    return(pnts)
pnts_in = sample_elev(raster_name, point)



# %%
# pretty manually coded for the Cosumnes River
def manual_slope_fix(pnts, geom, dline=10):
    
    pnts['slope'] = 0.002
    for i in np.arange(1,len(pnts)):
        if pnts.z.values[i] >= pnts.z.values[i-1]:
            # if strtop is greater than previous strtop use previous elevation minus the average slope
            slope = ((np.max(pnts.z) - np.min(pnts.z))/geom.length)*dline
            if pnts.index[i] < 800:
                slope = 0.002
            elif pnts.index[i] < 2700:
                slope = 0.0003
            elif pnts.index[i] < 3200:
                slope = 0.001
            else:
                slope = 0.0003
            pnts.z.values[i] = pnts.z.values[i-1] - slope*dline
            pnts.slope.values[i] = slope
    return(pnts)
    
pnts = manual_slope_fix(pnts_in.copy(), geom)


# %%
from scipy.stats import gmean

# %%
# trying to incorporate a more fluid version
# will need to do quite a bit more work to make functional
# the hope was to use this for deer creek.
dline=10
delr = 200
def clean_profile(pnts, dline=10, delr=200, min_slope=1E-6, window=10):
    # find minimum value in XS related to thalweg
    pnts['z_m_min'] = pnts.z

    #roling mean of 6 window centered removes any negative slope
    pnts['z_m_min_cln'] = pnts.z_m_min.rolling(6,center=False).mean()
    
    # calculate slope and fill NAs, fill slope with nearby
    z_cln_diff = pnts.z_m_min_cln.diff().bfill()
    pnts['slope'] = z_cln_diff.abs()/dline
    # correct slope less than 1E-4
    pnts.loc[pnts.slope<min_slope,'slope'] = min_slope
    # rolling mean of slope to clean up slope for manning's
    pnts['slope_raw'] = pnts.slope.copy()
    # issue with the rolling fill is it overweights the high slope values
    # need to use a geometric mean
    pnts.slope = pnts.slope.rolling(window, center=True, min_periods=1).apply(gmean).bfill().ffill()
    
    # fix str bot so all is downward sloping
    for i in pnts.index[-2::-1]:
    # fill NAs due to rolling mean, with backward filling
        if np.isnan(pnts.loc[i,'z_m_min_cln']):
            pnts.loc[i,'z_m_min_cln'] = pnts.loc[i+1,'z_m_min_cln'] + pnts.loc[i,'slope']*delr
    
    for i in pnts.index[:-1]:
        if pnts.loc[i+1,'z_m_min_cln'] >= pnts.loc[i,'z_m_min_cln']:
            pnts.loc[i+1,'z_m_min_cln'] = pnts.loc[i,'z_m_min_cln'] - pnts.loc[i,'slope']*delr
    
    # calculate the elevation if we use slope only
    pnts['z_m_slope'] = pnts.z.max() - (dline*pnts.slope).cumsum()
    
    avg_slope = (pnts.z_m_min_cln.max() - pnts.z_m_min_cln.min())/(dline*len(pnts))
    # new column for easier modflow consistency
    pnts['z_m_final'] = pnts.z_m_slope.copy()
    return(pnts)
pnts_cln = clean_profile(pnts_in.copy(), window=500)

# %%
# pnts.slope.plot()
# pnts_cln.slope.plot()

# %%
fig,ax = plt.subplots(figsize = (6.5,3), dpi=300)

pnts_in.z.plot(ax=ax, label='Original')
pnts.z.plot(ax=ax, label='Manual')
pnts_cln.z_m_final.plot(ax=ax, label='Auto')
plt.legend()

# %%
pnts['Point_order'] = pnts.index

# %%
pnts_gpd = gpd.GeoDataFrame(pnts, geometry = gpd.points_from_xy(pnts.easting, pnts.northing))
pnts_gpd.crs = 'epsg:32610'

# %%
# # Samples the points every 10 meters to match with the 100 meter grid
grid_sfr = gpd.sjoin(grid_p, pnts_gpd, how = "inner", predicate= "intersects")


# # Set reach length for each reach based on the separation used to create the points from the line object
# # dline is 10 meters
grid_sfr['length_m'] = dline

# Dissolve the points again but using sum this time to get the total length of each reach
length_m = grid_sfr.dissolve(by = 'node', aggfunc = 'sum').length_m.values
# approximate the thalweg of each segment
min_z = grid_sfr.dissolve(by = 'node', aggfunc = 'min').z.values

# Dissolves the points every 10 meters to the 200 meter spacing, using mean because the interested component is elevation
grid_sfr = grid_sfr.dissolve(by = 'node', aggfunc = 'mean')
grid_sfr['length_m'] = length_m
grid_sfr['z_min'] = min_z

grid_sfr = grid_sfr.sort_values(by = 'Point_order')
grid_sfr['reach'] = np.arange(1,len(grid_sfr)+1)
grid_sfr['dist_m'] = grid_sfr.length_m.cumsum()

# only keep cells with segments greater than 100m long to reduce computation time
# and as short segments should contribute less seepage and it creates a more continuous line
grid_sfr_final = grid_sfr[grid_sfr.length_m>100]

# %%
fig,ax=plt.subplots(figsize=(6,4))
grid_sfr.plot(ax=ax, color='red')
grid_sfr_final.plot(ax=ax)

# %%
# grid_sfr.to_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')
# grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')


# %% [markdown]
# # Extras

# %%
# floodplain sensor info
fp_sensors = gpd.read_file(gwfm_dir+"/LAK_data/floodplain_logger_metadata.csv", header = True)
fp_sensors.Northing = fp_sensors.Northing.astype(np.float64)
fp_sensors.Easting = fp_sensors.Easting.astype(np.float64)
fp_sensors.geometry = gpd.points_from_xy(fp_sensors.Easting, fp_sensors.Northing)
fp_sensors.crs = 'epsg:32610'

# %%
# od_breach is the sensor location where the breach was made in the levees for flow to leave the river
od_breach = fp_sensors[fp_sensors['Logger Location']=='OD_Excavation']
od_breach


# %%
# first need to try with just adding extra cross section
# append the extra reach to the list of reaches and resort by reach
# grid_sfr = grid_sfr.append(grid_breach).sort_values(by = 'reach')
# # next, need to relabel reaches to account for the added reach witha duplicate number
# grid_sfr.reach = np.arange(1,len(grid_sfr)+1)
# grid_sfr

# %%
# # Buffer the location of the breach sensor to have overlap with the river streamline
# just sjoin the geometry because the extra info is unnecessary
# spatial join breach sensor polygon with sfr grid locations to find match
grid_breach = gpd.sjoin(grid_sfr, 
                        gpd.GeoDataFrame(geometry = od_breach.geometry.buffer(25), crs = 'epsg:32610'), how = "inner", op= "intersects")
# add a reach to the overlap cell that will be used to divert flow (there will be two reaches in one cell)
grid_breach

# %%
XS8pt = pd.read_csv(sfr_dir+'8pointXS.csv')
XSlocs = gpd.read_file(sfr_dir+'8pointXS_locs\\8pointXS_locs.shp')
XSlocs.crs = 32610

XSg  = gpd.sjoin(grid_sfr, XSlocs, how = "inner", op= "contains", lsuffix = 'sfr',rsuffix = 'xs')
# Append the grid_breach location to the list of cross sections to split the segment
XSg = XSg.append(grid_breach).sort_values('reach')
# Copy the XS site name from the previous last site to the breach site to keep same XS
XSg.Site.iloc[-1] = XSg.Site.iloc[-2]
len(XSg), len(XS8pt.loc[0,:])/2

# %% [markdown]
# ## Define segment data

# %% [markdown]
# Median grain size (mm) ranges from 1 mm – 30 mm along surveyed sites, which gives a range of 0.026-0.035 for a stable channel
# Moderate channel irregularity due to channel scouring and pools alternating, range of 0.006-0.010
# Gradual cross section change: 0.000 adjustment
# Effect of obstructions: minor due to occasional downed logs and debris in river, 0.005-0.015
# Amount of vegetation: large on banks due to willows and cottonwood trees, 0.025-0.050, and negligible in the channel
# Degree of meandering: minor due to levees, m = 1.0
#
# n = (nb+n1+n2+n3+n4)*m (b=base,1=surface irregularity, 2 = XS variation, 3 = obstructions, 4 = vegetation, m = correction for meandering)
# n = (0.03+0.08+0.01) = 0.048 in channel
# n = (0.048 +0.03) = 0.078 on banks
#
