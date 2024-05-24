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
# Cosumnes Model 
# @author: Andrew
#
# Newest version of script to identify the minimum reach elevations in a more fluid way. This works for both the Cosumnes River and Deer Creek.
# This script supercedes SFR_reach_input_setup.py which used a simpler elevation sampling directly along the NHD lines.

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
rivers = gpd.read_file(join(sfr_dir, "Sac_valley_rivers/Sac_valley_rivers.shp"))

rivers = rivers.to_crs('EPSG:32610')
rivers_clip = gpd.clip(rivers, m_domain)
rivers_clip.plot()
# rivers_clip.GNIS_Name.unique()

# %% [markdown]
# Need to update this script to take in label as an argument likely via a class to more easily produce the outputs

# %%
# print(rivers_clip.columns)
# Split into individual streams/creeks
cr_ind = rivers_clip[rivers_clip.GNIS_Name == 'Cosumnes River']
dc_ind = rivers_clip[rivers_clip.GNIS_Name == 'Deer Creek']
cc_ind = rivers_clip[rivers_clip.GNIS_Name ==  'Coyote Creek']
# Pull out data for each river/creek
cr = rivers_clip.loc[cr_ind.index,]
dc = rivers_clip.loc[dc_ind.index,]
cc = rivers_clip.loc[cc_ind.index,]

# geom = linemerge(cr.geometry.unary_union)
label = 'Cosumnes River'
xs_width = 300

label = 'Deer Creek'
xs_width = 100

r_ind = rivers_clip[rivers_clip.GNIS_Name == label]
geom = linemerge(r_ind.geometry.unary_union)
print(label)



# %% [markdown]
# # Cross-sections to find minimum
# Better approach than directly sampling the NHD lines because using a cross-section perpendicular to a cell should identify the true channel in a cell.

# %%
import fiona
from shapely.geometry import shape, mapping
from shapely.ops import linemerge

from shapely.geometry import LineString


# %%
def make_transects(geom, dline, xs_width):
    """ Create regular transects along a line"""
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


# %%
long_dline = 20
transg = make_transects(geom, dline=long_dline, xs_width = xs_width)

# %%
# transg.plot()
# example zoom in
# transg.loc[1800:2000].plot()

# %%
# # takes longer to run with points every 20 m steps
# # how often to interpolate a point in a XS
dline = 10

def create_points(transg, xs_width, dline=10):
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
xs_all = create_points(transg, xs_width=xs_width, dline=10)


# %%
def sample_points(xs_all):
    xs_all.crs='epsg:32610'
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
xs_pts = sample_points(xs_all.copy())

# %%
## center cross-sections on the thalweg
# need to drop XS where there are a lot of NAs 
max_xs_npts = int(xs_width/dline) # maximum points per XS
xs_num_full = xs_pts.xs_num.unique()[(xs_pts.groupby('xs_num').count()['z_m']>int(max_xs_npts*0.75)).values]
xs_pts = xs_pts[xs_pts.xs_num.isin(xs_num_full)]

xs_pts_width = int(0.5*xs_width/dline)
# clean xs to only center 200 m out from true thalweg rather than NHD line
xs_all_cln = gpd.GeoDataFrame()
# correct XS points by those that were chosen for final XS
for n in xs_pts.xs_num.unique():
    # mid = xs_pts.loc[xs_pts.xs_num==n,'z_m'].iloc[int(xs_width/2)-50:int(xs_width/2)+50].idxmin()
    mid_df = xs_pts.loc[xs_pts.xs_num==n,'z_m']
    mid_center = int(np.median(mid_df.index.values))
    mid = mid_df.loc[np.max((mid_center-xs_pts_width,0)):mid_center+xs_pts_width].idxmin()
    mid_vals = xs_pts.loc[np.max((mid- 2*xs_pts_width,0)):mid+ 2*xs_pts_width].copy()
    # correct distance from right
    mid_vals['dist_from_right_m'] -= mid_vals['dist_from_right_m'].min()
    xs_all_cln = pd.concat((xs_all_cln, mid_vals))


# %%
# convert XS to columns
xs_all_df = xs_all_cln.pivot_table(columns='xs_num',values='z_m',index='dist_from_right_m')
xs_all_df.iloc[:,1000:1100:10].plot()

# %%
# get minimum elevation for each cross-section (thalweg)
xs_num_min = xs_all_cln.groupby('xs_num')[['z_m']].min()
# xs_num_min.plot()

# %%
# join minimum values with geodataframe
xs_mins = xs_all_cln.merge(xs_num_min.reset_index(), how='inner')
# when multiple points have the same minimum drop the duplicates ( should be adjacent)
xs_mins = xs_mins.drop_duplicates('xs_num').reset_index(drop=True)


# %%
xs_mins.z_m.plot()

# %% [markdown]
# # points to reach

# %%
from scipy.stats import gmean, hmean


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
# pnts_cln = clean_profile(pnts_in.copy(), window=500)


# %%
# minimum slope from flopy SFR is 1E-4
xs_mins_cln = clean_profile(xs_mins, 'z_m', dline=long_dline, window=50, min_slope=1E-4, max_slope=1E-2)

# %%
xs_mins_cln.slope.describe()

# %%
# xs_mins_cln.plot(y='z_m', kind='line')
xs_mins_cln[['z_m', 'z_m_min_cln']].plot()
# xs_mins_cln['slope'].plot()

# %%
# clean XS before saving output
xs_mins_out = xs_mins_cln[['xs_num','z_m','z_m_min_cln','slope','geometry']].copy()

# %%
# # Samples the points every 10 meters to match with the 100 meter grid
grid_sfr = gpd.sjoin(grid_p, xs_mins_out, how = "inner", predicate= "intersects")

# # Set reach length for each reach based on the separation used to create the points from the line object
# # dline is 10 meters
grid_sfr['length_m'] = long_dline

# Dissolve the points again but using sum this time to get the total length of each reach
length_m = grid_sfr.dissolve(by = 'node', aggfunc = 'sum').length_m.values
# approximate the midpoint elevation of each cell
mean_z = grid_sfr.dissolve(by = 'node', aggfunc = 'mean').z_m_min_cln.values

# Dissolves the points every 10 meters to the 200 meter spacing, using mean because the interested component is elevation
grid_sfr = grid_sfr.dissolve(by = 'node', aggfunc = 'mean')
grid_sfr['length_m'] = length_m
grid_sfr['z_min'] = mean_z

grid_sfr = grid_sfr.sort_values(by = 'xs_num')
grid_sfr['reach'] = np.arange(1,len(grid_sfr)+1)
grid_sfr['dist_m'] = grid_sfr.length_m.cumsum()


# %%
# # code to confirm that a minimum length doesn't disconnect reaches
min_length = 10
length_step = 10
diff_bool = np.array([True])
while all(diff_bool==True):
    grid_sfr_temp = grid_sfr[grid_sfr.length_m>min_length].copy()
    # determine the difference between the current cell and next
    cell_diff = grid_sfr_temp[['row','column']].diff().abs().iloc[1:]
    # should not allow more than a 1 cell difference in either
    diff_bool = (cell_diff.row<2)&(cell_diff.column<2)
    # increase step in minimum length
    min_length += length_step
# final min length is from previous step
min_length -= length_step*2
print('Final min length is %i m' %min_length)

# only keep cells with segments greater than 100m long to reduce computation time
# and as short segments should contribute less seepage and it creates a more continuous line
grid_sfr_final = grid_sfr[grid_sfr.length_m>min_length]

# %%
fig,ax=plt.subplots(figsize=(6,4))
grid_sfr.plot(ax=ax, color='red')
grid_sfr_final.plot(ax=ax)
# grid_sfr_final.plot(x='dist_m',y='z_min', kind='line')

# %%
grid_sfr_final.to_file(sfr_dir+'/final_grid_sfr/'+label.replace(' ','_')+'_sfr.shp')
# grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')

