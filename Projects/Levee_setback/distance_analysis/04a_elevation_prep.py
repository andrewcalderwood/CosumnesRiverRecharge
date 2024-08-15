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

# %%
# standard python utilities
import os
from os.path import dirname, basename, exists, join
import sys
import glob
import pandas as pd
import numpy as np
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
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
    
git_dir = join(doc_dir,'GitHub','CosumnesRiverRecharge')
# dir of all gwfm data
gwfm_dir = join(dirname(doc_dir),'Box/research_cosumnes/GWFlowModel')




# %%
def add_path(fxn_dir):
    if fxn_dir not in sys.path:
        sys.path.append(fxn_dir)
        
add_path(doc_dir+'/GitHub/CosumnesRiverRecharge/python_utilities')

# %%
from map_cln import gdf_bnds, pnt_2_tup, lab_pnt, plt_cln


# %%
# set box directory for output figures and data
box_dir = gwfm_dir+'/Levee_setback/levee_setback_distance_analysis/'

fig_dir = box_dir+'figures/'
data_dir = box_dir+'data_output/'

chan_dir = box_dir+'channel_data/'
gis_dir = chan_dir+'GIS/'

# %%
ext_dir = 'F:/WRDAPP'
c_dir = 'C:/WRDAPP'

if os.path.exists(ext_dir):
    loadpth = ext_dir 
elif os.path.exists(c_dir):
    loadpth = c_dir 

loadpth = loadpth +'/GWFlowModel/Cosumnes/levee_setback/'
model_ws = loadpth+'flood_depth_analysis'

nrow = 100
ncol = 230
nlay = 1
delr = 200
delc = 200

# %%
# burn river shapefile into the 10 meter dem and then read it out to find the cells where it is
# # Full size dem of northern sac valley
raster_name = gwfm_dir+"/DEM_data/USGS_ten_meter_dem/modeldomain_10m_transformed.tif"



# %%
# grid_sfr = gpd.read_file(gwfm_dir+'/SFR_data/final_grid_sfr/grid_sfr.shp')
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')

m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')
sfr_dir = join(gwfm_dir, 'SFR_data')
grid_sfr = gpd.read_file(join(sfr_dir,'final_grid_sfr/grid_sfr.shp'))
# load sacramento river, creeks
rivers = gpd.read_file(join(sfr_dir,'Sac_valley_rivers/Sac_valley_rivers.shp'))
cr = gpd.overlay(rivers.loc[rivers.GNIS_Name=='Cosumnes River'].to_crs('epsg:32610'), m_domain)

# %% [markdown]
# There is no good way to adjust the stream line because in the end the channel will be off from the cross-sections. This means it is best to choose at a coarse interval (e.g., 2000 m) and then correct the center alignment.

# %%
from shapely.geometry import MultiLineString, LineString, Point, shape, mapping
from shapely.ops import linemerge
import fiona

cr_line = MultiLineString(cr.geometry.values)

# %%
setbacks = np.arange(0,3400,200)
import h5py
f = h5py.File(join(chan_dir, 'setback_locs.hdf5'), "r")
local_str_setbacks = f['setbacks']['local'][:]
str_setbacks = f['setbacks']['regional'][:]

f.close()


# %% [markdown]
# # Prepare Elevation data inputs

# %% [markdown]
# ## Parallel XS

# %%


geom = shape(cr_line)

# # length of the LineString
length = int(geom.length)
sfr_sp = gpd.GeoDataFrame(pd.DataFrame(np.zeros((length,1)),columns=['id']))
sfr_sp['point'] = shapely.geometry.Point(0,0)
sfr_sp = sfr_sp.set_geometry('point')
# enumerate keeps track of the count in addition to the object being iterated on
for i, distance in enumerate(range(0, length, 1)):
    point = geom.interpolate(distance)
    sfr_sp.loc[i,'geometry'] = point
    sfr_sp.loc[i,'id'] = i

# %%
sfr_sp['Easting'] = sfr_sp.geometry.x.values
sfr_sp['Northing'] = sfr_sp.geometry.y.values

point = sfr_sp.loc[:,['Easting','Northing']].values

with rasterio.open(raster_name) as src:
    sfr_sp['z_ft'] = [sample[0] for sample in src.sample(point)]
    
sfr_sp['z_m'] = sfr_sp.loc[:,'z_ft']*0.3048


# %% [markdown]
# Find XS every 1000 meters avoids too much overlap in XS,  after calculating depth there is some discontinuity but WSE is uniform with slope, so the numebr could be reduced to 2000m which would also aid muskingum-cunge routing requirements (coarser time step allowed).   
# I also realized that I should simplify the river feature to linearize it and avoid XS occuring on oddities where there is a large turn in the channel for a short distance. Adding a 500m tolerance simplify additionally reduces overlap and produces good XS locations/angles.
# - Although the XS are wide enough, there are some that are offset from the channel center and if we include an extra 200 m then we can adjust for the true channel center

# %%

def make_transects(geom, dline = 2000, xs_width = 3305*2):
    """ 
    Function to make regular transects for a line object (geom : shapely polyline) 
    at some longitudinal distance (dline : distance in meters)
    max width needed is 3200 m +100m from original channel + 5 to allow full 3200 both direction
    """
    # tolerance is distance which simplified points can be from the original
    geom = geom.simplify(500)

    # how often to interpolate a point, 200 m matches model grid
    dline_long = np.copy(dline)
    # # length of the LineString
    length = int(geom.length)

    num_xs = np.floor(length/dline).astype(int)
    transects = pd.DataFrame(np.zeros((num_xs,1)), columns = ['line'])
    transects['geometry'] = LineString([(0,0),(0,1)]) #initiate LineString geometry column


    for i, distance in enumerate(range(0, int(length), dline)):
        short_line = LineString([geom.interpolate(distance),geom.interpolate(distance+dline)])
        geom_left = short_line.parallel_offset(xs_width/2,'left', resolution = 32, join_style = 2)
        geom_right = short_line.parallel_offset(xs_width/2,'right', resolution = 32, join_style = 2)
        perp_line = LineString([list(geom_left.coords)[0], list(geom_right.coords)[0]])
        transects.loc[i,'geometry'] = perp_line

    # save output
    transg = gpd.GeoDataFrame(transects)
    transg = transg.drop_duplicates('geometry')
    transg['line'] = np.arange(0,len(transg))
    transg.crs = 'epsg:32610'
    return(transg)



# linemerge is essential because it joins linesegments by matching (x,y) so the line is continuous
geom = linemerge(cr.geometry.unary_union)

# %%
xs_width = 3305*2 # original w/o adjustment for channel center
xs_width = 3705*2 # add 400 m to both sides to then cut off later for finding channel center

transg200 = make_transects(geom = linemerge(cr.geometry.unary_union), dline = 200, xs_width = xs_width)
# drop transects that are less than 25% in the domain
transg200 = transg200[gpd.overlay(transg200, m_domain).geometry.length > xs_width*0.25]
transg200['line'] = np.arange(0,len(transg200))
transg200.to_file(gis_dir+'/transect_lines_subsegments.shp')

transg = make_transects(geom = linemerge(cr.geometry.unary_union), dline = 2000, xs_width = xs_width)
# drop transects that are less than 25% in the domain
transg = transg[gpd.overlay(transg, m_domain).geometry.length > xs_width*0.25]
transg['line'] = np.arange(0,len(transg))
transg.to_file(gis_dir+'/transect_lines.shp')

# check cross section lines are parallel
fig,ax=plt.subplots(figsize=(6,6))
transg200.plot(color='red', ax=ax)
transg.plot(ax=ax)
cr.plot(ax=ax)
m_domain.plot(ax=ax, color='None')

# %%
# lines for plotting figures
transg = make_transects(geom = linemerge(cr.geometry.unary_union), dline = 2000, xs_width = 3300*2)
# drop transects that are less than 25% in the domain
transg = transg[gpd.overlay(transg, m_domain).geometry.length > 3300*2*0.25]
transg['line'] = np.arange(0,len(transg))
# transg.to_file(gis_dir+'/transect_lines_3300.shp')

# %%
# Combine all lines in gdf_lines_2 into a single geometry (unary_union)
splitter = transg.unary_union

# Function to split each line in gdf_lines_1 by the splitter
def split_line_by_splitter(line, splitter):
    if line.intersects(splitter):
        # Split the line by the splitter
        split_result = split(line, splitter)
        # If the result is a MultiLineString, split into separate LineStrings
        if isinstance(split_result, MultiLineString):
            return list(split_result.geoms)
        else:
            return [split_result]
    else:
        return [line]

# Apply the splitting function to all lines in gdf_lines_1



# %%
from shapely.ops import split

# Apply the splitting function to all lines in gdf_lines_1
split_geometries = transg['geometry'].apply(lambda x: split_line_by_splitter(x, splitter))

# Flatten the list of lists into a single list of geometries
split_lines = [geom for sublist in split_geometries for geom in sublist]

# Create a new GeoDataFrame with the split lines
gdf_split_result = gpd.GeoDataFrame(geometry=split_lines, crs=gdf_lines_1.crs)


# %%
# need to crop transects where they overlap

# using buffer and overlay doesn't work because the buffer doesn't sufficiently overlap the cross-sections
# need to use a line split and drop the that doesn't touch the river line
# create uni-direction buffer
transg_buf_1way = transg.copy()

# transg_buf_1way.geometry = transg_buf_1way.buffer(-400, single_sided=True)
# # first remove the original transect from the buffer
# transg_buf_1way = transg_buf_1way.overlay(transg, how='difference')
# # now create an udpated cross-section with the overlapping part removed
# transg_new = transg.copy()
# transg_new = transg_new.overlay(transg_buf_1way, how='difference')

fig,ax = plt.subplots()
transg_buf_1way.plot(ax=ax,color='red')
transg.plot(ax=ax)
transg_new.plot(color='limegreen', ax=ax)


# %%
def transect2pts(transg, dline = 10, xs_width = 3305*2):
    """ 
    Convert line transects (transg : geodataframe) to points 
    at a regular interval along the transect (dline : distance in meters)
    """

    xs_all = gpd.GeoDataFrame(pd.DataFrame(columns=['xs_num','dist_from_right_m','geometry']))

    for j in np.arange(0,len(transg)):
        xs = gpd.GeoDataFrame(pd.DataFrame(np.zeros((int(xs_width/dline),2)), columns=['xs_num','dist_from_right_m']))
        xs['point'] = Point([(0,0)])
        xs = xs.set_geometry('point')
        xs['xs_num'] = j

        # pick one geometry at a time
        geom = transg.iloc[j].geometry

        # # # length of the LineString
        length = int(geom.length)
        # create discrete points for each lien
        for i, distance in enumerate(range(0, int(length), dline)):
            point = geom.interpolate(distance)
            xs.loc[i,'geometry'] = point
            xs.loc[i,'dist_from_right_m'] = i
        # append individual cross section to all dataframe
    #     xs_all = xs_all.append(xs)
        xs_all = pd.concat([xs_all, xs], axis=0, join='outer', ignore_index=True)
    # remove na values
    xs_all = xs_all.dropna(subset=['xs_num'])
    return(xs_all)
xs_all = transect2pts(transg, dline = 10, xs_width = xs_width)


# %%
# Sample DEM for every XS point for an elevation

def sample_pts(xs_all, raster_name):
    """
    Given a geodataframe of points sample elevations from a raster
    """
    src = rasterio.open(raster_name)
    xmin, xmax = src.bounds[0], src.bounds[2]
    ymin, ymax = src.bounds[1], src.bounds[3]

    xs_all.crs='epsg:32610'
    xs_all['Easting'] = xs_all.geometry.x
    xs_all['Northing'] = xs_all.geometry.y
    # filter to points that truly overlap the DEM
    xs_all = xs_all.copy().cx[xmin:xmax, ymin:ymax]

    # point = xs_all.loc[:,['Easting','Northing']].values # old
    point = list(zip(xs_all.geometry.x, xs_all.geometry.y)) # same error

    print(os.path.basename(raster_name))
    with rasterio.open(raster_name) as src:
        xs_all['z_m'] = [sample[0] for sample in list(src.sample(point))]

    # convert distance from right from 1/10 of meters because points were every 10 meters
    # the distance was set with the index but really it should be 3200 meters not just 320 meters
    xs_all.dist_from_right_m *= 10
    # remove any NA values picked up from DEM raster
    xs_all.loc[xs_all['z_m'] == src.meta['nodata'], ['z_m']] = np.nan
    xs_all.index = np.arange(0,len(xs_all))
    return(xs_all)
xs_all = sample_pts(xs_all, raster_name)

# %% [markdown]
# ## XS Cleaning
# There isn't going to be a clear way to distinguish the levees and the XS shows bumpiness. The best solution will be to assume that levee setback will involve a cut and fill approach such that the ground surface elevation after setback is the mean of the elevation from before, but then this won't allow variable flooding based on elevation...
#
# It will take looking at the XS in different regions.
# It might be doable to fill in the channel as Sierra suggested which would raise the river above the levee and then install new "levee walls" at the desired distance. Filling in just the channel requires: 1. Going to the center line 2. go out some distance to account for channel width 3. set those values as a fraction of the levee height to insure more overbank flooding

# %%
# calculate elevations for subsegments (200 m)
xs_all200 = transect2pts(transg200, dline = 10, xs_width = 3305*2)
xs_all200 = sample_pts(xs_all200, raster_name)

xs_all200.drop(['geometry'],axis=1).to_csv(chan_dir + 'XS_point_subsegments_elevations.csv')

# convert to dataframe for easier plotting
xs_all_df200 = pd.DataFrame(xs_all200)
# pivot based on XS number and save only elevation in z_m
xs_all_df200 = xs_all_df200.pivot_table(index='dist_from_right_m',columns='xs_num',values='z_m')
xs_all_df200.to_csv(chan_dir + 'Elevation_by_XS_number_subsegments_meters.csv')

# %%
# elevations and XS for every 2 km
xs_all.drop(['geometry'],axis=1).to_csv(chan_dir+'XS_point_elevations.csv', index=False)

# convert to dataframe for easier plotting
xs_all_df = pd.DataFrame(xs_all)
# pivot based on XS number and save only elevation in z_m
xs_all_df = xs_all_df.pivot_table(index='dist_from_right_m',columns='xs_num',values='z_m')
# xs_all_df.to_csv(chan_dir+'Elevation_by_XS_number_meters.csv')

# %%
dline=10
roll_window=600
channel_middle = int(len(xs_all_df)/2)*dline
channel_bool = (xs_all_df.index >= channel_middle - (roll_window/2))& (xs_all_df.index <= channel_middle + (roll_window/2))    
min_idx =  xs_all_df[channel_bool].idxmin()
# min_idx

# %%

fig,ax = plt.subplots(5,1, figsize=(5,12),sharex=True)
for ng in np.arange(0,5):
    nseg = np.arange(ng*5,5*(ng+1))
    xs_chk = xs_all_df[nseg]
    xs_chk[channel_bool].plot(ax=ax[ng])
    for n, ns in enumerate(nseg):
        ax[ng].plot(min_idx.iloc[n], xs_chk.loc[min_idx.iloc[n], ns],marker='x', color='black')
    
fig.tight_layout(h_pad=0.1)


# %%
f_width = 3305 # final width to extract on each side of the new min
xs_all_fix = pd.DataFrame()
for nseg in np.arange(0,xs_all_df.shape[1]):
# for nseg in [2]:
    # index the widths
    xs_fix = xs_all_df.loc[min_idx.loc[nseg]-f_width:min_idx.loc[nseg]+f_width, nseg].copy()
    # assign the new distances
    xs_fix.index -= min_idx.loc[nseg]
    xs_all_fix = pd.concat((xs_all_fix, xs_fix), axis=1)
xs_all_fix.index.name='dist_from_center_m'

# %%
xs_all_fix.to_csv(chan_dir+'Elevation_by_XS_number_meters.csv')

# %%

# xs_all_fix.plot(legend=False)
# xs_all_df.plot()
# plt.legend(ncol=4, loc=(1,0.01), title='XS Number')


# %% [markdown]
# # XS adjustments

# %% [markdown]
# The process will be too complicated too directly edit the raster, so individual cross-sections will be edited by applying a rolling mean in the floodplain in the 200-400 m around the main channel to "flatten" any levee like structure. Then new 'levee walls' will be placed at the desired setback distance.
#
#  Helen suggested looking at the historical 100 year flood to estimate the levee height they would have designed for, or look at the physical characteristics to design it. Levee construction began in 1930s, biggest flood was in 1907 at 71,000 cfs (estiamted) with a gage height of 16 ft at Michigan Bar so levee's would most likely be designed to at least 16 ft.
#
#  For each levee setback we should imagine that the channel is partially filled in, but rather than uniformly removing all levee walls within the floodplain, we breach the levee wall in 200 m sections for every 2km to allow flood water to inundate the floodplain and then new levees are constructed at the specified distance to bring ground level to 16 ft above channel bottom. In real life this height could probably be steadily lowered as there is so much room in the floodplain to accomdate the volume that a lower levee height suffices.
#  1. Rolling mean is used to fill in channel with earth from levee or nearby.
#  2. New levees are continuously installed at specified distance
#
# **The issue with the rolling mean is that this doesn't translate properly to the DEM elevation comparison and creates a weird artifact where the depth increases with larger cross-sections which doesn't make sense conceptually.**

# %%
xs_all = pd.read_csv(chan_dir+'XS_point_elevations.csv')

xs_all_df = pd.read_csv(chan_dir+'Elevation_by_XS_number_meters.csv',index_col='dist_from_center_m')


def xs_channel_smooth(xs_all_df, roll_window = 400, dline = 10):
    """ 
    Apply a rolling mean to center window (roll_window) of channel (xs_all_df) 
    to represent incision infill and reduction of bank/levee height,
    40 observations is equal to 400 meters
    """
    xs_levee_smooth = xs_all_df.copy()
    xs_roll_mean = xs_levee_smooth.rolling(int(roll_window/dline), center=True).mean()
    channel_middle = int(len(xs_all_df)/2)*dline
    channel_bool = (xs_all_df.index >= channel_middle - (roll_window/2))& (xs_all_df.index <= channel_middle + (roll_window/2))
    xs_levee_smooth.loc[channel_bool,:] = xs_roll_mean.loc[channel_bool,:]
    return(xs_levee_smooth)

# %%
# not sure if needed anymore
# xs_levee_smooth200 = xs_channel_smooth(xs_all_df200, roll_window = 400, dline = 10)
# # save for flow-recharge notebook
# xs_levee_smooth200.to_csv(chan_dir+'xs_levee_smooth_subsegments.csv')


# %%

xs_levee_smooth = xs_channel_smooth(xs_all_df, roll_window = 400, dline = 10)
# # save for flow-recharge notebook
xs_levee_smooth.to_csv(chan_dir+'xs_levee_smooth.csv')

# %%
fig,ax = plt.subplots(figsize=(12,6))

# xs_levee_smooth.plot(ax=ax)
# ax.legend(ncol=4, loc=(1,0.01), title='XS Number')
n=26

xs_all_df.iloc[:,n].plot(ax=ax, label='XS')
xs_levee_smooth.iloc[:,n].plot(ax=ax, label='Smooth')
plt.legend()


# %% [markdown]
# ## Rating curves for each segment and setback (50-points)
# 50 point rating curves to match modflow format.

# %%
xs_all_df = pd.read_csv(chan_dir+'Elevation_by_XS_number_meters.csv', index_col='dist_from_center_m')

xs_levee_smooth = pd.read_csv(chan_dir+'xs_levee_smooth.csv', index_col='dist_from_center_m')

# %%
# load flood typology characteristics (based on daily data 1908 - 2014) - median values 
#"cms_pk" for peak discharge, "pk_loc" for time to peak, and "log_no_d" for duration
flood_type = pd.read_csv(join(box_dir, 'whipple_grp6_w97ftmedians.csv'),index_col='Group.1')


# %%
# from importlib import reload
# import muskingum_recharge
# reload(muskingum_recharge)
from muskingum_recharge import min_Q, mannings, calc_depth_arr, gridded_interpolation, xs_setback, mannings_v


# %%
def get_mins(xs_all_cln):
    # find minimum from channel center
    xs_mins = xs_all_cln.loc[-100:100].min(axis=0)
    xs_mins.index = xs_mins.index.astype(int)
    # xs_mins.interpolate(method='linear').plot()
    # 4 is required to avoid negatives
    slope = xs_mins.diff().rolling(4, center=True, closed='right').mean().bfill().ffill()/2000*-1
    adj_xs_mins = np.append(xs_mins[0], (xs_mins[0]-slope.cumsum()*2000))
    return(adj_xs_mins, slope)



# %%
## having issues with NAs in all cels of a XS seemed to be driven by the levee smooth mixed in

# %%
def make_rating_curve(xs_all_cln):
    adj_xs_mins, slope = get_mins(xs_all_cln)
    
    segs = np.arange(0, xs_all_cln.shape[1])
    xs_flow_all = pd.DataFrame(columns=['nseg','setback','depth_m','flow_cms']).set_index(['nseg','setback','depth_m'])
    # original code takes 5-10 seconds, the slowness is from the loops and mannings equation, not the dataframe setup
    n= 0.048
    # iterate over cross-section segments
    for nseg in segs:
        # nseg = 15
        # df = xs_all_cln[str(nseg)]
        # iterate over the stream segments
        for setback in setbacks:
        # setback = 3200
            # maximum depth is tallest height of cross-section minus lowest point, flow above will run out of the channel
            xs_elevs = xs_setback(xs_all_cln.iloc[:,nseg].copy(), setback, 30)
            dmax = np.min((xs_elevs.max()-xs_elevs.min(), 20)) # should not simulate above 20 m (30 ft depth)
            dmax = np.max((dmax, 10)) # require at least 10 m to ensure sufficient flow calculation
            # calculate more flow points than needed then chop off above extreme flows
            for d in np.linspace(0.01, dmax, 50):
                flow = mannings(d, xs_elevs, n, slope.iloc[nseg])
                xs_flow_all.loc[(nseg, setback, d),'flow_cms'] = flow
    # return calculated rating curves
    return xs_flow_all.reset_index()





# %%
xs_flow_all = make_rating_curve(xs_all_df)
xs_flow_all = xs_flow_all.dropna(subset='flow_cms')
# xs_flow_all.to_csv(join(chan_dir,'all_xs_50pt_rating_curves.csv'), index=False)

# xs_flow_all = make_rating_curve(xs_levee_smooth)
# xs_flow_all.to_csv(join(chan_dir,'all_xs_smooth_50pt_rating_curves.csv'), index=False)

# %%
# rolling mean to smooth out rating curves and create all unique flow-depths
xs_flow_all_mean = xs_flow_all.groupby(['nseg','setback']).rolling(window=6, min_periods=0, center=False).mean().reset_index(['nseg','setback'])
xs_flow_all_mean.to_csv(join(chan_dir,'all_xs_50pt_rating_curves.csv'), index=False)


# %%
# xs_elevs = xs_setback(xs_all_df.iloc[:,nseg].copy(), setback, 30)
# xs_elevs.plot()

# %%
nseg=9
xs_flow_all_mean[(xs_flow_all.nseg==nseg)&(xs_flow_all.setback==0)].max()

# %%
# abnormalities in the rating curve mean that flow can decrease with depth and create non-unique flow-depths
# apply smoothing to alleviate this
fig,ax=plt.subplots()
nseg = 9
setback=0

xs_flow_all[(xs_flow_all.nseg==nseg)&(xs_flow_all.setback==0)].plot('depth_m','flow_cms',ax=ax)
xs_flow_all_mean[(xs_flow_all.nseg==nseg)&(xs_flow_all.setback==0)].plot('depth_m','flow_cms',ax=ax)


# %% [markdown]
# ## Given flow return depth
#

# %%

def depth_match(seg_flow, flow):
    """ Given a XS (nseg, setback) return the expected depth (m) given a flow (cms)"""
    # find flows above and below the input flow
    flow_diff = (seg_flow.flow_cms-flow)
    f_high = flow_diff[flow_diff>0].argsort().index[0]
    f_low = flow_diff[flow_diff<0].argsort().index[-1]
    match_d = seg_flow.loc[[f_low, f_high]].sort_values('flow_cms')
    # linearly interpolate to calculate exact depth
    flow_slope = (match_d.iloc[1].flow_cms-match_d.iloc[0].flow_cms)/(match_d.iloc[1].depth_m-match_d.iloc[0].depth_m)
    out_depth = match_d.iloc[0].depth_m + (flow-match_d.iloc[0].flow_cms)/flow_slope
    return(out_depth)


nseg = 10
setback=3200
seg_flow = xs_flow_all[(xs_flow_all.nseg==nseg)&(xs_flow_all.setback==setback)]
# seg_flow
depth_match(seg_flow, flow=500)

# %% [markdown]
# # Identify XS to grid cells

# %%
xs_gdf = gpd.GeoDataFrame(xs_all, geometry = gpd.points_from_xy(xs_all.Easting, xs_all.Northing), crs='epsg:32610')

# %%
transg = gpd.read_file(gis_dir+'/transect_lines.shp')
def transect2arr(transg, grid_p):
    """
    Given transects (transg : geodataframe of polylines) identify the corresponding
    model grid cells (grid_p:geodataframe of polygons) to create an array mapping XS number
    """
    transg_box = transg.copy()
    # upstream buffer only, seepage connects upstream to downstream cross-section
    transg_box['geometry'] = transg_box.buffer(-3000, single_sided=True)

    # overlay transect with grid to identify XS number for each row, col
    grid_xs = gpd.sjoin(grid_p, transg_box)
    # dissolve and find minimum assuming upstream reach is already covered by previous XS
    grid_xs = grid_xs.dissolve('node','min')
    # map XS number to an array for identifying which recharge goes to which XS number
    xs_arr = np.full((nrow, ncol), np.NaN)
    xs_arr[grid_xs.row-1, grid_xs.column-1] = grid_xs.line
    return(xs_arr)

xs_arr = transect2arr(transg, grid_p)



# %%
plt.imshow(xs_arr)
print(np.unique(xs_arr))
np.savetxt(chan_dir+'XS_num_grid_reference.tsv', xs_arr, delimiter='\t')

# %% [markdown]
# # Calculate the XS minimum elevations

# %%
transg200 = gpd.read_file(gis_dir+'/transect_lines_subsegments.shp')

xs_arr200 = transect2arr(transg200, grid_p)

# %%
zs = gpd.read_file(gwfm_dir+'/DIS_data/grid_elevation_m_statistics.shp')
# columns with different quantiles 0 to 100% of elevation
q_cols = zs.columns[zs.columns.str.contains('perc')]
df_elevs = zs[q_cols]

# convert quantile dataframe to a 3D array
arr_elev = np.zeros((df_elevs.shape[1], zs.row.max(),zs.column.max()))
for n in np.arange(0,df_elevs.shape[1]):
    arr_elev[n, zs.row-1, zs.column-1] = df_elevs.iloc[:,n]

# %% [markdown]
# There needs to be a choice between having a linear, continuous water surface or having true ground elevations. Having a linear continuous elevation is a requirement with the kinematic wave when assuming continually downsloping channel and thus the slope being continually sloping. However in reality the cosumnes has an irregular flood elevation due to choke points, velocity changes, etc so it might be acceptable to assume a singular depth applies over irregular elevations with some discontinuities. 
# An intermediate option would be a mostly linear downsloping channel but the bottom elevations after cleaning align closely.
#
# An issue with creating a linearly sloping water surface is that when mapping out inundation there are cells that will be consistently inundated if any interpolated stream bed (xs_mins_arr) is above the original minimum (arr_elev 0th percentile).
# I updated the subsegments minimum to follow this rule.

# %%
# find the corresponding minimum segment elevations from the raster sampled minimums
arr_elev_mins = np.zeros(len(np.unique(xs_arr200)[:-1]))
for n in np.unique(xs_arr200)[:-1].astype(int):
    arr_elev_mins[n] = np.min(arr_elev[0,xs_arr200==n])


# %%
def linearize_xs_min(xs_mins, dline, window=2):
    """ 
    Clean up the XS minimums (xs_mins:dataframe) so they are always decreasing linearly,
    this ensures that the constant depth applied will lead to an approximately linear water surface elevation
    dline : distance between XS, used for slope 
    window : number of upstream or downstream points to use in rolling mean
    """
    # find minimum value in XS related to thalweg
    xs = pd.DataFrame(xs_mins, columns=['z_min'])
    #roling mean of 2 window centered removes any negative slope (manually tested)
    xs['z_min_cln'] = xs.rolling(window, center=False).mean()
    # if rolling mean set values greater than 5% above previous value, reset
    xs.loc[xs.z_min_cln > xs.z_min*1.01, 'z_min_cln'] =  xs.loc[xs.z_min_cln > xs.z_min*1.01, 'z_min']
    xs['z_min_cln'] = xs['z_min_cln'].rolling(window, center=False).mean()

    # calculate slope and fill NAs, fill slope with nearby
    z_cln_diff = xs.z_min_cln.diff().bfill()
    xs['slope'] = z_cln_diff.abs()/dline
    # # correct slope less than 1E-4
    xs.loc[xs.slope<1E-4,'slope'] = 1E-4
    # average out slope
    xs.loc[:,'slope'] = xs.loc[:,'slope'].rolling(window*4, center=False).min()
    xs.loc[:,'slope'] = xs.loc[:,'slope'].bfill().ffill()

    for i in xs.index[-2::-1]:
    # fill NAs due to rolling mean, with backward filling
        if np.isnan(xs.loc[i,'z_min_cln']):
            xs.loc[i,'z_min_cln'] = xs.loc[i+1,'z_min_cln'] + xs.loc[i+1,'slope']*dline
        
#     fix str bot so all is downward sloping
    for i in xs.index[:-1]:
        if xs.loc[i+1,'z_min_cln'] >= xs.loc[i,'z_min_cln']:
            xs.loc[i+1,'z_min_cln'] = xs.loc[i,'z_min_cln'] - xs.loc[i,'slope']*dline

    return(xs)
xs_mins = np.nanmin(xs_levee_smooth.loc[3100:3300], axis=0)
xs = linearize_xs_min(xs_mins, dline=2000, window=2)

fig,ax=plt.subplots()

xs_plt = xs.set_index(np.arange(0,2000*xs.shape[0],2000))

xs_mins = np.nanmin(xs_levee_smooth200.loc[3100:3300], axis=0)
xs = linearize_xs_min(xs_mins, dline=200, window=5)

xs_plt = xs.set_index(np.arange(0,200*xs.shape[0],200))
xs_plt.plot(y='z_min',label='Subseg Min',ax=ax)
xs_plt.plot(y='z_min_cln',label='Subseg Clean',ax=ax)

# version based on raster minimum samples by grid cell
r = linearize_xs_min(arr_elev_mins, dline=200, window=5)
r_plt = r.set_index(np.arange(0,200*r.shape[0],200))
r_plt.plot(y='z_min',label='Raster Min',ax=ax)
r_plt['z_min_adj'] = np.where(r_plt.z_min_cln>r_plt.z_min, r_plt.z_min, r_plt.z_min_cln)
r_plt.plot(y='z_min_adj',label='Raster Clean',ax=ax)

plt.legend()

# %%
# create an array of minimum eleations tied to the cleaning mins above spread out to each subsegment
xs_mins_arr200 = np.full((nrow,ncol),np.nan)
for n in np.arange(0, len(xs)):
#     xs_mins_arr200[xs_arr200==n] = xs.z_min_cln[n] # sampled XS minimum
    xs_mins_arr200[xs_arr200==n] = r_plt.z_min_adj.values[n] # sampled grid cell minimum
    
np.savetxt(chan_dir+'subsegments_xs_mins.tsv', xs_mins_arr200, delimiter='\t')

# %%
xs_test = xs_levee_smooth.copy()
# rather than saving a version of the xs for each setback just add this adjustment
n=5
# the channel should fall within the center 400m so that minimum can be used to set new levee height
thalweg = xs_test.loc[channel_middle-roll_window/2:channel_middle+roll_window/2].min()
# check to see XS elevation at setback distance to determine if it should be raised to needed levee height
levee_diff = xs_test.loc[3100-setbacks[n],:] - thalweg
# where XS height is less than 20 ft above channel bottom then raise to 20 ft above
xs_test.loc[3100-setbacks[n],:][levee_diff< 20*0.3048] = thalweg + 16*0.3048
levee_diff = xs_test.loc[3300+setbacks[n],:] - thalweg
xs_test.loc[3300+setbacks[n],:][levee_diff< 20*0.3048] = thalweg + 16*0.3048

xs_test.iloc[:,24].plot()


fig,ax=plt.subplots(figsize=(6,3))
xs_levee_smooth.min().plot(ax=ax, label='Adj XS')
xs_all_df.min().plot(ax=ax, color='red',label='XS')
grid_sfr.assign(river_km = np.cumsum(grid_sfr.length_m/dline_long)).plot(x='river_km',y='z',kind='line',ax=ax, label='Model')
plt.legend()

# dem data for cropping above land surface
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_linear.tsv')


