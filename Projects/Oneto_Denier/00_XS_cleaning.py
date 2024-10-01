# ---
# jupyter:
#   jupytext:
#     formats: py:percent
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

# %% editable=true slideshow={"slide_type": ""}
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

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %% editable=true slideshow={"slide_type": ""}
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study
proj_dir = gwfm_dir + '/Oneto_Denier/'
dat_dir = proj_dir+'Stream_level_data/'

sfr_dir = gwfm_dir+'/SFR_data/'
uzf_dir = gwfm_dir+'/UZF_data/'


# %% editable=true slideshow={"slide_type": ""}
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

# from mf_utility import get_layer_from_elev
from mf_utility import get_layer_from_elev, param_load

from map_cln import gdf_bnds, plt_cln


# %% [markdown]
# # SFR

# %%
model_nam = 'oneto_denier'

# %%
# write modelgrid to get updated row and col numbers specific to the child grid
grid_dir = join(gwfm_dir, 'DIS_data/streambed_seepage/grid')
grid_fn = join(grid_dir, 'inset_oneto_denier','rm_only_grid.shp')

# m.modelgrid.write_shapefile(grid_fn)
grid_p = gpd.read_file(grid_fn)
grid_p.crs = 'epsg:32610'

# %%
delr, delc = (100, 100)

# %% [markdown]
# ## XS pre-processing

# %%
# cross sections sampled using NHD lines at regular 100 m intervals (not aligned with any grid)
xs_all = pd.read_csv(dat_dir+'XS_point_elevations.csv',index_col=0)
xs_all = gpd.GeoDataFrame(xs_all,geometry = gpd.points_from_xy(xs_all.Easting,xs_all.Northing), crs='epsg:32610')

# find XS that are in the modeled domain by thalweg point
thalweg = xs_all[xs_all.dist_from_right_m==100]
thalweg = gpd.overlay(thalweg, grid_p)
# thalweg = gpd.sjoin_nearest(grid_p, thalweg, how='inner')
# thalweg = thalweg.cx[xmin:xmax, ymin:ymax]
# with XS every 100m I need to choose whether the first or second is used in a cell
# thalweg = thalweg.dissolve(by='node', aggfunc='first')

# pivot based on XS number and save only elevation in z_m
xs_all_df = pd.read_csv(dat_dir+'Elevation_by_XS_number_meters.csv',index_col=0)
xs_all_df = xs_all_df.dropna(axis=0,how='any')

# filter XS by those that are within the domain bounds
xs_all = xs_all[xs_all.xs_num.isin(thalweg.xs_num.values)]
xs_all_df = xs_all_df.loc[:, thalweg.xs_num.astype(str)]

# renumber XS
thalweg.xs_num = np.arange(0,thalweg.shape[0])
xs_all.xs_num = np.repeat(thalweg.xs_num.values,xs_all.dist_from_right_m.max()+1)
xs_all_df.columns = thalweg.xs_num

# with 100 m cells we should crop to a 100m XS instead of 200 m
center = int(np.mean(xs_all_df.index.values))
offset = int(xs_all_df.shape[0]/4)
xs_all_df = xs_all_df.loc[center-offset:center+offset]

# %%
# load channel points to identify sfr reaches
sfr_sp = pd.read_csv(join(proj_dir, 'Stream_level_data', 'stream_line_elevation.csv'))
sfr_sp = gpd.GeoDataFrame(sfr_sp, geometry = gpd.points_from_xy(sfr_sp.Easting, sfr_sp.Northing), crs = xs_all.crs)
sfr_grid_all = gpd.overlay(sfr_sp, grid_p)
# sfr_grid = gpd.sjoin_nearest(sfr_sp, grid_p, how='inner')
sfr_grid_all['rch_len'] = 1 # spacing of points was 1 meter

# %%
sfr_grid_all.loc[sfr_grid_all.z_m<-1E3, 'z_m'] = np.nan


# %%
grp_cols = ['node','row','column']
sfr_grid = sfr_grid_all.groupby(grp_cols).mean(numeric_only=True).reset_index()
sfr_grid['rch_len'] = sfr_grid_all.groupby(grp_cols).sum(numeric_only=True).reset_index()['rch_len']
sfr_grid = grid_p.merge(sfr_grid)
sfr_grid.z_m = sfr_grid.z_m.interpolate('linear').copy()

# %%
# switching to only segments longer than 50 m 
min_len = 50
sfr_reach = sfr_grid[sfr_grid.rch_len>min_len]
print('Final reaches', sfr_reach.shape[0], 'Reaches dropped',sfr_grid[sfr_grid.rch_len<=min_len].shape[0])
# sfr_grid[sfr_grid.rch_len>min_len].plot()

# %%
# could use the error between cross-sections as a way to approximate for how long a cross-section is valid
# it's best to use each cross-section statically against all others
from sklearn.metrics import mean_squared_error

xs_mse = np.zeros((xs_all_df.shape[1], xs_all_df.shape[1]))
for k in np.arange(0, xs_all_df.shape[1]):
    for n in np.arange(0, xs_all_df.shape[1]):
        # xs_diff = xs_all_df.iloc[:,n] - xs_all_df.iloc[:,n+1]
        xs_mse[k, n] = np.sqrt(mean_squared_error(xs_all_df.iloc[:,k], xs_all_df.iloc[:,n]))

# %%
# plotting an individual cross-section minus all cross-sections shows for how long it is valid (RMSE >2 m)
for n in np.arange(0,xs_all_df.shape[1])[:20]:
# for n in np.arange(0,xs_all_df.shape[1], 20):
    plt.plot(xs_mse[n])

# it shows on avg for how far a transect is valid, the middle transects are valid for much of the domain
print('Avg RMSE <1 is %.i segments' %np.median(np.sum(xs_mse<2, axis=1)))
# plt.plot(np.sum(xs_mse<1, axis=1))

# identify cross-sections that have less than 1.5 m error as a group
# then filter for continuous segments with an allowance for 1 cross-section segment that didn't fit

# %%
# group cross-sections into 1 km reaches to uspcale
xs_num = xs_all_df.columns.values
seg_len = 4
nseg_grp = int(len(xs_num)/seg_len)


# %% [markdown]
# The cross-section mean shows a slight incision toward the downstream segments. Averaging the cross-sections every 10 segments may not be appropriate because it seems to lead to less incision (less rectangular) and a more sloping cross-section. Every 5 segments kept more uniform cross-section but still there was a loss of channel area, every 4 does slightly better as well. Stick with every 4 with a simple mean unless I create a more in-depth upscaling method. 
# - final choice: average every 4 segments -> this didn't lead to a uniform down slope of the water surface which means we need to use a longer segment.
# - **the water depth isn't uniform because the slope varies quite a bit between reaches**
# - Using every 20 segments wasn't any better because slope was still variable and run-time wasn't improved. Go back to every 4.

# %%

# %%
# even plotting all XS they show the same triangular shape
xs_all_avg = pd.DataFrame()
fig,ax = plt.subplots(nseg_grp, 1, figsize=(6.5, nseg_grp*2))
for n in np.arange(0, nseg_grp):
    for ns in np.arange(n*seg_len,(n+1)*seg_len):
        if ns in xs_all_df.columns:
            xs_all_df.loc[:,ns].plot(x='dist_from_right_m',y='z_m', ax=ax[n],legend=False)
    # mean of segments
    XS_grp = xs_all_df.loc[:,xs_all_df.columns.isin(np.arange(n*seg_len,(n+1)*seg_len))]
    mean_XS = pd.DataFrame(XS_grp.mean(axis=1), columns = [n])
    mean_XS.plot(ax=ax[n], color='black', linestyle='--', label='Avg')
    xs_all_avg = pd.concat((xs_all_avg, mean_XS),axis=1)

# %%
# xs_all_df

# %%
from shapely.geometry import LineString
def make_8pt(xs_all_df):
    """
    Given a dataframe of channel cross-sections simplify down to 8 points
    Input:
    xs_all_df: dataframe with an index of distance from right bank and columns of cross-section number
    Output:
    lp: geodataframe with upscale cross-section points
    """
    i = 0
    # Number of cross sections
    numxs = int(len(xs_all_df.columns))
    # i is the cross-section number
    lp = pd.DataFrame(np.linspace(1,int(numxs),int(numxs)))
    lp['geometry'] = LineString([(0,0),(0,1)])
    
    for i in np.arange(0,numxs): #numxs
        # Number of points in each cross section
        numl = np.sum(pd.notna(xs_all_df.iloc[:,i]))
        # Create empty array to fill with coordinates
        lines = np.zeros((numl,2))
        # j is the number of points in each individual cross-section
        lm = LineString(list(zip(xs_all_df.index.values, xs_all_df.iloc[:,i].values)))
        tol = 0.6
        deltol = 0.1
        count = 0
        lms = LineString(lm).simplify(tolerance = tol)
        while len(list(lms.coords))>8:
            if len(list(lms.coords)) <5:
                deltol = 0.001
            temp = lms
            lms = LineString(lm).simplify(tolerance = tol)
            tol += deltol
    #         if count drops below 8 then reduce deltol
    #         if len(list(lms.coords)) <6:
    #             lms = temp
    #             tol -= deltol
    #             deltol *= 0.5     
            count += 1
    
        print(i,':',len(list(lms.coords)),end = ' - ') #count, 
        lp.geometry.iloc[int(i)] = LineString(lms)
    # return final geometry
    return(lp)
        
# some segments will never be able to match the ideal number of points despite very fine loops
lp_all = make_8pt(xs_all_df)

# %%
# 8-point cross-sections for averaged segments
lp_avg = make_8pt(xs_all_avg)

# %%
# the thalweg data is the reach data and isn't affected in this section

# create summary of XS for creating SFR inputs
xs_wide = xs_all.pivot_table(index='dist_from_right_m',values='z_m',columns='xs_num')
thalweg_pts = xs_wide.idxmin().values.astype(int)
xs_mins = xs_all.set_index(['dist_from_right_m','xs_num']).loc[list(zip(thalweg_pts, xs_wide.columns))]
XSg_in = xs_mins.reset_index('dist_from_right_m')
XSg_in.crs = xs_all.crs

# join segment data to grid
# XSg_in = gpd.sjoin(XSg_in, grid_p, predicate='within', how='inner')
XSg_in = gpd.sjoin_nearest(XSg_in, grid_p, how='inner')
# # if multiple points in one cell take first, not a big deal since there are points every 100 m
XSg_in = XSg_in.reset_index().groupby(['row','column'], as_index=False).first()
# clean up labels
# XSg_in = XSg_in.sort_values('xs_num').reset_index(drop=True)
# XSg_in['reach_order'] = np.arange(0, len(XSg_in))
# XSg_in.crs = xs_all.crs

# identify the cross-sections for each segment (introduces duplicates)
xs_join = xs_mins.reset_index()[['xs_num', 'z_m', 'geometry']].rename(columns={'z_m':'xs_z_m'})
# assign averaged segment number to the cross-sections
xs_join['iseg'] = np.nan
xs_join.loc[xs_join.xs_num<nseg_grp*seg_len, 'iseg'] = np.repeat(np.arange(1, nseg_grp+1), seg_len)
xs_join.iseg = xs_join.iseg.ffill().astype(int).values # forward fill to end
# join reaches to cross-sections
XSg_in = gpd.sjoin_nearest(sfr_reach, xs_join, how='left').drop(columns=['index_right'])
# drop duplicates introduced
XSg_in = XSg_in.drop_duplicates(['row','column'], keep='first')
# clean up order and fill in elevation
XSg_in = XSg_in.sort_values('id')
XSg_in.xs_z_m = XSg_in.xs_z_m.interpolate('linear').copy()
XSg_in['reach_order'] = np.arange(0, len(XSg_in))


# %%
# # here the segment, reach numbers changes between cross-section usage
# # for 1 cross-section each segment
# seg_num = pd.Series(np.arange(1, XSg_in.shape[0]+1)) # add the segment that corresponds to each cross section

# # make segment numbers
# seg_num = pd.Series(np.full((XSg_in.shape[0]), np.nan))
# # segment numbers based on groups
# seg_num.iloc[:nseg_grp*seg_len] = np.repeat(np.arange(1, nseg_grp+1), seg_len)
# seg_num = seg_num.ffill().astype(int).values # forward fill to end

# %%
# create segment numbers, starting at 1 to allow for first segment defined by michigan bar criteria
# XSg_in['iseg'] = seg_num
# update xs_num to represent new cross-sections
XSg_in['xs_num_full'] = XSg_in.xs_num.copy()
XSg_in['xs_num'] = XSg_in.iseg.copy()-1

# create reach numbers
XSg_in['ireach'] = 1
for n in XSg_in.iseg.unique():
    XSg_in.loc[XSg_in.iseg==n, 'ireach'] = np.arange(1, np.sum(XSg_in.iseg==n)+1)


# %%
# filter cross sections to those that matched in the grid
# xs_all = xs_all[xs_all.xs_num.isin(XSg_in.xs_num_full)]
# xs_all_df = xs_all_df.loc[:, XSg_in.xs_num_full]

# %%
def clean_8pt(numxs, lp):
    """
    Create a clean dataframe with the 8-point cross-sections given a geodataframe geometry
    Fill in end points if needed
    
    """
    XS8pt = pd.DataFrame(np.zeros((numxs*8, 3)), columns=['xs_num','dist_from_right_m','z_m'])
    XS8pt.xs_num = np.repeat(np.arange(0,numxs), 8)
    
    # lpg = gpd.GeoDataFrame(lp[:])
    xscoords = np.zeros((8, numxs))
    filler = np.zeros(2)
    filler[:] = np.nan
    for i in np.arange(0, numxs):
        coordtemp = np.array(list(lp.geometry.iloc[i].coords))
        coordtemp = coordtemp[~np.isnan(coordtemp[:,0])]
        # if missing points add to make 8
        while len(coordtemp) < 8:
            endfill = np.copy(coordtemp[-1,:]) # take last and add new point
            endfill[0] += 1 # offset with different x
            coordtemp = np.vstack((coordtemp, endfill))
        # reset distance from right to start at 0
        coordtemp[:,0] -= coordtemp[0,0]
        XS8pt.loc[XS8pt.xs_num==i,['dist_from_right_m','z_m']] = coordtemp   
        # print(i)
    return(XS8pt)



# %%
# original, one cross-section per reach (1 segment, 1 reach)
# numxs = int(len(xs_all_df.columns))
# XS8pt = clean_8pt(numxs, lp_all)

# upscaled, multiple reaches per cross-sections
numxs = int(len(xs_all_avg.columns))
XS8pt = clean_8pt(numxs, lp_avg)



# %%
# # filter for XS in final segments
XS8pt = XS8pt.loc[XS8pt.xs_num.isin(XSg_in.xs_num)]
# XS8pt = XS8pt.loc[XS8pt.xs_num.isin(XSg_in.iseg)]
XS8pt.to_csv(join(proj_dir, 'SFR',  '8pointXS_'+model_nam+'.csv'), index = False)
XS8pt = XS8pt.set_index('xs_num')

# %%
# # even plotting all XS they show the same triangular shape
fig,ax = plt.subplots()
# for n in XS8pt.index.unique()[::10]:
for n in XS8pt.index.unique():
    XS8pt.loc[n].plot(x='dist_from_right_m',y='z_m', ax=ax,legend=False)

# %% [markdown]
# Need to look at filtering down the number of cross-sections.
# - groups of 10 seems to fit pretty well
# - also need to cut down from 200 m wide to 100 m and to have more focus on the main channel

# %%
# XSg_z = XSg_in.copy().set_index('iseg')
XSg_z = XSg_in.copy().set_index('reach_order')

# find minimum value in XS related to thalweg
# XSg_z['z_m_min'] = xs_all.dissolve('xs_num','min').z_m
XSg_z['z_m_min'] = XSg_z.xs_z_m

#roling mean of 6 window centered removes any negative slope
XSg_z['z_m_min_cln'] = XSg_z.z_m_min.rolling(6,center=False).mean()

# calculate slope and fill NAs, fill slope with nearby
z_cln_diff = XSg_z.z_m_min_cln.diff().bfill()
XSg_z['slope'] = z_cln_diff.abs()/XSg_z.rch_len
# correct slope less than 1E-4
XSg_z.loc[XSg_z.slope<1E-4,'slope'] = 1E-4
# rolling mean of slope to clean up slope for manning's
XSg_z['slope_raw'] = XSg_z.slope.copy()
# XSg_z.slope = XSg_z.slope.rolling(10, center=True, min_periods=1).mean()
XSg_z.slope = XSg_z.slope.rolling(10, center=True, min_periods=10).mean().bfill().ffill()

# fix str bot so all is downward sloping
for i in XSg_z.index[-2::-1]:
# fill NAs due to rolling mean, with backward filling
    if np.isnan(XSg_z.loc[i,'z_m_min_cln']):
        XSg_z.loc[i,'z_m_min_cln'] = XSg_z.loc[i+1,'z_m_min_cln'] + XSg_z.loc[i,'slope']*delr

for i in XSg_z.index[:-1]:
    if XSg_z.loc[i+1,'z_m_min_cln'] >= XSg_z.loc[i,'z_m_min_cln']:
        XSg_z.loc[i+1,'z_m_min_cln'] = XSg_z.loc[i,'z_m_min_cln'] - XSg_z.loc[i,'slope']*delr

# calculate the elevation if we use slope only
XSg_z['z_m_slope'] = XSg_z.z_m.max() - (XSg_z.rch_len*XSg_z.slope).cumsum()

avg_slope = (XSg_z.z_m_min_cln.max() - XSg_z.z_m_min_cln.min())/(XSg_z.rch_len.sum())

# new column for easier modflow consistency
# XSg_z['z_m_final'] = XSg_z.z_m_min_cln.copy()
XSg_z['z_m_final'] = XSg_z.z_m_slope.copy()

# %%
plt.plot(XSg_z.z_m_min_cln.diff().bfill().multiply(-1)/100 ,label='Z diff')
# plt.plot(XSg_z.slope.rolling(10, center=False, min_periods=10).mean().bfill().ffill(), label='roll mean slope')
# plt.plot(XSg_z.slope.rolling(10, center=True, min_periods=10).mean().bfill().ffill(), label='roll mean slope center')
plt.plot(XSg_z.slope, label='roll mean slope center')
plt.plot(XSg_z.slope_raw, label='slope')
plt.axhline(avg_slope,label= 'avg slope', color='black')
plt.legend()

# %% [markdown]
# Applying a rolling mean to the slope doesn't have a big impact on the elevations if we only correct those that have a downstream rise. If we apply the slope method to all segments then we see a more uniform channel with perhaps more fall then there should be by the outlet.
# - the rolling mean (window 10, centered) helped create more regularly down sloping conditions

# %%
x = XSg_z.index.values

# plt.plot(x, XSg_in.xs_z_m, label='XS Str Top')
plt.plot(x, XSg_z.z_m, label='Cell Riv Mean')
plt.plot(x, XSg_z.z_m_min_cln, label='Str Top')
plt.plot(x, XSg_z.z_m_min, label='min')

plt.plot(x, XSg_z.z_m_slope, label='slope driven')
plt.legend()

# %%
XSg_z['easting'] = XSg_z.geometry.centroid.x
XSg_z['northing'] = XSg_z.geometry.centroid.y
XSg_z.drop(columns=['geometry']).to_csv(join(proj_dir, 'SFR', 'reach_data_cleaned.csv'))

# %%
# XSg_z.to_file(join(proj_dir, 'SFR', 'reach_data_cleaned.shp'))
