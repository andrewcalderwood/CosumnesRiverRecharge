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

# %%
# standard python utilities
import os
import sys
import pandas as pd
import numpy as np
import calendar
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import pyproj # for converting proj4string
import shapely
import geopandas as gpd
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %%
import geostatspy.GSLIB as GSLIB                          # GSLIB utilities, viz and wrapped functions
import geostatspy.geostats as geostats                    # GSLIB converted to Python


# %%
## Set up directory referencing
# Package data
gwfm_dir = os.path.dirname(os.path.dirname(os.getcwd()))
ghb_dir = gwfm_dir+'/GHB_data'

# %%
# Original model domain, 44.7 deg angle
# m_domain = gpd.read_file(gwfm_dir+'\\GWModelDomain_UTM10N\\GWModelDomain_Rec_UTM10N.shp')
# New model domain 52.9 deg
m_domain = gpd.read_file(gwfm_dir+'\\NewModelDomain\\GWModelDomain_52_9deg_UTM10N_WGS84.shp')

# %%
f2013 = gpd.read_file(ghb_dir+'/Fall_2013_Elevation_Points/F2013_WSEL_Points_20150720_090332.shp')
f2013 = f2013.to_crs('epsg:32610')

# %%
fig,ax= plt.subplots(figsize =(8,8))
m_domain.plot(ax=ax,edgecolor = 'black', color = 'none')
buffer = 2e4
m_domain.buffer(buffer).plot(ax=ax,edgecolor = 'black', color = 'none')

# buffer domain in order to include boundary observations to influence water table elevations at the model edge
domain_buffered = gpd.GeoDataFrame(data = m_domain, geometry = m_domain.buffer(buffer))

# %%
fig,ax= plt.subplots(figsize =(8,8))
m_domain.plot(ax=ax,edgecolor = 'black', color = 'none')

f2013 = gpd.sjoin(f2013,domain_buffered, op = 'within', how = 'inner')
f2013.plot('WSEL', ax=ax, legend = True)

# %%
# f2013.WSEL.hist()
# scale by 100 to make all positive values
f2013['logWSEL'] = np.log(f2013.WSEL.values+100)
f2013.logWSEL.hist()

# %%
import scipy.spatial as sp  # for fast nearest neighbor search

def invdist(
    df,
    xcol,
    ycol,
    vcol,
    tmin,
    tmax,
    nx,
    xmn,
    xsiz,
    ny,
    ymn,
    ysiz,
    ndmin,
    ndmax,
    radius,
    power
):
    """ Inverse Distance to Python by Michael Pyrcz, the University of Texas at
    Austin (April, 2020).  Based on modification of the GSLIB kb2d program by Deutsch and Journel (1997)
    :param df: pandas DataFrame with the spatial data
    :param xcol: name of the x coordinate column
    :param ycol: name of the y coordinate column
    :param vcol: name of the property column
    :param tmin: property trimming limit
    :param tmax: property trimming limit
    :param nx: definition of the grid system (x axis)
    :param xmn: definition of the grid system (x axis)
    :param xsiz: definition of the grid system (x axis)
    :param ny: definition of the grid system (y axis)
    :param ymn: definition of the grid system (y axis)
    :param ysiz: definition of the grid system (y axis)
    :param ndmin: minimum number of data points to use for kriging a block
    :param ndmax: maximum number of data points to use for kriging a block
    :param radius: maximum isotropic search radius
    :param power: the inverse distance power
    :return:
    """
    
# Constants
    UNEST = -999.
    EPSLON = 1.0e-10
    VERSION = 0.1
        
# Load the data
    df_extract = df.loc[(df[vcol] >= tmin) & (df[vcol] <= tmax)]    # trim values outside tmin and tmax
    nd = len(df_extract)
    ndmax = min(ndmax,nd)
    x = df_extract[xcol].values
    y = df_extract[ycol].values
    vr = df_extract[vcol].values
    
# Allocate the needed memory:   
    xa = np.zeros(ndmax)
    ya = np.zeros(ndmax)
    vra = np.zeros(ndmax)
    dist = np.zeros(ndmax)
    nums = np.zeros(ndmax)
    s = np.zeros(ndmax)
    estmap = np.zeros((nx,ny))
    
# Make a KDTree for fast search of nearest neighbours   
    dp = list((y[i], x[i]) for i in range(0,nd))
    data_locs = np.column_stack((y,x))
    tree = sp.cKDTree(data_locs, leafsize=16, compact_nodes=True, copy_data=False, balanced_tree=True)

# Summary statistics for the data after trimming
    avg = vr.mean()
    stdev = vr.std()
    ss = stdev**2.0
    vrmin = vr.min()
    vrmax = vr.max()

# Initialize accumulators:
    rad2 = radius*radius

# MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
    nk = 0
    ak = 0.0
    vk = 0.0
    for iy in range(0,ny):
        yloc = ymn + (iy-0)*ysiz  
        for ix in range(0,nx):
            xloc = xmn + (ix-0)*xsiz
            current_node = (yloc,xloc)
        
# Find the nearest samples within each octant: First initialize
# the counter arrays:
            na = -1   # accounting for 0 as first index
            dist.fill(1.0e+20)
            nums.fill(-1)
            dist, nums = tree.query(current_node,ndmax) # use kd tree for fast nearest data search
            # remove any data outside search radius
            nums = nums[dist<radius]
            dist = dist[dist<radius] 
            nd = len(dist)        
            
# Is there enough samples?
            if nd < ndmin:   # accounting for min index of 0
                est  = UNEST
#                print('UNEST at ' + str(ix) + ',' + str(iy))  # option to include this error
            else:

# Put coordinates and values of neighborhood samples into xa,ya,vra:
                for ia in range(0,nd):
                    jj = int(nums[ia])
                    xa[ia]  = x[jj]
                    ya[ia]  = y[jj]
                    vra[ia] = vr[jj]
                    
# Solve for weights
                dist = np.sqrt((xa-xloc)*(xa-xloc) + (ya-yloc)*(ya-yloc)) 
                s = 1/((dist + EPSLON)**power)        # calculate inverse weights
                s = s / np.sum(s)             # constrain sum of the weights to 1.0 for unbiasedness
                est = 0.0                
                for ia in range(0,nd):
                    est = est + s[ia] * vra[ia]
                # this line causes issues if ny!=nx
            estmap[ny-iy-1,ix] = est

# Track the estimates
            if est > UNEST:
                nk = nk + 1
                ak = ak + est
                vk = vk + est*est

# END OF MAIN LOOP OVER ALL THE BLOCKS:

    if nk >= 1:
        ak = ak / float(nk)
        vk = vk/ float(nk) - ak*ak
        print('  Estimated   ' + str(nk) + ' blocks ')
        print('      average   ' + str(ak) + '  variance  ' + str(vk))

    return estmap



# %%
f2013['x'] = f2013.geometry.x
x = f2013.x
f2013['y'] = f2013.geometry.y
y = f2013.y

z = f2013.WSEL.values
gridx = np.arange(np.min(x), np.max(x), 200)
gridy = np.arange(np.min(y), np.max(y), 200)

# %%
xmin = f2013.x.min()
xmax = f2013.x.max()
ymin = f2013.y.min()
ymax = f2013.y.max()

cmap = plt.cm.plasma                    # color map
gw_min = f2013.WSEL.min()
gw_max = f2013.WSEL.max()

# %%
xlen = xmax-xmin
ylen = ymax-ymin
xsiz = 100
ysiz = 100
nx = xlen/xsiz
ny = ylen/ysiz
nx = np.ceil(nx).astype(int)
ny = np.ceil(ny).astype(int)
if nx>ny:
    ny = nx
else:
    nx=ny
# readjust xmax, ymax based on nx, ny
xmax = xmin + xsiz*nx
ymax = ymin + ysiz*ny

# %%
xmin, xmax, ymin, ymax

# %%
tmin = -9999.9; tmax = 9999.9; 
# nx = 100; xsiz = 10.0; xmn = mon_gpd.Easting.min()
# ny = 100; ysiz = 100.0; ymn = mon_gpd.Northing.min()

# nx = 100; xsiz = 200.0; xmn = f2013.x.min(); xmin = xmn - 0.1 * xsiz; xmax = xmin + nx * xsiz
# ny = 100; ysiz = 200.0; ymn = f2013.y.min(); ymin = ymn - 0.1 * ysiz; ymax = ymin + ny * ysiz

ndmin = 1; ndmax = 10
radius = 10000.0; power = 2

invdist_est = invdist(f2013,'x','y','WSEL',tmin,tmax,nx,xmin,xsiz,ny,ymin,ysiz,ndmin,ndmax,radius,power)
print('The single estimate at x = ' + str(xmin) + ', y = ' + str(ymin) + ' is ' + str(invdist_est[0][0]) + '.')

# %%
GSLIB.locpix_st(invdist_est,xmin,xmax,ymin,ymax,xsiz,gw_min,gw_max,
                f2013,'x','y','WSEL','Power = 2','Easting (m)','Northing (m)','Feature',cmap)


# %%
Z = invdist_est

# %%
from rasterio.transform import Affine
res = 100
transform = Affine.translation(np.min(x) - res / 2, np.max(y) - res / 2) * Affine.scale(res, -res)
transform

# %%
new_dataset = rasterio.open(ghb_dir+'/interpolated_data/f2013_idw.tif',
                           'w',driver = 'GTiff',
                           height=Z.shape[0], width = Z.shape[1],
                           count = 1, dtype = Z.dtype,
                           crs = f2013.crs,
                            nodatavals = -999.,
                           transform = transform)

new_dataset.write(Z,1)
new_dataset.close()

# %% [markdown]
# ## Crop and prepare for use in modflow

# %%
from flopy.utils import Raster
import flopy


# %%
loadpth = gwfm_dir+'/JupyterNotebooks/SFR/data'
m = flopy.modflow.Modflow.load('MF.nam', model_ws=loadpth, 
                                exe_name='mf-NWT', version='mfNWT')


# %%

# %%
def idw_crop(filename):
    # krige_crop, reduces the geotiff to the model domain by resampling
    # additionally it filters out the data needed just for the NW,SE BCs
    df_rio = Raster.load(filename)
    
    # 2 seconds for linear, 0.2 seconds for nearest when 200m spacing is used
    # 9 seconds for linear when 100m grid is used for kriging
    # t0 = time.time()
    df_grid = df_rio.resample_to_grid(m.modelgrid.xcellcenters,
                                    m.modelgrid.ycellcenters,
                                    band=df_rio.bands[0],
                                    method="linear")
    # resample_time = time.time() - t0
    # print("Resample time, nearest neighbor: {:.3f} sec".format(time.time() - t0))

   
    np.savetxt(ghb_dir+'/final_fallWSEL_arrays/2013_idw_WSEL.tsv',
               df_grid, delimiter = '\t')

# %%
filename = ghb_dir+'/Fall_interpolated_data/f2013_idw.tif'
idw_crop(filename)

# %%
df_grid = df_rio.resample_to_grid(m.modelgrid.xcellcenters,
                                    m.modelgrid.ycellcenters,
                                    band=df_rio.bands[0],
                                    method="linear")

# %%
df_grid

# %%
fig, ax = plt.subplots(figsize=(9,9))

levels = np.arange(-50, 200, 5)

ax=df_rio.plot(ax=ax)
cbar = plt.colorbar(ax.images[0], values = np.arange(-70,150), cmap = 'viridis')

