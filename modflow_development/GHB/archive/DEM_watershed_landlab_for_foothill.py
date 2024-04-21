# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
# standard python utilities
import os
import glob
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
import landlab

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

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
# m_domain = m_domain.to_crs('epsg:3310')

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
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
from rasterio.transform import Affine


# %%
# define function to interpolate/krige water surface elevations
# part 1/2
def lin_krige(filename, domain_buffered, res):
    # filename of shapefiles with water surface elevation points
    # m_domain shapefile to clip grid
    # raster resolution for kriging and output
    
    # Kriging set up
    df = gpd.read_file(filename)
    df = df.to_crs('epsg:32610')
    
    df = gpd.sjoin(df,domain_buffered, op = 'within', how = 'inner')

    x = df.geometry.x.values
    y = df.geometry.y.values
    z = df.WSEL.values
    res = 100
    gridx = np.arange(np.min(x), np.max(x), res)
    gridy = np.arange(np.min(y), np.max(y), res)
    
    # Kriging
    # linear, gaussian, spherical, exponential, hole-effect and power
    OK = OrdinaryKriging(
        x,
        y,
        z,
        variogram_model="linear",
        verbose=True,
        enable_plotting=True,
    #     exact_values = False,
        enable_statistics = True,
        nlags = 50
    )

    
    # z is the kriged grid and ss is the variance grid (sigma ^2)
    z, ss = OK.execute("grid", gridx, gridy)
    # flip data because np sets 0,0 in top left while raster is bottom left
    Z  = np.flip(z.data,axis = 0)
    SS = np.flip(ss.data,axis = 0)
    
    
    # Output file creation
    new_dataset1 = rasterio.open(ghb_dir+'/interpolated_data/' +str(year)+'_kriged.tif',
                           'w',driver = 'GTiff',
                           height=Z.shape[0], width = Z.shape[1],
                           count = 1, dtype = Z.dtype,
                           crs = f2013.crs,
                           transform = transform)

    new_dataset1.write(Z,1)
    new_dataset1.close()

    new_dataset2 = rasterio.open(ghb_dir+'/interpolated_data/' +str(year)+'_variance.tif',
                               'w',driver = 'GTiff',
                               height=Z.shape[0], width = Z.shape[1],
                               count = 1, dtype = Z.dtype,
                               crs = f2013.crs,
                               transform = transform)

    new_dataset2.write(SS,1)
    new_dataset2.close()

# %%
strtyear = 2013
endyear = 2019

for year in np.arange(strtyear,endyear):
    filename = glob.glob(ghb_dir+'/Fall_'+str(year)+'_Elevation_Points/*.shp')[0]
    lin_krige(filename, domain_buffered, res=100)

# %%
# x = f2013.geometry.x.values
# y = f2013.geometry.y.values
# z = f2013.WSEL.values
# res = 100
# gridx = np.arange(np.min(x), np.max(x), res)
# gridy = np.arange(np.min(y), np.max(y), res)

# %%
# # linear, gaussian, spherical, exponential, hole-effect and power
# OK = OrdinaryKriging(
#     x,
#     y,
#     z,
#     variogram_model="linear",
#     verbose=True,
#     enable_plotting=True,
# #     exact_values = False,
#     enable_statistics = True,
#     nlags = 50
# #     anisotropy_angle = -52.9
# )
# # anisotropy_angle is 52.9 degrees cw, and 

# %%
# z is the kriged grid and ss is the variance grid (sigma ^2)
# z, ss = OK.execute("grid", gridx, gridy)

# %%
# plt.imshow(z)
# plt.colorbar()
# plt.show()

# plt.imshow(ss)
# plt.imshow(np.sqrt(ss))
# plt.colorbar()

# %%
# from rasterio.transform import Affine
# # res is defined above when the grid is defined
# transform = Affine.translation(np.min(x) - res / 2, np.max(y) - res / 2) * Affine.scale(res, -res)
# transform

# %%
# Z  = np.flip(z.data,axis = 0)
# SS = np.flip(ss.data,axis = 0)

# %%
# new_dataset1 = rasterio.open(ghb_dir+'/interpolated_data/f2013_kriged.tif',
#                            'w',driver = 'GTiff',
#                            height=Z.shape[0], width = Z.shape[1],
#                            count = 1, dtype = Z.dtype,
#                            crs = f2013.crs,
#                            transform = transform)

# new_dataset1.write(Z,1)
# new_dataset1.close()

# new_dataset2 = rasterio.open(ghb_dir+'/interpolated_data/f2013_var.tif',
#                            'w',driver = 'GTiff',
#                            height=Z.shape[0], width = Z.shape[1],
#                            count = 1, dtype = Z.dtype,
#                            crs = f2013.crs,
#                            transform = transform)

# new_dataset2.write(SS,1)
# new_dataset2.close()

# %% [markdown]
# ## Load with flopy to crop down to domain

# %%
from flopy.utils import Raster
import flopy


# %%
loadpth = gwfm_dir+'/JupyterNotebooks/SFR/data'
m = flopy.modflow.Modflow.load('MF.nam', model_ws=loadpth, 
                                exe_name='mf-NWT', version='mfNWT')


# %%
def krige_crop(filename):
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

   
    np.savetxt(ghb_dir+'/final_fallWSEL_arrays/'+str(year)+'_kriged_WSEL.tsv',
               df_grid, delimiter = '\t')

# %%
strtyear = 2013
endyear = 2019

for year in np.arange(strtyear,endyear):
    filename = glob.glob(ghb_dir+'/Fall_interpolated_data/'+str(year)+'_kriged.tif')[0]
    krige_crop(filename)

# %%
# f2013_rio = Raster.load(ghb_dir+'/interpolated_data/f2013_kriged.tif')

# %%
# # 2 seconds for linear, 0.2 seconds for nearest when 200m spacing is used
# # 9 seconds for linear when 100m grid is used for kriging
# t0 = time.time()
# f2013_grid = f2013_rio.resample_to_grid(m.modelgrid.xcellcenters,
#                                 m.modelgrid.ycellcenters,
#                                 band=f2013_rio.bands[0],
#                                 method="linear")
# resample_time = time.time() - t0
# print("Resample time, nearest neighbor: {:.3f} sec".format(time.time() - t0))

# %%
# ghb_dir = gwfm_dir+'/GHB_data'
# fig,ax= plt.subplots(figsize =(8,8))

# m_domain.plot(ax=ax,edgecolor = 'black', color = 'none')

# ax = f2013_rio.plot(ax=ax)
# cbar = plt.colorbar(ax.images[0])

# %%

# %% [markdown]
# ## Geostatus Guy - Professor Michael Pyrcz version based on Deutsch and Journel 1992

# %%
import geostatspy.GSLIB as GSLIB                          # GSLIB utilities, viz and wrapped functions
import geostatspy.geostats as geostats                    # GSLIB converted to Python


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
f2013.describe().transpose()

# %%
f2013['NWSEL'], tvWSEL, tnsWSEL = geostats.nscore(f2013, 'WSEL')

# %%
plt.subplot(121)
plt.title('NScore WSEL')
f2013.NWSEL.hist()
plt.subplot(122)
plt.title('WSEL')
f2013.WSEL.hist()
plt.show()

cmap = plt.cm.plasma                  
GSLIB.locmap_st(f2013,'x','y','NWSEL',xmin,xmax,ymin,ymax,-3,3,
                'Nscore Porosity - All Facies','X (m)','Y (m)','Nscore Porosity',cmap)

# %%
# Calculate Sample Data Isotropic Variograms
tmin = -9999.; tmax = 9999.; 
lag_dist = 100.0; lag_tol = 100.0; nlag = 50; bandh = 9999.9; azi = 0; atol = 90.0; isill = 1

lag, wsel_gamma, wsel_npair = geostats.gamv(f2013,"x","y","NWSEL",tmin,tmax,
                                                    lag_dist,lag_tol,nlag,azi,atol,bandh,isill)

# %%
plt.plot(lag,wsel_gamma,'x',color = 'black',label = 'All')


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

tmin = -9999.9; tmax = 9999.9; 


# %%
xmin,xmax,ymin,ymax,xsiz

# %%
skmean_wsel = f2013.WSEL.mean()      # simple kriging mean (used if simple kriging is selected below)
ktype = 1                                  # kriging type, 0 - simple, 1 - ordinary
radius = 10000                               # search radius for neighbouring data
nxdis = 100; nydis = 100                       # number of grid discretizations for block kriging (not tested)
ndmin = 0; ndmax = 10                      # minimum and maximum data for an estimate

wsel_vario = GSLIB.make_variogram(nug=0.0,nst=1,it1=1,cc1=1.0,azi1=0,hmaj1=100,hmin1=100) # wsel variogram

# %%
wsel_kmap, wsel_vmap = geostats.kb2d(f2013,'x','y','WSEL',tmin,tmax,nx,xmin,xsiz,ny,ymin,ysiz,nxdis,nydis,
         ndmin,ndmax,radius,ktype,skmean_wsel,wsel_vario)

plt.subplot(211)
GSLIB.locpix_st(wsel_kmap,xmin,xmax,ymin,ymax,xsiz,-40,170,f2013,'x','y','WSEL','Kriging Estimate','X(m)','Y(m)','WSEL (m)',cmap)

plt.subplot(212)
GSLIB.pixelplt_st(wsel_vmap,xmin,xmax,ymin,ymax,xsiz,0.0,1000,'Kriging Variance','X(m)','Y(m)','WSEL $(m^2)$',cmap)

# %%
