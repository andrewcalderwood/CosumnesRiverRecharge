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

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %%
doc_dir = os.getcwd()
while os.path.basename(doc_dir) != 'Documents':
    doc_dir = os.path.dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = os.path.dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
gwfm_dir
ghb_dir = gwfm_dir+'/GHB_data'

# %%
# Original model domain, 44.7 deg angle
# m_domain = gpd.read_file(gwfm_dir+'\\GWModelDomain_UTM10N\\GWModelDomain_Rec_UTM10N.shp')
# New model domain 52.9 deg
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')

# %%
# m_domain = m_domain.to_crs('epsg:3310')

# %%
fig,ax= plt.subplots(figsize =(8,8))
m_domain.plot(ax=ax,edgecolor = 'black', color = 'none')
buffer = 2e4

# buffer domain in order to include boundary observations to influence water table elevations at the model edge
domain_buffered = gpd.GeoDataFrame(data = m_domain.copy(), geometry = m_domain.buffer(buffer))
domain_buffered.plot(ax=ax,edgecolor = 'black', color = 'none')


# %%
f2013 = gpd.read_file(ghb_dir+'/Fall_Spring_GWE_Points/fall2013_GWE_Points.shp')

f2013 = f2013.to_crs('epsg:32610')

# %%
fig,ax= plt.subplots(1,1,figsize =(4,4))
m_domain.plot(ax=ax,edgecolor = 'black', color = 'none')
domain_buffered.plot(ax=ax,edgecolor = 'black', color = 'none')
# units are in feet
f2013.plot('GWE', ax=ax, legend = True)
# fig.tight_layout()

# %%
# f2013.WSEL.hist()
# scale by 100 to make all positive values
# f2013['logWSEL'] = np.log(f2013.WSEL.values+100)
# f2013.logWSEL.hist()

# %%
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
from rasterio.transform import Affine


# %%
# define function to interpolate/krige water surface elevations
# part 1/2
def lin_krige(filename, domain_buffered, res, season,year,z_name):
    # filename of shapefiles with water surface elevation points
    # m_domain shapefile to clip grid
    # raster resolution for kriging and output
    
    # Kriging set up
    df = gpd.read_file(filename)
    df = df.to_crs('epsg:32610')
    # remove NA values
    df=df.dropna(subset=[z_name])
    
    df = gpd.sjoin(df,domain_buffered, predicate = 'within', how = 'inner')

    x_in = df.geometry.x.values
    y_in = df.geometry.y.values
    z_in = df[z_name].values
    res = 100
    gridx = np.arange(np.min(x_in), np.max(x_in), res)
    gridy = np.arange(np.min(y_in), np.max(y_in), res)
    
    # Kriging
    # linear, gaussian, spherical, exponential, hole-effect and power
    OK = OrdinaryKriging(
        x_in,
        y_in,
        z_in,
        # gaussian overweights low values causing all data to look the same, power looks okay with high lag
        # linear still seems best if I can avoid the singular matrix issue
        variogram_model="linear", 
        verbose=True,
        enable_plotting=True,
    #     exact_values = False,
        enable_statistics = True,
        nlags = 50, # if lags is too low, then higher values seem to dominate?
        pseudo_inv=True
    )

    # z is the kriged grid and ss is the variance grid (sigma ^2)
    z, ss = OK.execute("grid", gridx, gridy)
    # flip data because np sets 0,0 in top left while raster is bottom left
    Z  = np.flip(z.data,axis = 0)
    SS = np.flip(ss.data,axis = 0)

    transform = Affine.translation(np.min(x_in) - res / 2, np.max(y_in) - res / 2) * Affine.scale(res, -res)
    
    # Output file creation
    new_dataset1 = rasterio.open(ghb_dir+'/interpolated_data/'+season +str(year)+'_kriged.tif',
                           'w',driver = 'GTiff',
                           height=Z.shape[0], width = Z.shape[1],
                           count = 1, dtype = Z.dtype,
                           crs = df.crs,
                           nodata = -9999,
                           transform = transform)

    new_dataset1.write(Z,1)
    new_dataset1.close()

    new_dataset2 = rasterio.open(ghb_dir+'/interpolated_data/'+season+str(year)+'_variance.tif',
                               'w',driver = 'GTiff',
                               height=Z.shape[0], width = Z.shape[1],
                               count = 1, dtype = Z.dtype,
                               crs = df.crs,
                                nodata = -9999,
                               transform = transform)

    new_dataset2.write(SS,1)
    new_dataset2.close()


# %% [markdown]
# ## Kriging with shapefiles from DWR periodic database
# **The units of GWE are in feet**

# %%

def load_gwe_pts(y,s):
    name = s+str(y)
    fn = ghb_dir+'/Fall_Spring_GWE_Points/'+name+'_GWE_Points.shp'
    df = gpd.read_file(fn)
    df = df.to_crs('epsg:32610')
    df = df.dropna(subset=['GWE'])
#     df = df.loc[:,]
    return(df)

# spring2020 = load_gwe_pts(2020,'spring')
# spring2020.plot('GWE')

# fall2019 = load_gwe_pts(2019,'fall')
# fall2019.plot('GWE')

# diff = fall2019.loc[:,['SITE_CODE','GWE','month']].join(spring2020.set_index('SITE_CODE'),
#                                                             on='SITE_CODE',rsuffix='s20')
# diff['gwe_diff'] = diff.GWEs20 - diff.GWE
# # diff.gwe_diff.plot()
# # diff.plot('gwe_diff',legend=True)

# diff.gwe_diff.mean()


# %%
for y in np.arange(2010,2022):
# for y in np.arange(2000,2001):
    for s in ['spring','fall']:
        df = load_gwe_pts(y,s)
        print(s+' '+str(y)+' mean GWE (ft) '+ str(df.GWE.mean()))


# %%
s='spring'
y=2000

name = s+str(y)
fn = ghb_dir+'/Fall_Spring_GWE_Points/'+name+'_GWE_Points.shp'
        
df = gpd.read_file(fn)
df = df.to_crs('epsg:32610')
df = df.dropna(subset=['GWE'])
df = gpd.sjoin(df,domain_buffered, op = 'within', how = 'inner')
x_in = df.geometry.x.values
y_in = df.geometry.y.values
z_in = df['GWE'].values
res = 100
gridx = np.arange(np.min(x_in), np.max(x_in), res)
gridy = np.arange(np.min(y_in), np.max(y_in), res)

# Kriging
# linear, gaussian, spherical, exponential, hole-effect and power
OK = OrdinaryKriging(
    x_in,
    y_in,
    z_in,
    variogram_model="linear", # use linear with high number of lags to cause more averaging
    verbose=True,
    enable_plotting=True,
#     exact_values = False,
    enable_statistics = True,
    # 20*100 = 2,000 m lag distance = fair for typical well spacing, 50 is too much for some years (not enough data)
    nlags = 10,
    pseudo_inv=True # can lead to more stability but takes longer
)
# 2014 requires nlags = 5, 2012 nlags=30, with exponential you can go farther
# z is the kriged grid and ss is the variance grid (sigma ^2)
z, ss = OK.execute("grid", gridx, gridy)
# flip data because np sets 0,0 in top left while raster is bottom left
Z_2019  = np.flip(z.data,axis = 0)
SS = np.flip(ss.data,axis = 0)



# %%
plt.imshow(Z_2019)
plt.colorbar()
plt.show()
# plt.imshow(Z_2020)

# %%
plt.plot(Z_2019[0,:],label='fall 2019')
# plt.plot(Z_2020[0,:])
# plt.plot(Z_2019[99,:])
# plt.plot(Z_2020[99,:])
plt.legend()

# %%
# Z_gauss, Z_exp, Z_power, Z_linear # cR =  757, 522, 1053, DNE

# %%
# strt_year = 2000
# end_year = 2010
strt_year = 2010
end_year = 2022
# linear with 50 lags works from 2013:2021
for s in ['spring','fall']:
    for y in np.arange(strt_year,end_year):
#     for y in np.arange(2000,2001):
        name = s+str(y)
        fn = ghb_dir+'/Fall_Spring_GWE_Points/'+name+'_GWE_Points.shp'
        lin_krige(fn, domain_buffered, res=100, season=s,year=y, z_name='GWE')
        print(y)


# %% [markdown]
# # Crop kriged data to array and at distance of 1000 m
# Using zonal_stats functions sample the entire grid for the mean WSEL for each cell and output as an array with previous format. Additionally use zonal_stats for the bufferline 500m, 1,000m from the grid edge.

# %%
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')

m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')


# %%
from rasterstats import zonal_stats


# %%
s='fall'
y = 2009
name = s+str(y)

fn = ghb_dir+'/interpolated_data/'+name+'_kriged.tif'

# takes several minutes
zs = zonal_stats(grid_p, raster=fn, stats=['min', 'max', 'mean', 'std'])
# convert to dataframe
zs_df = pd.DataFrame(zs)
# join zone stats of DEM to parcel data
zs_df = grid_p.join(zs_df)

# %%
zs_df.plot('mean')

# %% [markdown]
# ## Load with flopy to crop down to domain

# %%
flopy_dir = doc_dir+'/GitHub/flopy/'
if flopy_dir not in sys.path:
    sys.path.append(flopy_dir)
# sys.path
import flopy 


from flopy.utils import Raster


# %%
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
def krige_crop(filename, season,year):
    # krige_crop, reduces the geotiff to the model domain by resampling
    # additionally it filters out the data needed just for the NW,SE BCs
    df_rio = Raster.load(filename)
    
    # 2 seconds for linear, 0.2 seconds for nearest when 200m spacing is used
    # 9 seconds for linear when 100m grid is used for kriging
    # t0 = time.time()
    df_grid = df_rio.resample_to_grid(modelgrid,
                                    band=df_rio.bands[0],
                                    method="linear")
    # resample_time = time.time() - t0
    # print("Resample time, nearest neighbor: {:.3f} sec".format(time.time() - t0))

   
    np.savetxt(ghb_dir+'/final_WSEL_arrays/'+season+str(year)+'_kriged_WSEL.tsv',
               df_grid, delimiter = '\t')

# %%
strt_year

# %%
for s in ['spring','fall']:
    for y in np.arange(strt_year,end_year):
        name = s+str(y)
        fn = ghb_dir+'/interpolated_data/'+name+'_kriged.tif'
        krige_crop(fn, season=s,year=y)

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
