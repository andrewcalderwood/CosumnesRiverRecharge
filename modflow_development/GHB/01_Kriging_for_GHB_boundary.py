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
from os.path import join, basename, dirname, expanduser
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
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')

# dir of all gwfm data
gwfm_dir = os.path.dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
gwfm_dir
ghb_dir = gwfm_dir+'/GHB_data'

# %%
# Original model domain, 44.7 deg angle
# m_domain = gpd.read_file(gwfm_dir+'\\GWModelDomain_UTM10N\\GWModelDomain_Rec_UTM10N.shp')
# New model domain 52.9 deg
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')

buffer = 2e4

# buffer domain in order to include boundary observations to influence water table elevations at the model edge
domain_buffered = gpd.GeoDataFrame(data = m_domain.copy(), geometry = m_domain.buffer(buffer))

# %% [markdown]
# # Specify start and end year
# These could be set up to auto update based on spring, fall points available

# %%
# strt_year = 2000
# end_year = 2010
strt_year = 2022
end_year = 2024

# %%
fn = ghb_dir+'/Fall_Spring_GWE_Points/*shp'
pt_files = pd.DataFrame([basename(f) for f in glob.glob(fn)], columns=['fn'])
pt_files['year'] = pt_files.fn.str.extract(r'(\d{4})').astype(int)
pt_files.year.min(), pt_files.year.max()

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
z_name = 'gwe'


# %%

def load_gwe_pts(y,s, z_name):
    name = s+str(y)
    fn = ghb_dir+'/Fall_Spring_GWE_Points/'+name+'_GWE_Points.shp'
    df = gpd.read_file(fn)
    df = df.to_crs('epsg:32610')
    df = df.dropna(subset=[z_name])
#     df = df.loc[:,]
    return(df)



# %%
for y in np.arange(strt_year,end_year):
# for y in np.arange(2000,2001):
    for s in ['spring','fall']:
        df = load_gwe_pts(y,s, z_name)
        print(s+' '+str(y)+' mean GWE (ft) '+ str(df[z_name].mean()))


# %%
# Z_gauss, Z_exp, Z_power, Z_linear # cR =  757, 522, 1053, DNE

# %%

# linear with 50 lags works from 2013:2021
for s in ['spring','fall']:
    for y in np.arange(strt_year,end_year):
#     for y in np.arange(2000,2001):
        name = s+str(y)
        fn = ghb_dir+'/Fall_Spring_GWE_Points/'+name+'_GWE_Points.shp'
        lin_krige(fn, domain_buffered, res=100, season=s,year=y, z_name=z_name)
        print(y)

