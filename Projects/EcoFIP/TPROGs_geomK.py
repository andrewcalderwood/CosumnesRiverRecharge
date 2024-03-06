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
git_dir =join(doc_dir,'GitHub') 
py_dir = join(git_dir,'CosumnesRiverRecharge/python_utilities')
add_path(py_dir)
add_path(join(git_dir, 'flopy'))

from mf_utility import get_layer_from_elev
from report_cln import base_round

add_path(join(git_dir, 'CosumnesRiverRecharge','tprogs_utilities'))
import tprogs_cleaning
reload(tprogs_cleaning)
from tprogs_cleaning import tprogs_cut_elev, tprogs_arr_cut_elev, int_to_param, get_tprogs_for_elev

# %% [markdown]
# # Load data

# %%
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')

# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_elev = gpd.read_file(join(gwfm_dir,'DIS_data','grid_elevation_m_statistics.shp'))

# %%
nrow,ncol = (100,230)

# %%
grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')

# %%
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')

# %% [markdown]
# ## water level data

# %%
grid_id= 'parcel'
season='spring'
gdf_gw_long = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 
                         'GW_elevations_long_'+grid_id+'_'+season+'.shp'))
gdf_gw_long.year = gdf_gw_long.year.astype(int)

# %%
ghb_dir = join(gwfm_dir, 'GHB_data')
season='spring'
# year=2020
# load water surface data for recent years
# units are in feet
wse_yrs = np.arange(2000,2021)
wse_arr = np.zeros((len(wse_yrs), nrow,ncol))
for n,year in enumerate(wse_yrs):
    wse_arr[n] = np.loadtxt(join(ghb_dir, 'final_WSEL_arrays', season+str(year)+'_kriged_WSEL.tsv'), delimiter='\t')
# covert from feet to meters
wse_arr *= 0.3048
# correct where above land surface to 1 m below
for n,year in enumerate(wse_yrs):
    wse_arr[n, wse_arr[n] >dem_data] = dem_data[wse_arr[n]>dem_data] - 1 


# %% [markdown]
# # TPROGs Cosumnes

# %%
tprogs_name = 'tprogs_final'

mf_tprogs_dir = join(gwfm_dir,'UPW_data', tprogs_name)

tprogs_files = glob.glob(mf_tprogs_dir+'*')


# %%
from scipy.stats import gmean

tprogs_info = [80, -80, 320]

# Fleckenstein
params = pd.read_csv(upw_dir+'/ZonePropertiesInitial.csv',index_col='Zone')
# convert from m/s to m/d
params['K_m_d'] = params.K_m_s * 86400   

# params = pd.read_csv(upw_dir+'/ZonePropertiesInitial.csv',index_col='Zone')

# Maples had slightly different params that might be worth testing
# lower sand/gravel (360/120 to 60/40 m/day) and much smaller mud (0.5 to 0.006 m/day)
params_new = pd.read_csv(upw_dir+'/ZonePropertiesInitial_Maples.csv',index_col='Zone')



# %%
# grp_keys

# %%
import h5py
tprogs_fn = join(upw_dir, tprogs_name+'.hdf5')

# loads very quickly now! 
with h5py.File(tprogs_fn, mode='r') as f:
    grp = f['tprogs']
    grp_keys = f.keys()
    dset = grp['r001']
    arr = dset[:]

# %% [markdown]
# Switching from loading the tprogs line to using an h5py array sped up the load time from 30 min to a 1.5 minutes.

# %%
soil_K = np.zeros((100,nrow,ncol))
soil_K_new = np.zeros((100,nrow,ncol))
# a thicker soil layer reduces recharge
# from 1 to 2 to 3 m showed about 30-50% reduction each time
# using 2 meters would be more conservative than 1 m
soil_thick = 20
m_top = np.full((nrow,ncol), 80)

top = dem_data
bot_arr = dem_data - soil_thick
fn = join(proj_dir, 'geology','tprogs_geomK_'+str(soil_thick)+'m_depth.csv')

# takes about 30 minutes to run through all 100 realizations
# since loading tprogs data is the slowest, it would be efficient to pre-process it into hydraulic conductivity arrays
if not os.path.exists(fn):
    tic_all = time.time()
    for r in np.arange(0,100):
        tic = time.time()
        # tprogs_line = np.loadtxt(tprogs_files[r])
        # masked_tprogs= tprogs_cut_elev(tprogs_line, m_top, tprogs_info)
        # safer to open and close file to avoid issues with crashes
        with h5py.File(tprogs_fn, mode='r') as f:
            tprogs_arr = f['tprogs']['r'+str(r).zfill(3)][:].astype(float)
        masked_tprogs= tprogs_arr_cut_elev(tprogs_arr, m_top, tprogs_info)
        # these are rates in m/d
        K, Sy, Ss= int_to_param(masked_tprogs, params)
        soil_tprogs = get_tprogs_for_elev(K, top, bot_arr, tprogs_info)
        soil_K[r,:,:] = gmean(soil_tprogs,axis=0) 
        K, Sy, Ss= int_to_param(masked_tprogs, params_new)
        soil_tprogs = get_tprogs_for_elev(K, top, bot_arr, tprogs_info)
        soil_K_new[r,:,:] = gmean(soil_tprogs,axis=0) 
        toc = time.time()
        print(r, end=' ')
    #     print('Realization', t,'done in %.2f sec' %(toc-tic), end=' ')
    toc_all = time.time()
    print('\nTotal time is %.2f minutes' %((toc_all-tic_all)/60))

    soil_K_out = np.reshape(soil_K, (100*nrow, ncol))
    np.savetxt(fn, soil_K_out, delimiter=',')
else:
    soil_K_out = np.loadtxt(fn, delimiter=',')
    soil_K = np.reshape(soil_K_out, (100, nrow, ncol))

# %%
# os.remove(fn)

# %%
r=0
soil_thick = 10
m_top = np.full((nrow,ncol), 80)
top = dem_data
bot_arr = dem_data - soil_thick

with h5py.File(tprogs_fn, mode='r') as f:
    tprogs_arr = f['tprogs']['r'+str(r).zfill(3)][:].astype(float)
masked_tprogs= tprogs_arr_cut_elev(tprogs_arr, m_top, tprogs_info)# these are rates in m/d
K, Sy, Ss= int_to_param(masked_tprogs, params)
soil_tprogs = get_tprogs_for_elev(K, top, bot_arr, tprogs_info)
soil_tprogs = gmean(soil_tprogs,axis=0) 
K, Sy, Ss= int_to_param(masked_tprogs, params_new)
soil_tprogs_new = get_tprogs_for_elev(K, top, bot_arr, tprogs_info)
soil_tprogs_new = gmean(soil_tprogs_new,axis=0) 

# %% [markdown]
# A few nan values in the matrix tprogs is below land surface. Using either parameter set does not appear to effect the spatial distribution of the results, although there is a sharper contrast between low and high K

# %%

fig,ax = plt.subplots(1,2, figsize=(6.5,4),dpi=300)
im = ax[0].imshow(soil_tprogs, norm=mpl.colors.LogNorm())
fig.colorbar(im, ax=ax[0], shrink=0.3)
im = ax[1].imshow(soil_tprogs_new, norm=mpl.colors.LogNorm())
fig.colorbar(im, ax=ax[1], shrink=0.3)


# %%
r = 0
percentiles=[0,25,50,75,100]
p = np.nanpercentile(soil_tprogs, percentiles)
p_new = np.nanpercentile(soil_tprogs_new, percentiles)
for n in np.arange(0,len(percentiles)):
    print('%.2E' %p[n], '%.2E' %p_new[n])

# %% [markdown]
# Switching from Fleckenstein to Maples parameters does significantly down shift the hydraulic conductivity. Most noticeably in moderate to low K. The model would still need an additional vertical anisotropy because a rate of 67.5 m/day is much too high (could scale with soils data).
# - we need to use the Maples et al. 2019 parameters to use the geologic proxy parameter regression to get recharge potential in cm/day
# - 

# %%
gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'parcels_elevation.shp'))
# gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'sq_10ac_elevation.shp'))

# simpler geodataframe to bring dataframe to geodataframe
gdf_id = gdf_elev[['Region','geometry']].copy()

# %%
# calculate mean water surface elevation (2011-2020)
wse_mean = np.nanmean(wse_arr[11:], axis=0)

# %%
r=0
m_top = np.full((nrow,ncol), 80)
top = dem_data
# soil_thick = 10
# bot_arr = dem_data - soil_thick
bot_arr = wse_mean

with h5py.File(tprogs_fn, mode='r') as f:
    tprogs_arr = f['tprogs']['r'+str(r).zfill(3)][:].astype(float)
masked_tprogs= tprogs_arr_cut_elev(tprogs_arr, m_top, tprogs_info)# these are rates in m/d
# use new paramter set
K, Sy, Ss= int_to_param(masked_tprogs, params_new)
# sample tprogs between land surface and the water table
soil_tprogs_depth = get_tprogs_for_elev(K, top, bot_arr, tprogs_info)
soil_tprogs_new = gmean(soil_tprogs_depth,axis=0) 

# %%
# plt.imshow(soil_tprogs_depth[:,30])

# %%
plt.imshow(soil_tprogs_new, norm = mpl.colors.LogNorm())
plt.colorbar(shrink=0.5, label='$K_{geom}$ (m/d)')


# %%
## join the tprogs data to the grid shapefile to spatial join with the EcoFIP analysis units
# convert to shapefile
tprogs_gdf = grid_p.copy()
tprogs_gdf['K_m_d'] = soil_tprogs_new[tprogs_gdf.row-1, tprogs_gdf.column-1]
tprogs_gdf['wse_m_arr'] = wse_mean[tprogs_gdf.row-1, tprogs_gdf.column-1]
# spatial join to EcoFIP
tprogs_grid = gpd.overlay(tprogs_gdf, gdf_id)

# %%
# need to take the average vertical geometric mean for the analysis unit
# tprogs_K = tprogs_grid.groupby('Region').mean(numeric_only=True)[['K_m_d','wse_m_arr']].reset_index()
# average may underestimate
tprogs_K = tprogs_grid.groupby('Region').max(numeric_only=True)[['K_m_d','wse_m_arr']].reset_index()

# %%
# sample the recent period for the analysis
gdf_gw = gdf_gw_long[gdf_gw_long.year>2010].dropna(subset='dtw_m').groupby('Region').mean(numeric_only=True)
gdf_gw = gdf_id.merge(gdf_gw.reset_index())


# %%
# preliminary test of the recharge potential based on Steve's work
rech_est = gdf_gw.merge(tprogs_K)

rech_est['geomK_dtw'] = rech_est.K_m_d*rech_est.dtw_m

# 30 day average recharge estimate (cm/day)
rech_est['rch_cm_d'] = rech_est.geomK_dtw*0.0376+5.288

# %%
print('Std Dev between raster GWE to parcel and array GWE to parcel is %.2f m' %(rech_est.gwe_m - rech_est.wse_m_arr).std())

# %%
hr_area = rech_est[rech_est.rch_cm_d>25].geometry.area.sum()
lr_area = rech_est[rech_est.rch_cm_d<10].geometry.area.sum()
tot_area = rech_est.geometry.area.sum()
print('High recharge (>25 cm/day) covers %.1f %% of the estimated area' %(100*hr_area/tot_area))
print('Low recharge (<10 cm/day) covers %.1f %% of the estimated area' %(100*lr_area/tot_area))
print('Steve noted only 6% of area would have > 25 cm/day and 84% would be <10 cm/day')

# %%
fig,ax_n = plt.subplots(figsize=(6.5,5.5),dpi=300)
# log scale doesn't improve things much
rech_est.plot('rch_cm_d', ax = ax_n, legend=True,
             # norm = mpl.colors.LogNorm(vmin = rech_est.rch_cm_d.min(), vmax = rech_est.rch_cm_d.max())
              legend_kwds = {'shrink':0.7,'label':'30-day average recharge (cm/day)'}
             )

ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
# gdf_bnds(gdf_id, ax=ax_n, buf=2E3)

ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
                crs='epsg:26910', attribution=False)

# %% [markdown]
# If we use the average of the $ K_{geom} $ for a parcel then we don't see any parcels with recharge rates >25 cm/d. If we switch to the maximum of the rates then we see more distinction visually, 1.4% >25 cm/d and 94.3% < 10 cm/d instead of 99.8%
# - likely a result of the smaller high K formations in the Cosumnes that are predicted by TPROGs so the muds tend to overwhelm them.
#
# **It appears that the Maples et al. 2020 recharge regression may work but simply have lower high recharge proportions because there is not a direct IVF analog in the TProGS realizations. The realizations may not be useful because of their randomness.**
