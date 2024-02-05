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
from rasterstats import zonal_stats

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
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

import map_cln
reload(map_cln)

from mf_utility import get_layer_from_elev
from report_cln import base_round
from map_cln import plt_cln

# %% [markdown]
# # Load data

# %%
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')

# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_elev = gpd.read_file(join(gwfm_dir,'DIS_data','grid_elevation_m_statistics.shp'))

# %%
grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')

# %%
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
rivers = gpd.read_file(join(sfr_dir,"Sac_valley_rivers/Sac_valley_rivers.shp"))
rivers = rivers.to_crs('EPSG:32610')

rivers_clip = gpd.clip(rivers, m_domain)
rivers_clip = rivers_clip[['GNIS_ID','GNIS_Name','LengthKM', 'geometry']]
cr = rivers_clip[rivers_clip.GNIS_Name=='Cosumnes River']
mr = rivers_clip[rivers_clip.GNIS_Name=='Mokelumne River']


# %% [markdown]
# ## ecofip reference

# %%
# ecofip parcel scale analysis
ecofip_grid_name = 'accumAvg_Hist_regionRanks-g0'
sq_grid = gpd.read_file(join(gis_dir, ecofip_grid_name+'_10AcGrid',ecofip_grid_name+'.shp'))
parcel_grid = gpd.read_file(join(gis_dir, ecofip_grid_name+'_Parcels',ecofip_grid_name+'.shp'))
grid_cols = ['Region','FullArea','AnlysArea','RchIds','geometry']
# simplify to parcel ids only
sq_id = sq_grid[grid_cols]
parcel_id = parcel_grid[grid_cols]


# %%
def elev_stats(raster_name, gdf):
    # takes several minutes
    zs_parcels = zonal_stats(gdf, raster=raster_name, stats=['min', 'max', 'mean', 'std'])
    # convert to dataframe
    zs_df = pd.DataFrame(zs_parcels)
    # join zone stats of DEM to parcel data
    zs_df = gdf.join(zs_df)
    return zs_df
    


# %%
dem_raster_name = gwfm_dir+"/DEM_data/USGS_ten_meter_dem/modeldomain_10m_transformed.tif"
parcel_elev = elev_stats(dem_raster_name, parcel_id.to_crs(raster_crs))
parcel_elev.to_file(join(gis_dir, 'analysis_unit_reference', 'parcels_elevation.shp'))

sq_elev = elev_stats(dem_raster_name, sq_id.to_crs(raster_crs))
sq_elev.to_file(join(gis_dir, 'analysis_unit_reference', 'sq_10ac_elevation.shp'))


# %% [markdown]
# ## load interpolated groundwater data


# %%
ghb_dir = join(gwfm_dir, 'GHB_data')
# elevations are in feet
ghb_fp = join(ghb_dir, 'interpolated_data')

# os.listdir(ghb_fp)
y=2014
raster_name = join(ghb_fp, 'fall'+str(y)+'_kriged.tif')
with rasterio.open(raster_name) as r:
    # print(r.meta)
    raster_crs = r.crs

# %%

# %%
parcel_gw = parcel_elev.copy()
parcel_gw = parcel_gw[grid_cols+['mean']].rename(columns={'mean':'gse_m'})
# gw_df.plot('mean', legend=True)
# sampling is pretty quick for an individual year
# sample the groundwater elevation to the parcel level 
for y in np.arange(2014,2021):
    raster_name = join(ghb_fp, 'fall'+str(y)+'_kriged.tif')
    gw_df = elev_stats(raster_name, parcel_id.to_crs(raster_crs))
    parcel_gw['gwe'+str(y)] =  gw_df['mean'] * 0.3048

# %%
# parcel_gw
# calculate the mean groudnwater elevation for the period to represent average conditions
gw_cols = 'gwe'+pd.Series(np.arange(2014,2021).astype(str)).values
parcel_gw['gwe_mean'] = parcel_gw[gw_cols].mean(axis=1)
# calculate depth to water
parcel_gw['dtw_mean'] = parcel_gw.gse_m - parcel_gw.gwe_mean

# %%
# parcel_gw.plot('dtw_mean', legend=True)

# %% [markdown]
# ## AEM data

# %%
aem_folder = 'statewide_aem_survey_coarse_fraction_depth_slices_and_average'
aem_domain_f = join(upw_dir, aem_folder, 'domain_aem_data.shp')

# %%
aem_depth = gpd.read_file(join(upw_dir, aem_folder, 'aem_depth.shp'))
aem_depth['thick'] = aem_depth.Depth_Bott - aem_depth.Depth_Top_ 


# %%
# slow to read because of high cell count
# aem_grid = gpd.read_file(join(upw_dir, aem_folder, 'aem_grid.shp'))


# %%
# # default of 10 bins+1 has 40 m thickness
# interval = pd.IntervalIndex(pd.cut(aem_grid.Elev_Mid_, bins=20))
# aem_grid['bin_lower'] = interval.left.values
# aem_grid['bin_higher'] = interval.right.values

# %% [markdown]
# # Visualize

# %%
vmin = aem_depth.PC_avg.min()
vmax = aem_depth.PC_avg.max()
depths =np.sort(aem_depth.Depth_Mid_.unique())[:6]
# depths =np.sort(aem_depth.Depth_Mid_.unique())[:-2:2]
# depths =np.sort(aem_grid.bin_higher.unique())[6:-2]
nx = 3
ny = int(np.ceil(len(depths)/nx))
fig,ax = plt.subplots(ny, nx, sharex=True, sharey=True, figsize=(6.5, 6.5*(ny/nx)), dpi=300)
for n, d in enumerate(depths):
    ax_n = ax[int(n/nx), n%nx]
    ax_n.annotate(str(d)+' m', (0.1,0.9), xycoords='axes fraction')
    # aem_plt = aem_grid[aem_grid.bin_higher==d]
    aem_plt = aem_depth[aem_depth.Depth_Mid_==d]
    aem_plt.plot('PC_avg', ax=ax_n, vmin = vmin, vmax=vmax)

fig.tight_layout(pad=-0.25)
plt.ticklabel_format(style='plain')

for n in np.arange(0,ny):
    ax[n,0].set_yticklabels(labels=ax[n,0].get_yticklabels(), rotation=90, verticalalignment = "center");
    ax[n,0].locator_params(axis='y', nbins=1);


# %%
# fig,ax= plt.subplots(1,2, sharey=True, figsize=(6.5,6.5), dpi=300)

# aem_depth.plot('K_m_d', ax=ax[0], legend=True, legend_kwds={'shrink':0.5})
# ax[0].set_title('Geometric mean K_m_d')
# aem_depth.plot('PC_avg', ax=ax[1], legend=True, legend_kwds={'shrink':0.5})
# ax[1].set_title('Geometric mean PC')
# cr.plot(ax=ax[0], color='black', linestyle='--', alpha=0.6)
# cr.plot(ax=ax[1], color='black', linestyle='--', alpha=0.6)

# %% [markdown]
# Estimate the vertical conductance for the unsaturated zone.
# - quick test just go down to 30 m
#
# - An important question to address with cbec is do they want seepage rates or recharge rates.
# - The averaging also needs to be scaled by the depth of each layer: multiply or exponent to account for additional impact.
#     - harmonic mean: $\frac{\Sigma b}{\frac{b_1}{x_1}+\frac{b_2}{x_2}}$
#     - geometric mean: $ \sqrt[n]{x_1^{b_1} x_2^{b_2}}$

# %%
# fill in hydraulic parameters
# params = pd.read_csv(join(upw_dir, 'ZonePropertiesInitial.csv'))
# params['K_m_d'] = params.K_m_s *86400
# maples parameters
params = pd.read_csv(join(upw_dir, 'ZonePropertiesInitial_Maples.csv'))

params.Lithology = params.Lithology.str.lower()
aem_depth['K_m_d'] = np.nan
for n, lith in enumerate(params.Lithology):
    aem_depth.loc[aem_depth.facies==lith, 'K_m_d'] = params.K_m_d.iloc[n]

# %%
from scipy.stats import gmean, hmean

# %%
# pull out the polygons for the AEM data
aem_surf = aem_depth[aem_depth.Depth_Top_==0].copy()
aem_surf = gpd.overlay(aem_surf, parcel_id.to_crs(avg_K.crs))

aem_depth_poly = aem_surf[['Id','Region','gridcode','geometry']]

# %%
# quick check, inaccurate
# avg_K = aem_depth.copy()
# avg_K = gpd.overlay(avg_K, parcel_id.to_crs(avg_K.crs))
# geom_K = avg_K[cols].groupby(grp_col).aggregate(gmean).reset_index()

# %%
# avg_K = aem_depth.copy()
# # crop on the parcel level
# avg_K = gpd.overlay(avg_K, parcel_id.to_crs(avg_K.crs))
# avg_K.plot()

# %%
avg_K = aem_depth.copy()
# crop on the parcel level
avg_K = gpd.overlay(avg_K, parcel_id.to_crs(avg_K.crs))
# preliminary cut for unsaturated zone (improve with kriged DTW)
avg_K = avg_K[avg_K.Depth_Bott<35]
cols = ['Id', 'Region', 'thick','Depth_Mid_', 'PC_avg','K_m_d']
grp_col = ['Id','Region']

# calculate geometric mean while accounting for thickness
val_cols = ['PC_avg', 'K_m_d']
avg_K[val_cols] = avg_K[val_cols].pow(avg_K.thick.values,axis=0)
geom_K = avg_K[grp_col+val_cols].groupby(grp_col).aggregate(np.prod).reset_index()
thick_sum = avg_K[cols].groupby(grp_col).sum()['thick'].values
geom_K[val_cols] = geom_K[val_cols].pow((1/thick_sum), axis=0)


# %%
def calc_gmean(K, thick):
    """ Given hydraulic conductivity and layer thickness
    calculate the geometric mean while accounting for thickenss """
    K_scale = K**thick
    thick_sum = np.sum(thick)
    geom_K = np.prod(K)**(1/thick_sum)
    return geom_K

# %%

# avg_K['K_scale'] = avg_K.K_m_d**avg_K.thick.values
# geom_K = avg_K[cols+['K_scale']].groupby(grp_col).aggregate(np.prod)
# thick_sum = avg_K[cols].groupby(grp_col).sum()['thick']
# geom_K['K_final'] = (geom_K.K_scale)**(1/thick_sum)
# geom_K = geom_K[['Depth_Mid_','K_final']].rename(columns={'K_final':'K_m_d'})
# the changes are with +- 30% between the original and updated scaling
# (geom_K.K_m_d/geom_K.K_final).plot()

# %%


# %%
# join to geodataframe for plotting
# geom_K_gdf = aem_depth_poly.merge(geom_K_scale.reset_index())
geom_K_gdf = aem_depth_poly.merge(geom_K)


# %%
geom_K_gdf.to_file(join(proj_dir, 'GIS', 'Kgeometric_mean.shp'))

# %%
fig,ax= plt.subplots(1,2, figsize=(6.5, 6.5), dpi=300, sharey=True)

# percent coarse has similar results as using vertical conductivity
# just the order of magnitude difference helps
geom_K_gdf.plot('PC_avg', ax=ax[0], legend=True, legend_kwds={'shrink':0.3},
               # norm = mpl.colors.LogNorm(vmin=1, vmax=100)
               )

cr.plot(ax=ax[0], color='black', linestyle='--', alpha=0.6, linewidth=0.5)
# geom_K_gdf.plot('PC_avg', ax=ax[1], legend=True, legend_kwds={'shrink':0.5})
geom_K_gdf.plot('K_m_d', ax=ax[1], legend=True, legend_kwds={'shrink':0.3},
                norm = mpl.colors.LogNorm(vmin = geom_K_gdf.K_m_d.min(), vmax = geom_K_gdf.K_m_d.max())
)
cr.plot(ax=ax[1], color='black', linestyle='--', alpha=0.6, linewidth=0.5)

ax[0].set_title('Percent Coarse')
ax[1].set_title('$K_{geom}$ (m/day)')
# plt_cln(ax=ax[0])

# %% [markdown]
# The shallow upscaling does help show some hot spots near rooney and another mid-way up the Cosumnes.

# %% [markdown]
# ## Recharge estimate  
# Use regression from Steve's work to estimate the potential recharge for the sites  
# $ y = 0.0376 x +5.288 $  returns the 30-day average recharge rate in cm/day
# - initial estimate will use the geometric hydraulic conductivity to 30 m. Final results should subset by unsaturated zone
# - also the current work is going to use the 2014-2020 average dtw since this period includes a mix of dry and wet years, but the final could use the most recent contours available since this is the starting condition

# %%
# preliminary test of the recharge potential based on Steve's work
rech_est = geom_K_gdf.merge(parcel_gw)
rech_est = parcel_gw.merge(geom_K)
rech_est['geomK_dtw'] = rech_est.K_m_d*rech_est.dtw_mean

# 30 day average recharge estimate (cm/day)
rech_est['rch_cm_d'] = rech_est.geomK_dtw*0.0376+5.288

# %%
rech_est.to_file(join(gis_dir,'recharge_results','recharge_estimate.shp'))

# %%
# log scale doesn't improve things much
rech_est.plot('rch_cm_d',legend=True,
             # norm = mpl.colors.LogNorm(vmin = rech_est.rch_cm_d.min(), vmax = rech_est.rch_cm_d.max())
              legend_kwds = {'shrink':0.7,'label':'30-day average recharge (cm/day)'}
             )

# %% [markdown]
# The 30 day recharge estimate looks like it will be helpful in characterizing additionally simply by making the cutoff with depth to water and hydraulic conductivity simpler to look at in one image with the knowledge that an unsat model predicted similar recharge rates.
# - it makes sense the middle reaches are highest because they have the deepest DTW and higher Kgeom values
# - not sure about the gap between rooney and teichert is about

# %% [markdown]
# # WCR
# If possible it might be valuable to reduce the number of WCRs or focus on those considered more trustworthy

# %%
# load WCR summary data
wcr_list = pd.read_csv(join(gwfm_dir, 'WEL_data','wellcompletionreports.csv'))

# %%
# well completion reports
wcr = pd.read_csv(join(proj_dir, 'WCR_all_cleaned.csv'), index_col=0)


# %%

# %%
wcr_path = join(dirname(gwfm_dir), 'Large_TPROGS_run')
wcr_loc_file = join(wcr_path, 'All_locations_Large_Cosumnes_run.csv') # has all locations
wcr_file = join(wcr_path, 'TPROGS_update_full_run_12_2020_minus_NaN_and_duplicates.csv') # has all locations and depths

id_col = 'WCR Number'
wcr_df_in = pd.read_csv(wcr_file)
wcr_gdf = pd.read_csv(wcr_loc_file)
wcr_gdf_full = gpd.GeoDataFrame(wcr_gdf, 
                                geometry=gpd.points_from_xy(wcr_gdf.Easting, wcr_gdf.Northing), crs='epsg:3310')

# %%
from_col = wcr_df_in['Elevation (m)'] - wcr_df_in['Relative Elevation (m)']
to_col = from_col +0.3048
from_to = pd.concat((from_col, to_col),axis=1)[wcr_df_in['from'].isna()].values
wcr_df_in.loc[wcr_df_in['from'].isna(), ['from','to']] = from_to
# wcr_df_in.loc[wcr_df_in['from'].isna(), 'to'] = pd.concat((from_col, to_col),axis=1)[wcr_df_in['from'].isna()]
wcr_pre_adj =  wcr_df_in.copy()
# from_col

# %% [markdown]
# Try filtering by date, whether perforated interval/total depth exists, or planned use
# - if total completed depth is below -80 m then there might have been less care about shallower aquifer material
# - planned use: urban wells likely most accurate
# - newer wells likely more precise than older
#     - filtering by age didn't change the visual perception significantly, still high density in spots
# - depth
#     - \>1000 wells were \<100 m, only 78 wells were \<50 m deep, 560 were \<75 m deep
#     - 75 m depth was probably best cutoff in removing excess wells in wilton area from domestic use
# - filter by age and depth didn't make things clearer either
#
# **doesn't seem like WCRs by themselves can help identify hot spots for recharge**

# %%

wcr_list_domain = wcr_list[wcr_list.WCRNumber.isin(wcr_gdf_full['WCR Number'])].copy()
wcr_list_domain.columns
print('Num wells', len(wcr_list_domain))
# planned use
# wcr_list_domain.TopOfPerforatedInterval
# wcr_list_domain.TotalCompletedDepth
# wcr_list_domain.PlannedUseFormerUse
wcr_list_domain.DateWorkEnded = pd.to_datetime(wcr_list_domain.DateWorkEnded, errors='coerce')
min_date = '1980-10-1'
wcr_list_domain = wcr_list_domain[wcr_list_domain.DateWorkEnded>min_date]
print('Num wells', len(wcr_list_domain), 'post',min_date)
max_depth=100/0.3048 # focus on wells that are likely accurate in the domain
# 50 is way to short of a cutoff with only 7
# max_depth = 75/0.3048 # focus on wells in shallow region
wcr_list_domain = wcr_list_domain[wcr_list_domain.TotalCompletedDepth< max_depth]
print('Num wells', len(wcr_list_domain), 'shallower than %.2f ft'%max_depth)

wcr = wcr_pre_adj[wcr_pre_adj['WCR Number'].isin(wcr_list_domain.WCRNumber)].copy()

# %%
# wcr = wcr_pre_adj.copy()

# insert hydraulic conductivity value to well completion reports
wcr['K_m_d'] = np.nan
for n, z in enumerate(params.Zone):
    wcr.loc[wcr.classification==z, 'K_m_d'] = params.K_m_d.iloc[n]

# %%
# convert to geodataframe
# wcr_gdf = gpd.GeoDataFrame(wcr, geometry = gpd.points_from_xy(wcr.Easting, wcr.Northing), crs='epsg:3310')
# wcr_unique = wcr.drop_duplicates('WCR Number')
# wcr_gdf = gpd.GeoDataFrame(wcr_unique, 
#                            geometry = gpd.points_from_xy(wcr_unique.Easting, wcr_unique.Northing),
#                            crs = 'epsg:3310')
wcr_gdf = wcr_gdf_full.copy()
wcr_gdf = wcr_gdf[['WCR Number', 'geometry']]

# %%
# sample shallow wells
wcr_avg = wcr[wcr.to < 90].copy() # this accidentally removes wells
wcr_cols= ['WCR Number','from','to','classification','K_m_d']
wcr_arith = wcr_avg[wcr_cols].groupby('WCR Number').mean()
wcr_hmean = wcr_avg[wcr_cols].groupby('WCR Number').aggregate(hmean)
wcr_gmean = wcr_avg[wcr_cols].groupby('WCR Number').aggregate(gmean)
wcr_gmean['kind'] = 'geometric'
wcr_hmean['kind'] = 'harmonic'
wcr_arith['kind'] = 'arithmetic'

# %%
wcr_mean_all = pd.concat((wcr_gmean, wcr_hmean, wcr_arith))
wcr_mean_gdf = wcr_gdf.merge(wcr_mean_all.reset_index())


# %%
vmin, vmax = wcr_mean_all.K_m_d.quantile([0,1])

# %%
k= 'geometric'
# k='arithmetic'
fig, ax = plt.subplots(figsize=(3,3), dpi=300)
wcr_mean_gdf[wcr_mean_gdf.kind==k].plot('K_m_d', ax=ax, 
                                        legend=True, markersize=0.5,
                                       norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
                                       )
ctx.add_basemap(source=ctx.providers.Esri.WorldImagery, ax=ax, alpha = 0.6,
                crs='epsg:3310', attribution=False)



# %% [markdown]
# Laura noted that we should give them all of the possible datasets and ask that they do a sensitivity analysis in their methodology. Or we can do the sensitivity analysis by comparing the 30 day recharge estimated.
