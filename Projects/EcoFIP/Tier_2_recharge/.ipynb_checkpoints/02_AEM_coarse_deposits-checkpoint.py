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
# gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'


# %%
lwa_dir = join(usr_dir, 'LWA Dropbox','01_Project-Teams')
proj_dir = join(lwa_dir, '669.03 - DWR Cosumnes Floodplain Recharge')
# gis_dir = join(main_concept_dir,'GIS')

# should start using files for EcoFIP stored on Dropbox instead of Box
gwfm_dir = join(proj_dir, 'data','GWFlowModel')
gwfm_proj_dir = join(gwfm_dir, 'Projects','EcoFIP')

# %%
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
add_path(join(doc_dir, 'GitHub','flopy'))

import map_cln
reload(map_cln)

from mf_utility import get_layer_from_elev
from report_cln import base_round
from map_cln import gdf_bnds, lab_pnt, plt_cln, arr_lab, xy_lab, make_multi_scale, dir_arrow, plt_arrow


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
# grid_id= 'parcel'
# grid_id = 'sq'
grid_id = 'mf'

if grid_id=='parcel':
    gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'parcels_elevation.shp'))
elif grid_id == 'sq':
    gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'sq_10ac_elevation.shp'))
elif grid_id == 'mf':
    gdf_elev = gpd.read_file(join(gwfm_dir, 'DIS_data', 'grid','grid.shp'))
    # rename node to Region
    gdf_elev = gdf_elev.rename(columns={'node':'Region'})

# simpler geodataframe to bring dataframe to geodataframe
gdf_id = gdf_elev[['Region','geometry']].copy()


# %%
fig_dir = join(proj_dir, 'figures', 'AEM_coarse', grid_id)
os.makedirs(fig_dir, exist_ok=True)

# %%
from map_cln import plt_arrow, make_multi_scale
# function for within this script to clean up a plot
def clean_plot(ax_n):
    ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

    gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
    # gdf_bnds(gdf_id, ax=ax_n, buf=2E3)
    
    ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
                    crs='epsg:26910', attribution=False)
    plt_arrow(ax_n, 0.925, 0.2)
    make_multi_scale(ax_n, 0.65, 0.1, dist = 2E3, scales = [4,2,1], units='km')

# %% [markdown]
# ## load interpolated groundwater data
# Spring data is better to use as it will be more conservative in terms of the available space for recharge


# %%
season='spring'


# %%
if grid_id=='mf':
        gdf_gw_long.to_csv(join(gis_dir, 'analysis_unit_reference', 
                         'GW_elevations_long_'+grid_id+'_'+season+'.csv'))
else:
    gdf_gw_long = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 
                             'GW_elevations_long_'+grid_id+'_'+season+'.shp'))
gdf_gw_long.year = gdf_gw_long.year.astype(int)

# %%

# %% [markdown]
# For fall groundwater elevations that parcels have a mean standard deviation of 2.64 m across all years (2000-2020), for spring that goes down to 2.28 m.
#
# Looking at sorted values from different time periods it seems like we are getting more extremes in the basin.

# %%
# compare groundwater elevation distribution by period
gw_early_mean = gdf_gw_long[gdf_gw_long.year<=2010].groupby('Region').mean(numeric_only=True)
gw_late_mean = gdf_gw_long[gdf_gw_long.year>2010].groupby('Region').mean(numeric_only=True)
gw_recent_mean = gdf_gw_long[gdf_gw_long.year>=2018].groupby('Region').mean(numeric_only=True)
plt.plot(gw_early_mean.dtw_m.sort_values().values, label='2000-2010')
plt.plot(gw_late_mean.dtw_m.sort_values().values,label='2011-2020')
plt.plot(gw_recent_mean.dtw_m.sort_values().values, label='2018-2020')
plt.legend()

plt.close()

# %%
# sample the recent period for the analysis
gdf_gw = gdf_gw_long[gdf_gw_long.year>2010].dropna(subset='dtw_m').groupby('Region').mean(numeric_only=True)
gdf_gw = gdf_id.merge(gdf_gw.reset_index())

fig,ax_n = plt.subplots(figsize=(6.5,5.5),dpi=300)
gdf_gw.plot('dtw_m', ax=ax_n, legend=True, 
            legend_kwds = {'label':'Depth to water (m)','shrink':0.6})

# plt_arrow(ax_n, 0.925, 0.20)
# make_multi_scale(ax_n, 0.6, 0.1, dist = 2E3, scales = [4,2,1], units='km')
# ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

# gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
# # gdf_bnds(gdf_id, ax=ax_n, buf=2E3)

# ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
#             crs='epsg:26910', attribution=False)

clean_plot(ax_n)
plt.savefig(join(fig_dir, 'depth_to_water_m.png'), bbox_inches='tight')


# %% [markdown]
# ## AEM data

# %%
aem_folder = 'statewide_aem_survey_coarse_fraction_depth_slices_and_average'
# aem_domain_f = join(upw_dir, aem_folder, 'domain_aem_data.shp')

# %%
aem_depth = gpd.read_file(join(upw_dir, aem_folder, 'aem_depth.shp'))
aem_depth['thick'] = aem_depth.Depth_Bott - aem_depth.Depth_Top_ 


# %%
# pull out the polygons for the AEM data
aem_surf = aem_depth[aem_depth.Depth_Top_==0].copy()
aem_surf = gpd.overlay(aem_surf, gdf_id.to_crs(aem_surf.crs))

aem_depth_poly = aem_surf[['Id','Region','gridcode','geometry']]

# %%
# use the 0 top depth poly is the same as at 5
# aem_surf = aem_depth[aem_depth.Depth_Top_ ==5]
# aem_surf = gpd.overlay(aem_surf, gdf_id.to_crs(aem_surf.crs))
# aem_depth_poly2 = aem_surf[['Id','Region','gridcode','geometry']]
# fig,ax = plt.subplots()
# aem_depth_poly2.plot(ax=ax,color='red')
# aem_surf.plot(ax=ax)


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
# # convert aem to wide format for quick plotting in shapefile
# aem_shp_out = aem_depth[aem_depth.Depth_Mid_<30]
# aem_shp_out = aem_shp_out.pivot(columns='Depth_Mid_', index='Id', values='PC_avg')
# aem_shp_out = aem_surf[['Id','geometry']].merge(aem_shp_out.reset_index())
# aem_shp_out.columns = aem_shp_out.columns.astype(str)
# aem_shp_out.to_file(join(gis_dir,'shallow_percent_coarse_by_depth_cols.shp'))

# %%
np.sort(aem_depth.Depth_Bott.unique())

# %%
vmin = aem_depth.PC_avg.min()
vmax = aem_depth.PC_avg.max()
depths =np.sort(aem_depth.Depth_Mid_.unique())[:6]
# depths =np.sort(aem_depth.Depth_Mid_.unique())[:-2:2]
# depths =np.sort(aem_grid.bin_higher.unique())[6:-2]
nx = 3
ny = int(np.ceil(len(depths)/nx))
fig,ax = plt.subplots(ny, nx, sharex=True, sharey=True, figsize=(6.5, 5.5*(ny/nx)), dpi=300)

for n, d in enumerate(depths):
    ax_n = ax[int(n/nx), n%nx]
    ax_n.annotate(str(d)+' m', (0.1,0.9), xycoords='axes fraction')
    # aem_plt = aem_grid[aem_grid.bin_higher==d]
    aem_plt = aem_depth[aem_depth.Depth_Mid_==d]
    aem_plt.plot('PC_avg', ax=ax_n, vmin = vmin, vmax=vmax)
    # Esri.WorldStreetMap
    ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
                crs='epsg:26910', attribution=False)
    ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

    gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
    gdf_bnds(gdf_id, ax=ax_n, buf=2E3)

# fig.tight_layout(pad=-0.25) # for square plot
fig.tight_layout(h_pad=0.25, w_pad=0.25) # rectangular
plt.ticklabel_format(style='plain')

# for n in np.arange(0,ny):
    # ax[n,0].set_yticklabels(labels=ax[n,0].get_yticklabels(), rotation=90, verticalalignment = "center");
    # ax[n,0].locator_params(axis='y', nbins=1);
ax_n = ax[0,0]
# arrow and scale sizing isn't workign
# plt_arrow(ax_n, 0.925, 0.15)
# make_multi_scale(ax_n, 0.6, 0.1, dist = 2E3, scales = [4,2,1], units='km')



# %%
def plt_arrow(ax, xoff = 0.7, yoff=0.15):
    x, y, arrow_length = xoff, yoff, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20, 
                xycoords=ax.transAxes)
    return None

# %%

# import h5py
# with h5py.File(join(upw_dir, aem_folder, 'aem_array.hdf5'), 'r') as f:
#     arr = f['facies']['facies_array'][:]
#     arr_cf = f['percent_coarse']['pc_array'][:]

# %%
# an option represent vertical connectivity

# sfr_profile = arr[:,grid_sfr.row.astype(int)-1, grid_sfr.column.astype(int)-1]
# sfr_profile_cf = arr_cf[:,grid_sfr.row.astype(int)-1, grid_sfr.column.astype(int)-1]
# fig,ax = plt.subplots(2,1,sharex=True, figsize=(6.5,6))
# ax[0].imshow(sfr_profile, cmap='viridis_r')
# ax[1].imshow(sfr_profile_cf, cmap='viridis')
# ax[0].set_aspect(1/5)
# ax[1].set_aspect(1/5)
# fig.tight_layout(h_pad=-4)

# %%
# # # quick way to make a figure legend, to crop in for powerpoint
# fig,ax_n = plt.subplots()
# im = aem_plt.plot('PC_avg', ax=ax_n, legend=True, legend_kwds={'label':'Percent Coarse'},
#                   vmin = vmin, vmax=vmax)

# plt_arrow(ax_n, 0.925, 0.15)
# make_multi_scale(ax_n, 0.7, 0.1, dist = 2E3, scales = [4,2,1], units='km')


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


# %%
# [0, 40, 50, 80,100] > and <= percent for each facies
# mud is 0-10, sandy mud 40-50, sand is 70-80, gravel is 90-100
params_pc = params.set_index('Lithology').loc[['mud','sandy mud','sand','gravel']].reset_index().copy()
# assign each facies to a percent coarse average
params_pc['PC_avg'] = [5, 45, 65, 95]
params_pc = params_pc.set_index('PC_avg').reindex(np.sort(aem_depth.PC_avg.unique()))
# create a log K to have a log step between K parameters
params_pc['log_K_m_d'] = np.log10(params_pc.K_m_d)
# linearly interpolate K for each percent coarse
params_pc = params_pc.interpolate('linear').reset_index()
# calculate the log interpolated K value in normal space
params_pc['K_m_d_intrp'] = 10**(params_pc.log_K_m_d)
# compare the standard and log interpolated K
params_pc.plot(x='PC_avg', y=['K_m_d','K_m_d_intrp'])
plt.yscale('log')
plt.close()
# the biggest impact is in the low K facies Mud to sandy mud and sandy mud to sand

# %%
aem_depth['K_m_d_old'] = np.nan
aem_depth['K_m_d'] = np.nan
# original method uses exact K value for a subset of percent coarse
for n, lith in enumerate(params.Lithology):
    aem_depth.loc[aem_depth.facies==lith, 'K_m_d_old'] = params.K_m_d.iloc[n]

# updated method uses interpolated K between percent coarse
for n, pc_avg in enumerate(params_pc.PC_avg):
    aem_depth.loc[aem_depth.PC_avg==pc_avg, 'K_m_d'] = params_pc.K_m_d_intrp.iloc[n]
    # need to include specific yield for pore space filling
    aem_depth.loc[aem_depth.PC_avg==pc_avg, 'Sy'] = params_pc.Sy.iloc[n]

# %%
frac_coarse = aem_depth[aem_depth.K_m_d>10].area.sum()/aem_depth.area.sum()
print('Fraction coarse %.2f' %frac_coarse)
frac_coarse = aem_depth[aem_depth.K_m_d_old>10].area.sum()/aem_depth.area.sum()
print('Fraction coarse %.2f' %frac_coarse)


# %% [markdown]
# In the updated method with interpolated K by percent coarse, I found that it had reduce the area >10 m/day to only 4%. This with the sand set as 75% coarse, lowering to 65% increased this 10%. Lowering sandy mud to 35% didn't make a difference for coarse (log scale interpolation kept 55 too low). Setting sand at 55% brought percent coarse to 19%.
# - the issue is that with sand at 75% coarse then 0% of the domain ends up with recharge >25 cm/day but at 55% it ends up with 20% >25 cm/day and only 50% < 10 cm/day when it should be nearly 84%
# - using sand at 65% seems most reasonable since only 2% are >25 cm/day and %74 are <10 cm/day which is about what it should be for the cosumnes

# %%
from scipy.stats import gmean, hmean


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
def calc_unsat_thickness(avg_K):
    # any aem with their bottom above depth to water is fully useable
    avg_K_unsat = avg_K[avg_K.Depth_Bott<avg_K.dtw_m].copy()
    # data with water depth in the layer should have thickness subset to unsat thickness
    avg_K_part_unsat = avg_K[(avg_K.Depth_Top_ < avg_K.dtw_m)&(avg_K.Depth_Bott>avg_K.dtw_m)].copy()
    # update the thickness based on the unsaturated thickenss
    avg_K_part_unsat['thick'] = avg_K_part_unsat.dtw_m - avg_K_part_unsat.Depth_Top_
    # the line above could be done for all cells and then replace thickness smaler than the original
    # and drop negative thicknesses
    avg_K = pd.concat((avg_K_unsat, avg_K_part_unsat))
    return avg_K


# %%
# 

# %%
avg_K = aem_depth.copy()
# crop on the parcel level
avg_K = gpd.overlay(avg_K, gdf_id.to_crs(avg_K.crs))

# # join K data with dtw to subset
avg_K = avg_K.merge(gdf_gw[['Region','dtw_m']], how='inner')
# # preliminary cut for unsaturated zone (improve with kriged DTW)
avg_K = avg_K[avg_K.Depth_Bott<35]
# the unsat_thickness calc removes areas if they are all saturated
avg_K = calc_unsat_thickness(avg_K)
cols = ['Id', 'Region', 'thick','Depth_Mid_', 'PC_avg','K_m_d']
# Id is for AEM, Region is for EcoFIP
grp_col = ['Id','Region']

# calculate geometric mean while accounting for thickness
val_cols = ['PC_avg', 'K_m_d', 'Sy']
avg_K[val_cols] = avg_K[val_cols].pow(avg_K.thick.values,axis=0)
geom_K = avg_K[grp_col+val_cols].groupby(grp_col).aggregate(np.prod).reset_index()
thick_sum = avg_K[cols].groupby(grp_col).sum()['thick'].values
geom_K[val_cols] = geom_K[val_cols].pow((1/thick_sum), axis=0)

# %%
# join to geodataframe for plotting
# geom_K_gdf = aem_depth_poly.merge(geom_K)


# %%
# the maximum K_m_d for each region should be selected
# need to drop Id
geom_K_avg = geom_K.groupby('Region').mean().reset_index().drop(columns=['Id'])
#geom_K_avg.plot(y='K_m_d')
geom_K_gdf = aem_depth_poly.merge(geom_K_avg)
# geom_K_gdf.plot('K_m_d')


# %%
# geom_K_gdf.plot('Sy',  legend=True, legend_kwds={'shrink':0.3},
#                # norm = mpl.colors.LogNorm(vmin=1, vmax=100)
#                )
# the value range of Sy is small only 0.25 to 0.3
# params

# %%
fig,ax= plt.subplots(1,2, figsize=(6.5, 6.5), dpi=300, sharey=True)

# percent coarse has similar results as using vertical conductivity
# just the order of magnitude difference helps
geom_K_gdf.plot('PC_avg', ax=ax[0], legend=True, legend_kwds={'shrink':0.3},
               # norm = mpl.colors.LogNorm(vmin=1, vmax=100)
               )

# geom_K_gdf.plot('PC_avg', ax=ax[1], legend=True, legend_kwds={'shrink':0.5})
geom_K_gdf.plot('K_m_d', ax=ax[1], legend=True, legend_kwds={'shrink':0.3},
                norm = mpl.colors.LogNorm(vmin = geom_K_gdf.K_m_d.min(), vmax = geom_K_gdf.K_m_d.max() )
)

ax[0].set_title('Percent Coarse')
ax[1].set_title('$K_{geom}$ (m/day)')
# plt_cln(ax=ax[0])

ax_n = ax[0]
# plt_arrow(ax_n, 0.925, 0.3)
# make_multi_scale(ax_n, 0.6, 0.1, dist = 2.5E3, scales = [4,2,1], units='km')
for n in [0,1]:
    ax_n = ax[n]
    cr.plot(ax=ax_n, color='black', linestyle='--', alpha=0.6, linewidth=0.5)

    ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)
    
    gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
    # gdf_bnds(gdf_id, ax=ax_n, buf=2E3)
    
    ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
                crs='epsg:26910', attribution=False)
    
fig.tight_layout()
plt.savefig(join(fig_dir, 'percent_coarse_K_geometric_mean.png'), bbox_inches='tight')

# %%
geom_K_gdf.to_file(join(proj_dir, 'GIS', 'Kgeometric_mean_'+grid_id+'.shp'))

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
# rech_est = geom_K_gdf.merge(gdf_gw)
# join by grid Region only
# rech_est = gdf_gw.merge(geom_K) # old version had multiple polygons per grid
rech_est = gdf_gw.merge(geom_K_avg)
# aggregate to the grid region level
rech_est['geomK_dtw'] = rech_est.K_m_d*rech_est.dtw_m

# 30 day average recharge estimate (cm/day)
rech_est['rch_cm_d'] = rech_est.geomK_dtw*0.0376+5.288

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
             norm = mpl.colors.LogNorm(vmin = rech_est.rch_cm_d.min(), vmax = rech_est.rch_cm_d.max()),
              legend_kwds = {'shrink':0.7}
             )
plt.title('30-day average recharge (cm/day)')
# ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

# gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
# # gdf_bnds(gdf_id, ax=ax_n, buf=2E3)

# ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
#                 crs='epsg:26910', attribution=False)
clean_plot(ax_n)
# plt.savefig(join(fig_dir, 'recharge_30d_avg.png'), bbox_inches='tight')

# %%
rech_est.to_file(join(gis_dir,'recharge_results','recharge_estimate_'+grid_id+'.shp'))

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
