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

# %%
# standard python utilities
import os
from os.path import basename, dirname, join, exists
import sys
from importlib import reload
import glob
import pandas as pd
import numpy as np
import calendar
import time

# standard python plotting utilities
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# standard geospatial python utilities
# import pyproj # for converting proj4string
# import shapely
import geopandas as gpd
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator

# import flopy
# import flopy.utils.binaryfile as bf
from importlib import reload
import h5py



# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

bas_dir = join(gwfm_dir, 'BAS6')
proj_dir = join(gwfm_dir, 'Projects','EcoFIP')
plt_dir = join(proj_dir,'figures/')
out_dir = join(proj_dir,'output')
gis_dir = join(proj_dir, 'GIS')


# %%
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

from map_cln import gdf_bnds, plt_cln, xy_lab
from report_cln import base_round

from mf_utility import get_dates, clean_hob


# %%
from map_cln import plt_arrow, make_multi_scale
# function for within this script to clean up a plot
def clean_plot(gdf_id, ax_n):
    gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)

    ax_n.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

    ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
                    crs='epsg:26910', attribution=False)
    plt_arrow(ax_n, 0.925, 0.2)
    make_multi_scale(ax_n, 0.65, 0.1, dist = 2E3, scales = [4,2,1], units='km')


# %%
# method to set default plot parameters, no longer need to specify 300 dpi each time, may not need to specify dimensions either
plt.rcParams.update({"figure.dpi": 600})

# %% [markdown]
# # Load data
#
# 1. The quick method to estimate recharge rates for the areas is to do a nearest spatial join with the sfr grid cells 
# 2. While we only have stream seepage data for the existing stream channel, we could do a regression like Steve did to extrapolate to nearby cells what stream seepage rates we might expect based on the streambed conductance and/or stream depth. 
#

# %% [markdown]
# ## spatial reference

# %%
gdf_sfr = gpd.read_file(join(gis_dir, 'sfr_reach_reference.shp'))

# %%
line_ref = gpd.read_file(join(gis_dir, 'key_road_lines_reference.shp'))
# line_ref.plot()

# %%
grid_id = 'parcel'
grid_id = 'sq'

if grid_id=='parcel':
    gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'parcels_elevation.shp'))
elif grid_id == 'sq':
    gdf_elev = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 'sq_10ac_elevation.shp'))

# simpler geodataframe to bring dataframe to geodataframe
gdf_id = gdf_elev[['Region','geometry']].copy()
gdf_id['grid_area'] = gdf_id.area

# %%
season='spring'
gdf_gw_long = gpd.read_file(join(gis_dir, 'analysis_unit_reference', 
                         'GW_elevations_long_'+grid_id+'_'+season+'.shp'))
gdf_gw_long.year = gdf_gw_long.year.astype(int)

# %%
fig_dir = join(proj_dir, 'figures', 'recharge', grid_id)
os.makedirs(fig_dir, exist_ok=True)

# %% [markdown]
# # comparison
# Compare the hydrologic soil group maps, the recharge potential estimates and simulated seepage from the riverbed.

# %% [markdown]
# ## Hydrologic soil group

# %%
hsg = gpd.read_file(join(gis_dir, 'Combined_County_Soil_HydGrp_ClippedCos','Combined_County_Soil_HydGrp_Clipped.shp'))
# csv with Ksat in feet/day for each hydrologic soil group
hsg_K = pd.read_csv(join(gis_dir, 'HSG_K_ft_d.csv'))
# add K to geodataframe with merge by hydrologic soil group
hsg = hsg.merge(hsg_K.rename(columns={'HSG':'Hyd_group'}))
# transfer to project crs
hsg = hsg.to_crs(gdf_id.crs)

# %%
# overlay analysis unit onto hsg
hsg_grid = gpd.overlay(hsg, gdf_id)

# need to get average or representative rate for each cell (area based mean or majority)
hsg_grid['frac_area'] = hsg_grid.area/hsg_grid.grid_area
hsg_grid['frac_K'] = hsg_grid.K_ft_d*hsg_grid.frac_area
hsg_grid_agg = hsg_grid.groupby('Region').agg({'Hyd_group':pd.Series.mode,'frac_K':'sum','K_ft_d':'mean', 'frac_area':'sum'})
# the aggregation by area or simple mean doesn't appear to impact spatial distribution
hsg_grid_agg = gdf_id.merge(hsg_grid_agg.reset_index())

# there are quite a few cells whose aggregated fractional area is well under 1
# this means that the average will be too low because it's missing area
print('%.i with fractional area less than 0.9' %hsg_grid_agg[hsg_grid_agg.frac_area<0.9].shape[0])
# correct for cells with fractional area less than 1 by dividing by sum of fraction area
hsg_grid_agg['K_ft_d_adj'] = hsg_grid_agg.frac_K/hsg_grid_agg.frac_area

# %%
# aggregating to the grid analysis does have an impact beause parcels with lower K are reduced from a high value
fig,ax = plt.subplots(2,1)
hsg_grid_agg.plot('K_ft_d_adj', ax=ax[0])
hsg_grid.plot('K_ft_d', ax=ax[1])
plt.close()


# %%
# create new hsg_grid for output
hsg_grid = hsg_grid_agg[['Region','Hyd_group','K_ft_d_adj','geometry']].rename(columns = {'K_ft_d_adj':'K_ft_d'})

# %%
fig,ax_n = plt.subplots()
# hsg.plot('K_ft_d', ax=ax_n, legend=True)
hsg_grid.plot('K_ft_d',ax=ax_n,legend=True,
              legend_kwds={'shrink':0.7}
             )
# gdf_id.plot(ax=ax_n, color='None', linewidth=0.1, alpha=0.7)
# # gdf_bnds(gdf_id, ax=ax_n, buf=2E3)

# ctx.add_basemap(source= ctx.providers.OpenStreetMap.Mapnik, ax=ax_n, alpha = 0.6,
#                 crs='epsg:26910', attribution=False)
clean_plot(gdf_id, ax_n)
plt.title('Hydrologic Soil Group Ksat (ft/day)')
plt.savefig(join(fig_dir, 'hyd_soil_group.png'), bbox_inches='tight')

# %% [markdown]
# - The hydrologic soil group has fairly large units and on the lower Cosumnes the floodplan between Deer Creek and the Cosumnes is sandy which results in a huge chunk of the floodplain predicting with rates of 4.5 ft/day. A little excessive.  
# - When averaging by parcel the low K soils end up reducing the higher K sands to average rate

# %%
# soil data cropped to the model domain
# units of Ksat are m/day
# the data here are from STATSGO (national survey in 2016)
# it's not clear exactly how statsgo calculates
soils = gpd.read_file(join(gwfm_dir, 'UZF_data', 'clean_soil_data', 'cleaned_spatial_soil.shp'))
soils['K_ft_d'] = soils.Ksat /0.3048
soils_grid = gpd.overlay(soils, gdf_id)

# %%
# soils.plot('HydGroup', legend=True)
var = 'K_ft_d'
soils_grid.plot(var, legend=True,
           # norm = mpl.colors.LogNorm(soils[var].min(), soils[var].max())
          )
plt.close()
# the statsgo data shows the same result as the hydrologic soil group

# %% [markdown]
# ## Modeled seepage

# %%
strm_seep = False
if strm_seep:
    q=50
    sfr_rch = gpd.read_file(join(out_dir, 'stream_loss_'+grid_id, 'stream_loss_p'+str(q)+'.shp'))
    # calculate rate in ft per day for consistency
    sfr_rch['Qaq_ft_d'] = sfr_rch.Qaq_rate/0.3048


# %%
if strm_seep:
    var = 'Qaq_ft_d'
    vmin = 1E-2
    vmax = 1E1
    sfr_rch_plt = sfr_rch.copy()
    sfr_rch_plt.loc[sfr_rch_plt[var]<0, var] = np.nan
    # vmin = sfr_rch_plt[var].min()
    # vmax = sfr_rch_plt[var].max()
    fig,ax_n = plt.subplots()
    sfr_rch_plt.plot(var, ax=ax_n, legend=True,
               norm = mpl.colors.LogNorm(vmin, vmax),
                     legend_kwds={'shrink':0.7, 'extend':'both'}, #
              )
    # there is one site that is over-estimating reality with rates of 17.5 ft/day
    clean_plot(gdf_id, ax_n)
    
    plt.title('Average of nearby stream seepage (ft/day)')
    plt.savefig(join(fig_dir, 'avg_stream_seepage.png'), bbox_inches='tight')


# %% [markdown]
# ## Minimum of soil map and stream seepage
# Identify where the soil map data is limiting

# %%
def get_min_rate(rate_comp, rate1_name, rate2_name):
    """ 
    Take a dataframe with two rates to create a new column with the minimum of the two rates
    Input:
    rate_comp: dataframe with columns rate1_name and rate2_name
    Output:
    rate_comp: dataframe with additional columns
        lim: identifies when rate2 was the minimum
        lim_rate: minimum of both rates supplied
    """
    # identify where the aquifer rates will be limiting
    rate_comp['lim'] = False
    lim_rate = rate_comp[rate2_name] < rate_comp[rate1_name]
    rate_comp.loc[lim_rate, 'lim']= True
    # create a composite rate based on the supplied rates
    rate_comp['lim_rate'] = rate_comp[rate1_name]
    rate_comp.loc[lim_rate, 'lim_rate'] = rate_comp.loc[lim_rate, rate2_name]
    return(rate_comp)


# %%
if strm_seep:
    # testing the using of left join shows the soil map data would result in bad overweighting
    # of soil data only zones
    hsg_sfr_comp = hsg_grid.merge(sfr_rch_plt[['Region','Qaq_ft_d']], how='inner')
    hsg_sfr_comp = get_min_rate(hsg_sfr_comp, 'K_ft_d', 'Qaq_ft_d').rename(columns={'lim_rate':'lim_K_ft_d'})

# %%
# # calculate fraction of area limited by aquifer vs soil 
# hsg_sfr_comp['area'] = hsg_sfr_comp.geometry.area
# hsg_sfr_comp.groupby('lim')['area'].sum()/hsg_sfr_comp.area.sum()

# %%
gdf_sfr_buf = gdf_sfr.copy()
gdf_sfr_buf.geometry = gdf_sfr_buf.buffer(500)

# %%
if strm_seep:
    # hsg_sfr_comp.plot('Qaq_lim', legend=True)
    fig,ax_n = plt.subplots(figsize=(6.5,5.5))
    var = 'lim_K_ft_d'
    vmin = 1E-2
    vmax = 1E1
    # sfr_rch_plt = sfr_rch.copy()
    # sfr_rch_plt.loc[sfr_rch_plt[var]<0, var] = np.nan
    # vmin = hsg_sfr_comp[var].min()
    vmax = hsg_sfr_comp[var].max()
    hsg_sfr_comp.plot(var, ax=ax_n, legend=True,
                     legend_kwds={'shrink':0.7},
                     # norm = mpl.colors.LogNorm(vmin, vmax),
                     )
    clean_plot(gdf_id, ax_n)
    # gdf_sfr_buf.plot(ax=ax_n, color='none',edgecolor='red', linewidth=0.1)
    plt.title('Minimum of soil Ksat and seepage (ft/day)')
    # plt.savefig(join(fig_dir, 'min_soil_Ksat_seepage.png'), bbox_inches='tight')

# %% [markdown]
# ## recharge potential
# Repeat for steve's recharge potential.

# %%
Kgeom_gdf = gpd.read_file(join(proj_dir, 'GIS', 'Kgeometric_mean_'+grid_id+'.shp'))

# %%
hsg_grid.shape

# %%
# Kgeom_gdf

# %%
rech_est = gpd.read_file(join(gis_dir,'recharge_results','recharge_estimate_'+grid_id+'.shp'))

# %%
rech_est['rch_ft_d'] = rech_est['rch_cm_d']/100/0.3048
rech_est['K_ft_d'] = rech_est.K_m_d/0.3048

# %%
var = 'rch_ft_d'
rech_est[var].quantile([0,.5, 0.95, .99, .999, 1])


# %% [markdown]
# ### unsaturated recharge potential

# %%
fig,ax_n = plt.subplots(figsize=(6.5,5.5))
# log scale doesn't improve things much
var = 'K_ft_d'

vmin = 0.1
vmax = 1
vmin = rech_est[var].min()
vmax = rech_est[var].max()
rech_est.plot(var,  ax = ax_n, legend=True,
             norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax),
              vmin = vmin, vmax = vmax,
              legend_kwds = {'shrink':0.7}
             )

clean_plot(gdf_id, ax_n)
plt.title('Geometric mean of K (ft/day)')
plt.savefig(join(fig_dir, 'Kgeom_unsat_zone.png'), bbox_inches='tight')

# %% [markdown]
# Check how much the geometric mean of conductivity is greater than the soil data
# - cbec might prefer this more but realistically having large sections with 4.5 ft/day is unrealistic but for now is more conservative in keeping sites as potential concepts
# - there are a few spots where the AEM data is missing that would likely be classified as mud rather than coarse based on nearby data but the soil data are used which gives a high recharge rate
#     - need to update these sites by interpolating AEM out slightly or adjusting the soil data manually.
# - tier 3 will be more refined in analyzing sites
#
# When comparing the recharge to the soil Ksat it is even more extreme as the hot spots of sandy soil stick out, this would be fixed by correcting the previous step of soil and AEM ksat.
# - additionally the minimum recharge rate should be between the soil Ksat and recharge because the minimum recharge is 0.17 ft/day from the linear regression but there are many point that fall below that
# - the maximum saturate rate for all cells should be the maximum recharge rate from cells within the domain.

# %% [markdown]
# ### Minimum of Kgeom and soil K

# %% [markdown]
# wherever there isn't AEM data do a slight buffer to determine if it might be reasonable to have coarse, e.g., see if there are any other coarse within 500 m. It might make more sense to manually make this adjustment based on the knowledge of the delta region being lower K and the foothills being lower as well.

# %%
rech_est_comp = rech_est[['Region','K_ft_d']].rename(columns={'K_ft_d':'Kgeom_ft_d'}).copy()
hsg_Kgeom_comp = hsg_grid.merge(rech_est_comp, how='left')


# %%
def adjust_soil_max(hsg_comp, rate1_name, max_dist = 500):
    # this method may not work as well for the parcels
    # separate grid units with and without data
    hsg_filled =  hsg_comp.loc[~hsg_comp[rate1_name].isna(), [rate1_name, 'geometry']]
    hsg_filled = hsg_filled.rename(columns={rate1_name:rate1_name+'_nearby'})
    # use grid cells with data within 500 m
    # sjoin instead of nearest allows more cells to be interected to find nearby maximum
    hsg_filled.geometry = hsg_filled.buffer(max_dist)
    hsg_fill = gpd.sjoin(hsg_comp[hsg_comp[rate1_name].isna()],
        hsg_filled)
    # hsg_fill = gpd.sjoin_nearest(
    #     hsg_comp[hsg_comp[rate1_name].isna()],
    #     hsg_filled, 
    #     max_distance = 500)
    # identify the maximum of the nearby to set a supposed upper limit
    hsg_fill_max = hsg_fill.groupby('Region')[rate1_name+'_nearby'].max()
    
    # fill in NA values with maximum of nearby values
    hsg_comp = hsg_comp.set_index('Region')
    hsg_comp.loc[hsg_fill_max.index, rate1_name] = hsg_fill_max
    hsg_comp = hsg_comp.reset_index()
    return(hsg_comp)


# %%
hsg_Kgeom_comp = adjust_soil_max(hsg_Kgeom_comp, 'Kgeom_ft_d')

# %%
hsg_Kgeom_comp = get_min_rate(hsg_Kgeom_comp, 'K_ft_d', 'Kgeom_ft_d').rename(columns={'lim_rate':'lim_K_ft_d'})

# %%
hsg_Kgeom_comp.area.sum()/hsg_grid.area.sum()

# %%
# hsg_sfr_comp.plot('Qaq_lim', legend=True)
fig,ax_n = plt.subplots()
var = 'lim_K_ft_d'
vmin = 1E-2
vmax = 1E1
# vmin = hsg_Kgeom_comp[var].min()
vmax = hsg_Kgeom_comp[var].max()
hsg_Kgeom_comp.plot(var, ax=ax_n, legend=True,
                 legend_kwds={'shrink':0.7},
                 # norm = mpl.colors.LogNorm(vmin, vmax),
                 )
clean_plot(gdf_id, ax_n)

plt.title('Minimum of soil Ksat and Kgeom (ft/day)')
plt.savefig(join(fig_dir, 'min_soil_Ksat_Kgeom.png'), bbox_inches='tight')

# %%
# save output
hsg_Kgeom_out = hsg_Kgeom_comp[['Region', 'lim_K_ft_d', 'geometry']]
hsg_Kgeom_out = hsg_Kgeom_out.rename(columns={'lim_K_ft_d':'K_ft_d'})
lowK_frac = (hsg_Kgeom_out.K_ft_d<0.12).sum()*100/10E3
print('Only %.2f %% of the min K is less than 0.12 ft/day' %lowK_frac)
# important to note that is becuse there is a homogenizing done with the thicker AEM
# so we see fewer extreme lows which are brought up by the coarse
# hsg_Kgeom_out.plot('K_ft_d')

# %%
hsg_Kgeom_out.to_file(join(fig_dir, 'unsaturated_recharge_rate.shp'))

# %% [markdown]
# ## Saturated recharge potential

# %%
fig,ax_n = plt.subplots(figsize=(6.5,5.5))
# log scale doesn't improve things much
var = 'rch_ft_d'

vmin = 0.1
vmax = 1
vmin = rech_est[var].min()
vmax = rech_est[var].max()
rech_est.plot(var,  ax = ax_n, legend=True,
             # norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax),
              vmin = vmin, vmax = vmax,
              legend_kwds = {'shrink':0.7}
             )

clean_plot(gdf_id, ax_n)
plt.title('30-day average recharge (ft/day)')
plt.savefig(join(fig_dir, 'recharge_30d.png'), bbox_inches='tight')

# %%
peak_rech = rech_est[rech_est.rch_ft_d>1.4]
print('There are %.i' %peak_rech.shape[0], 'cells with more than 1.4 ft/day')

# %% [markdown]
# Note: the minimum possible value of recharge from Steve's results are 5 cm/day (0.17 ft/day).

# %% [markdown]
# ### Minimum of recharge potential and soil Ksat

# %%
rech_est_comp = rech_est[['Region','rch_ft_d']].copy()
# hsg_rech_comp = hsg_grid.merge(rech_est_comp, how='left')
hsg_rech_comp = hsg_Kgeom_out.rename(columns={'lim_K_ft_d':'K_ft_d'}).merge(rech_est_comp, how='left')



# %%
hsg_rech_comp = adjust_soil_max(hsg_rech_comp, 'rch_ft_d', max_dist = 500)

# %%
hsg_rech_comp = get_min_rate(hsg_rech_comp, 'K_ft_d', 'rch_ft_d').rename(columns={'lim_rate':'lim_rch_ft_d'})
# print(hsg_rech_comp.lim_rch_ft_d.max())


# %%

# for soils data set a maximum value for the recharge as the maximum of the recharge
rch_max = rech_est_comp.rch_ft_d.max() # shows extreme hot spots
# the choice of 90th or 99th percentile is based on visual interpolation
# the alternative would be to do a nearest spatial join of recharge again
rch_max = rech_est_comp.rch_ft_d.quantile(0.99)
# adjust where recharge is greater than max where there is soils data
soil_adj = (hsg_rech_comp.lim_rch_ft_d > rch_max)&(hsg_rech_comp.rch_ft_d.isna())
print('%.2f%% soil cells reduced' %(soil_adj.sum()*100/hsg_rech_comp.shape[0]))
hsg_rech_comp.loc[soil_adj, 'lim_rch_ft_d'] = rch_max

# %%
# hsg_sfr_comp.plot('Qaq_lim', legend=True)
fig,ax_n = plt.subplots()
var = 'lim_rch_ft_d'
vmin = 1E-2
vmax = 1E1
# vmin = hsg_rech_comp[var].min()
vmax = hsg_rech_comp[var].max()
hsg_rech_comp.plot(var, ax=ax_n, legend=True,
                 legend_kwds={'shrink':0.7},
                 # norm = mpl.colors.LogNorm(vmin, vmax),
                 )
clean_plot(gdf_id, ax_n)

plt.title('Minimum of soil Ksat and recharge (ft/day)')
plt.savefig(join(fig_dir, 'min_soil_Ksat_recharge.png'), bbox_inches='tight')

# %%
hsg_rech_out = hsg_rech_comp[['Region','lim_rch_ft_d','geometry']]
hsg_rech_out = hsg_rech_out.rename(columns={'lim_rch_ft_d':'rch_ft_d'})


# %%
hsg_rech_out.to_file(join(fig_dir, 'saturated_recharge_rate.shp'))

# %%
