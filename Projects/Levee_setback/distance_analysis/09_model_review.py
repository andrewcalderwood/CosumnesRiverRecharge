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

# %% [markdown]
# # Model review of results
# Following review by Graham and questions on modeling assumptions it is valuable to answer his concerns:
# Do you have justification to use a simple vertical recharge rate at the surface?
# - this is appropriate if there is a disconnected surface from groundwater
#     - answer: quantify volume of recharge and calculate if this would fill the unsat zone
# - the surface is the controlling mechanism of recharge or is the same rate for the lower body of HCPs
#     - answer: refer to steve's work that the geometric of Ksat is most representative, but also refer to the correlation/fit for the surface K. Also look at recharge rates compared to Steve's.
#     - By adding a vertical anisotropy of 100 for all I was essentially doing this but it would likely be better to use geometric mean of K and use Steve's values for Ks.
#
# Use total recharge by model cell to answer the question on recharge rates and the volume of the unsat zone that will be filled.
#
# Notes:
# - steve used a specified head boundary condition of 10 cm above land surface to represent the recharge sites. His vertical gradients from the start will be smaller.

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

# tprogs_id = '' # original tprogs with conditioning data in output tsim
# tprogs_id = '_no_conditioning'
tprogs_id = '_no_cond_c3d'

data_dir = box_dir+ tprogs_id+'/data_output/'
fig_dir = box_dir+tprogs_id+'/figures/'

chan_dir = box_dir+'channel_data/'
gis_dir = chan_dir+'GIS/'

# %% [markdown]
# # Load data
# ## general model data

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
# grid_sfr = gpd.read_file(gwfm_dir+'/SFR_data/final_grid_sfr/grid_sfr.shp')
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')

m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')
sfr_dir = join(gwfm_dir, 'SFR_data')
grid_sfr = gpd.read_file(join(sfr_dir,'final_grid_sfr/grid_sfr.shp'))
# load sacramento river, creeks
rivers = gpd.read_file(join(sfr_dir,'Sac_valley_rivers/Sac_valley_rivers.shp'))
cr = gpd.overlay(rivers.loc[rivers.GNIS_Name=='Cosumnes River'].to_crs('epsg:32610'), m_domain)

# %% [markdown]
# ## setback data

# %%
flow_percentile=6
hf_tot_in =  np.loadtxt(data_dir+'surface_highflow_by_realization_'+str(flow_percentile)+'.tsv',delimiter = '\t')
hf_tot = np.reshape(hf_tot_in, (100, nrow, ncol))
# hf_tot

# %%
setbacks = np.arange(0, 3400,200)
# original XS data
xs_all_cln = pd.read_csv(chan_dir+'Elevation_by_XS_number_meters.csv', index_col='dist_from_center_m')
# smoothed XS data used for setback analysis
# xs_all_cln = pd.read_csv(chan_dir+'xs_levee_smooth.csv', index_col='dist_from_center_m')
num_segs = xs_all_cln.shape[1]
# load array identifying row,col to XS id (1,28)
xs_arr = np.loadtxt(chan_dir+'XS_num_grid_reference.tsv')

# load flood typology characteristics (based on daily data 1908 - 2014) - median values 
#"cms_pk" for peak discharge, "pk_loc" for time to peak, and "log_no_d" for duration
flood_type = pd.read_csv(join(box_dir, 'whipple_grp6_w97ftmedians.csv'),index_col='Group.1')


# %%
soil_thick=2
fn = chan_dir+'/tprogs_geomK_'+str(soil_thick)+'m_depth_dem_mean.tsv'
soil_K_out = np.loadtxt(fn, delimiter='\t')
soil_K = np.reshape(soil_K_out, (100, nrow, ncol))


# %%
import h5py
f = h5py.File(join(chan_dir, 'setback_locs.hdf5'), "r")
local_str_setbacks = f['setbacks']['local'][:]
str_setbacks = f['setbacks']['regional'][:]

f.close()

# %%
region='regional'
ft_in = 1
# def load_hdf5_output(ft_in, region):
tic = time.time()
T_in = int(10**flood_type.loc[ft_in,'log_no_d'])
p_l_in = flood_type.loc[ft_in,'pk_loc']
tp_in = int(p_l_in*T_in)
rch_hf_all = np.zeros((100, len(setbacks),nrow,ncol))

# filter out for only those realizations that successfully ran
base_fn = join(data_dir, region, 'type'+str(ft_in))
r_out = pd.Series(os.listdir(base_fn)).str.extract(r'(\d{3})')[0].unique().astype(int)
# takes a 
for t in r_out[[0]]: # np.arange(0,100): #[0]:
    # load hdf5 files for each realization
    r_fn = join(base_fn,'r'+str(t).zfill(3)+'_')
    f = h5py.File(r_fn+'output.hdf5', "r")
    rch_hf = f['array']['rch_hf'][:]
    # depth is a little complicated to summarize, not so bad to back it out from
    f.close()

    # sum recharge across time to save storage space (breaks python at 25GB)
    # rch_hf = np.nansum(rch_hf, axis=0)
    # rch_hf_all[t] = np.copy(rch_hf)


# convert to m3/day and will have the total recharged after summing individual days
rch_hf_all = rch_hf_all*86400
# convert recharge to m/day
rch_hf_all /= 200*200

toc = time.time()
print('Loading',region,'for flow type',str(ft_in), 'took %.2f minutes' %((toc-tic)/60))
# return(Q_all, rch_hf_all, d_all, d_arr_all, cell_frac_all, d_xs_all)

# %%
print('ndays', T_in, ', Peak day', tp_in)

# %%
# convert from m3/s to m/day
rch_rate = rch_hf * 86400/(200*200)

# %%
# start with intermediate setback (6) then go across all
# want to see recharge rate across time
s = 7
t = tp_in
plt_arr = np.ma.masked_where(rch_rate[tp_in, s, :, :]==0, rch_rate[t, s, :, :])

plt.imshow(plt_arr, norm =mpl.colors.LogNorm())
plt.colorbar(shrink=0.5)


# %% [markdown]
# Very large, 1400m setback, peak flow:
# - The peak recharge rates under the peak day of flow is the max Kvert of 7.5 m/day but most of the spatial distribution is in the 10E-3 to 10E-1 range.
# - The percentiles confirm this as the 50th is 0.9 cm/day and the 75th is 17 cm/day, but the 95th is 80.8 and the 100th is ~750 cm/day which is likely higher than it should be since Steve's max was around 110 cm/day but his inundation depth was only 10 cm.

# %%
def arr_stats(plt_arr, units):
    stat_arr = np.ma.filled(plt_arr, np.nan)
    percentiles = [0,5, 25, 50,75,95,100]
    percs = np.nanpercentile(stat_arr, percentiles)
    print('Percentiles')
    for n in np.arange(0,len(percentiles)):
        print(percentiles[n], 'is %.2E ' %percs[n], units)

arr_stats(plt_arr*100, 'cm/day')

# %% [markdown]
# Thinking through Steve's recharge paper, his discussion really emphasizes it is the interconnected coarse that drive recharge deeper and a deeper pressure response. I think my original model of only comparing recharge in the interconnected coarse was better because then we have a consistent scheme where it is expected that all their recharge will fill the unsat zone and go deeper. This is better because we can expect the unsat zone to fill equally for these sites and the main hydraulic difference should be the flood depth.
# - if I'm able to fact check that the recharge estimates wouldn't fill the unsat zone then we can at least approximate the system isn't greatly impacted by unsat zone filling.
# - the decision is to decide what scale to consider for the unsat zone space (cell by cell, within an HCP?, 25 x 25 like steve's example, floodplain cross-section, or domain-wide in the floodplain)
#     - start cell by cell then 25x25 for simplicity
#     - the mean lengths of 350 m in the x by 1200 m in the y would be another good grouping unit

# %%
ma_arr = np.ma.masked_invalid(rch_rate[:, s, :, :])
tot_rech = np.sum(ma_arr, axis=0)
plt_arr = np.ma.masked_where(tot_rech==0, tot_rech)

plt.imshow(tot_rech, 
           # norm =mpl.colors.LogNorm()
          )
plt.colorbar(shrink=0.5)

# %% [markdown]
# For 1400 m.
# There are 73 cells where the total recharge is greater than 100 m which would be ridiculous.
# - for HCPs would could do a simple vertical accounting based on Steve's work to identify when the unsat zone fills up and then reduce the gradient. This won't work well because reducing the gradient would require specifying the head in the aquifer to the sub meter scale within the soil zone which we wouldn't be able to accurately do.
#
# In contrast, the recharge depth averaged across the setback domain is 6 m.
#
# **Steve's high K site has a cumulative change in storage of about 2500 cm or 25 m in the first year of the 90-day recharge period simulation**
# - and his sites are averaged across 144 cells here so within that there could be even greater peaks. His perturbation snapshots at the 180 day simulation show a max of 40 m increase.

# %%
np.where(plt_arr>100)[0].shape

# %%
arr_stats(plt_arr, 'm')

# %%
tot_rech[str_setbacks[s].astype(bool)].mean()

# %% [markdown]
#  Load old soil data to compare impact of expanding depth on geometric mean of K

# %%
soil_thick = 2
fn = chan_dir+'/tprogs_geomK_'+str(soil_thick)+'m_depth_dem_mean.tsv'
soil_K_out = np.loadtxt(fn, delimiter='\t')
soil_K_2 = np.reshape(soil_K_out, (100, nrow, ncol))


# %%
soil_thick = 10
fn = chan_dir+'/tprogs_geomK_'+str(soil_thick)+'m_depth_dem_mean.tsv'
soil_K_out = np.loadtxt(fn, delimiter='\t')
soil_K_10 = np.reshape(soil_K_out, (100, nrow, ncol))


# %%
r=0
fig,ax = plt.subplots(1,2, figsize=(6.5,4),dpi=300)
im = ax[0].imshow(soil_K_2[r], norm=mpl.colors.LogNorm())
fig.colorbar(im, ax=ax[0], shrink=0.3)
im = ax[1].imshow(soil_K_10[r], norm=mpl.colors.LogNorm())
fig.colorbar(im, ax=ax[1], shrink=0.3)



# %%
plt.imshow(np.std(np.log10(soil_K_2), axis=0))
plt.colorbar()

# %%
K_scale = soil_K_2/ soil_K_10
K_scale[:,~str_setbacks[-1].astype(bool)] = np.nan
r = 60
plt.imshow(K_scale[r],norm=mpl.colors.LogNorm())
plt.colorbar(shrink=0.5)
print('Mean reduction is %i x' %np.nanmean(K_scale[r]))
print('Old max %.2f' %np.nanmax(soil_K_2[r][~str_setbacks[-1].astype(bool)]))
print('new max %.2f' %np.nanmax(soil_K_10[r][~str_setbacks[-1].astype(bool)]))

# %% [markdown]
# For an example realization it looks like the 2 m K is on average about 5x greater than the 10 m. 
