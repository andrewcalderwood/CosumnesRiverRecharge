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


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

from importlib import reload

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

from mf_utility import get_layer_from_elev
from report_cln import base_round


# %%
import h5py


# %% [markdown]
# # Load data
# 1. Pre-processed AEM to array data
# 2. TPROGs arrays.
# 3. Reference GIS

# %%
grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')


# %%
def elev_2_tprogs_lay(elev):
    """
    Elevation (m) to tprogs layers
    """
    elev = base_round(elev, 0.5)
    lay = int((80 - elev)/0.5)
    return lay
    
grid_sfr['tprogs_lay'] = grid_sfr.z_min.apply(elev_2_tprogs_lay)

# %% [markdown]
# ## AEM

# %%
aem_folder = 'statewide_aem_survey_coarse_fraction_depth_slices_and_average'

with h5py.File(join(upw_dir, aem_folder, 'aem_array.hdf5'), 'r') as f:
    aem_arr = f['facies']['facies_array'][:]
    aem_cf = f['percent_coarse']['pc_array'][:]


# %%
aem_arr_ma = np.ma.masked_where(aem_arr==0,aem_arr)
aem_cf_ma = np.ma.masked_where(aem_cf==0,aem_cf)

k=160
fig,ax = plt.subplots(2,1, sharex=True)
ax[0].imshow(aem_arr_ma[k], cmap='viridis_r')
ax[1].imshow(aem_cf_ma[k], cmap='viridis')
fig.tight_layout(h_pad=0.1)

# %% [markdown]
# ## TPROGs

# %%
mf_tprogs_dir = join(upw_dir, 'tprogs_final')
tprogs_files = glob.glob(mf_tprogs_dir+'/*')


# %%
tprogs_fxn_dir = doc_dir+'/GitHub/CosumnesRiverRecharge/tprogs_utilities'
if tprogs_fxn_dir not in sys.path:
    sys.path.append(tprogs_fxn_dir)
# sys.path
import tprogs_cleaning as tc

# reload(tc)

# %%
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
t=43
tprogs_info = [80, -80, 320]
def load_tprogs(file, dem_data, tprogs_info=[80,-80,320]):
    tprogs_line = np.loadtxt(file)
    masked_tprogs= tc.tprogs_cut_elev(tprogs_line, dem_data, tprogs_info)
    # hide conditioning data
    tprogs_arr = np.where(masked_tprogs<0, masked_tprogs*-1, masked_tprogs)
    return tprogs_arr

tprogs_arr = load_tprogs(tprogs_files[t], dem_data)

# %%
plt.imshow(tprogs_arr[150], cmap = 'viridis_r')


# %% [markdown]
# # Compare TPROGS and AEM
#
# Do a cell by cell comparison and summarize the differences. The greater concern is aligning the coarse where the AEM suggests there is coarse (weighting scheme). It may be better to group into fine/coarse for this analysis because we aren't as worried about getting gravel vs sand as gravel vs mud.

# %%
def cf_comp(aem_arr_ma, tprogs_arr):
    aem_coarse = np.zeros(aem_arr_ma.shape)
    tprogs_coarse = np.zeros(tprogs_arr.shape)
    for f in [1,2]:
        aem_coarse[aem_arr_ma ==f] = 1
        tprogs_coarse[tprogs_arr ==f] = 1
    # mask to avoid comapring no-data
    aem_coarse = np.ma.masked_where(aem_arr_ma.mask, aem_coarse)
    # find where the coarse and fine are matched
    comp = aem_coarse==tprogs_coarse
    return comp
comp = cf_comp(aem_arr_ma, tprogs_arr)

 # %%
 # there is a mix of match and no-match with most match likely from the mud/sandy mud
plt.imshow(comp[150])
plt.colorbar(shrink=0.5)


# %%
def comp_stats(comp, aem_coarse):
    """
    Given the compared arrays return summary stats on matches
    """
    # use the AEM as true/observed and tprogs as the test
    match = 100*(comp).sum()/(aem_coarse>=0).sum()
    coarse_match = 100*(comp[aem_coarse==1]).sum()/(aem_coarse==1).sum()
    
    fine_match = 100*(comp[aem_coarse==0]).sum()/(aem_coarse==0).sum()
    return match, coarse_match, fine_match
    
match, coarse_match, fine_match = comp_stats(comp, aem_coarse)
print('Mathces %.1f %% of AEM cells' %match)
print('Matches %.1f %% coarse' %coarse_match,'and %.2f %% fine' %fine_match)

# %% [markdown]
# As expected the model matches the fines better than the coarse and has poor match of the coarse. Now this should be iterated over all 100 realizations to see how this varies and if some match the coarse better

# %%
sfr_row = grid_sfr.row.values.astype(int)-1
sfr_col = grid_sfr.column.values.astype(int)-1
sfr_lay = grid_sfr.tprogs_lay.values.astype(int)-1

# %%
# sfr_match = 100*(comp[aem_coarse==1]).sum()/(aem_coarse==1).sum()
comp_stats(comp, aem_coarse)


# %%

def subset_arr(arr, row, col, lay=None):
    """Mask non-indexed parts of the array"""
    arr_sub = np.full(arr.shape, np.nan)
    if lay is None:
        arr_sub[:, row, col] = arr[:,row, col]
    else:
        arr_sub[lay, row, col] = arr[lay,row, col]
    arr_sub = np.ma.masked_invalid(arr_sub)
    return(arr_sub)


# %%
# identify aem data below sfr cells only
# aem_sfr = subset_arr(aem_arr_ma, sfr_row, sfr_col)
# aem_coarse_sfr = subset_arr(aem_coarse, sfr_row, sfr_col)
aem_sfr = subset_arr(aem_arr_ma, sfr_row, sfr_col, sfr_lay)
aem_coarse_sfr = subset_arr(aem_coarse, sfr_row, sfr_col,sfr_lay)


# %%
# data frame to save stats on realizations
stats =['all','coarse','fine']
stats_df = pd.DataFrame(columns=stats)
stats_df['r'] = np.tile(np.arange(0,100), 2)
stats_df['grp'] = np.repeat(['all', 'sfr'],100)
stats_df = stats_df.set_index(['r','grp'])

for r in np.arange(0,100):
    # load next tprogs array
    tprogs_arr = load_tprogs(tprogs_files[r], dem_data)
    # perform comparison
    comp = cf_comp(aem_arr_ma, tprogs_arr)
    # assign stats to dataframe
    stats_df.loc[(r,'all')] =  comp_stats(comp, aem_coarse)
    comp_sfr = cf_comp(aem_sfr, tprogs_arr)
    stats_df.loc[(r,'sfr')] = comp_stats(comp_sfr, aem_coarse_sfr)
    print(r, end=', ')
# clean up
stats_df = stats_df.reset_index(level='grp')
stats_df[stats] = stats_df[stats].astype(float)


# %% [markdown]
# A quick test of the first 10 shows that the realizations have a range of 26-30% match of the coarse and a consistent 73% match of the fines. 
# - we may want to get more specific by reviewing the match in the river and near ground surface
# - adding the SFR subset increased the range from 18-33% so we are starting to get a bigger range. Going to the sfr layer explicitly increased the range to 10-39% so it would be best to filter on this level.

# %%
stats_df[stats_df.grp=='sfr'].boxplot(stats)
# stats_df[stats_df.grp=='all'].boxplot(stats)


# %% [markdown]
# There may not be a significant impact of this exercise if the AEM data is highlighting a feature not in the TPROGs model such as the American River channel deposits.
# - it might also be that the best realization isn't the top coarse match but one of the top 10 while having a greater fines match

# %%
# review best SFR matches
# sfr_r = stats_df[stats_df.grp=='sfr'].coarse.idxmax()
sfr_stats = stats_df[stats_df.grp=='sfr'].sort_values('coarse', ascending=False)
sfr_r = sfr_stats.iloc[:10].index.values
# stats_df.loc[sfr_r]

# %%
r=sfr_r[2] # realizations that match the coarse in the channel
r=43 # current realization used based on calibration fit
r = sfr_stats.index[50] # median
print(r)
tprogs_arr = load_tprogs(tprogs_files[r], dem_data)
comp = cf_comp(aem_arr_ma, tprogs_arr)


# %%
aem_arr_dem = tc.get_tprogs_for_elev(aem_arr_ma, dem_data-0.5, dem_data-1, tprogs_info)[0]
tprogs_arr_dem = tc.get_tprogs_for_elev(tprogs_arr_mask, dem_data-0.5, dem_data-1, tprogs_info)[0]
comp_dem = tc.get_tprogs_for_elev(comp, dem_data-0.5, dem_data-1, tprogs_info)[0]
aem_cf_dem = tc.get_tprogs_for_elev(aem_cf, dem_data-0.5, dem_data-1, tprogs_info)[0]

aem_arr_dem = np.ma.masked_where(aem_arr_dem==0, aem_arr_dem)
aem_cf_dem = np.ma.masked_where(aem_cf_dem==0, aem_cf_dem)
tprogs_arr_dem = np.ma.masked_where(aem_arr_dem.mask, tprogs_arr_dem)
comp_dem = np.ma.masked_where(aem_arr_dem.mask, comp_dem)


# %%
fig,ax = plt.subplots(2,2, sharex=True, sharey=True, figsize=(6.5, 3), dpi=300)
ax[0,0].imshow(aem_arr_dem, cmap = 'viridis_r')
ax[0,1].imshow(aem_cf_dem, cmap = 'viridis')
ax[1,0].imshow(tprogs_arr_dem, cmap = 'viridis_r')
ax[1,1].imshow(comp_dem)

# fig.tight_layout(h_pad=0.1)

# %% [markdown]
# - r5 When evaluating the fit visually, it isn't that bad although there are patches in the central region where the AEM would suggest a more continuous high K unit that tprogs intersperses low K units. Is this an acceptable lacking or should we delinate an additional unit?
# - r39 seems to have just as big if not bigger gaps in the middle where the high K would be. Same with r73
# - all realizations lack the high K region toward the foothills to the north east as well, not sure what that could be. Perhaps a sign of fractured bedrock? would need to look at AEM interpretation
