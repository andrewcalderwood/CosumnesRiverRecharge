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
from os.path import basename, dirname, join, exists
import sys
from importlib import reload
import glob
from importlib import reload
import time

# standard python plotting utilities
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# standard geospatial python utilities
import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio



# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

bas_dir = join(gwfm_dir, 'BAS6')
proj_dir = join(gwfm_dir,'EcoFIP')
out_dir = join(proj_dir, 'output')
fig_dir = join(proj_dir,'figures')


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

from map_cln import gdf_bnds, plt_cln
from report_cln import base_round

from flopy_utilities import zone_clean, reach_data_gdf
from mf_utility import get_dates, clean_hob
# from mf_utility import clean_wb # not yet because use GW_IN instead


# %%
run_dir = 'C://WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'
base_model_ws = join(loadpth, model_nam)
# model_nam = 'foothill_vani10'
# model_nam = 'strhc1_scale'
# model_nam = 'sfr_uzf'
# model_nam = 'parallel_realizations/realization005'

model_ws = loadpth+model_nam
print(model_nam)


# %%
load_only = ['DIS','BAS6','UPW','OC','SFR','LAK',
            'RCH', 'WEL'
            ]
m = flopy.modflow.Modflow.load('MF.nam', model_ws= model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only
                              )

if 'LPF' in m.get_package_list():
    gel_nam = 'LPF'
else:
    gel_nam = 'UPW'
gel = m.__getattr__(gel_nam)

# %%
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')
# round of the steady state period
dt_ref['dt'] = dt_ref.dt.dt.round('D')
dt_ref = dt_ref[~dt_ref.steady]

# %% [markdown]
# ## review geology

# %%
# quick check on horizontal to vertical conductivity
vka = gel.vka.array
ibound = m.bas6.ibound.array
hk = gel.hk.array
hk_ma = np.ma.masked_where( ~ibound.astype(bool), hk)

# %%
# quick check of the vertical anisotropy introduced by upscaling
vani = hk/vka
def plt_gel_arr(vani, hk, k):
    fig,ax= plt.subplots(2,1, sharex=True)
    ax[0].set_title('Horizontal Conductivity (m/day)')
    im = ax[0].imshow(hk[k])
    plt.colorbar(im, shrink=0.5)
    ax[1].set_title('Vertical Ansiotropy')
    im = ax[1].imshow(vani[k])
    plt.colorbar(im, shrink=0.5)
    plt.savefig(join(fig_dir, 'geology','layer'+str(k+1)+'_hk_vani.png'), bbox_inches='tight')
    plt.close()
    return None
# k=0
# for k in np.arange(0, m.dis.nlay):
#     plt_gel_arr(vani, hk, k)


# %% [markdown]
# The plot of vertical anisotropy demonstrates that there is a strong vertical anisotropy 100-150 times introduced in cells where there is sand/gravel such that hydraulic conductivity will be more reduced by the fine facies. This is consistent for all layers

# %%
rownum=50

fig, ax = plt.subplots(figsize=(6.5, 3.5)) 

mcs = flopy.plot.PlotCrossSection(model=m, line={'Row' : rownum})

linecollection = mcs.plot_grid(linewidth = 0.3)
ax.add_collection(linecollection)

# need to mask to show land surface
# mcs.plot_array(a=gel.hk.array, norm = mpl.colors.LogNorm())
mcs.plot_array(a=hk_ma, norm = mpl.colors.LogNorm())

# %% [markdown]
# # Load data

# %%
# realizations to present the results for
best10 = pd.read_csv(join(gwfm_dir, 'Regional','top_10_accurate_realizations.csv'))


# %%
sfr_dir = gwfm_dir+'/SFR_data/'
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_p['easting'] = grid_p.geometry.centroid.x
grid_p['northing'] = grid_p.geometry.centroid.y

lak_grid_clip = gpd.read_file(gwfm_dir+'/Levee_setback/lak_grid_clip/lak_grid_clip.shp')


# %%
# load sfr data 
vka = gel.vka.array

grid_sfr = reach_data_gdf(m.sfr, grid_p)
grid_sfr[['row','column']] = grid_sfr[['i','j']] +1 # convert to 1 based to match with SFR output
drop_iseg = grid_sfr[grid_sfr.strhc1==0].iseg.values
grid_sfr['vka'] = vka[grid_sfr.k, grid_sfr.i, grid_sfr.j]
vka_quants = pd.read_csv(join(base_model_ws, 'vka_quants.csv'))
grid_sfr['facies'] = ''
for p in vka_quants.index:
    facies = vka_quants.loc[p]
    grid_sfr.loc[(grid_sfr.vka< facies.vka_max)&(grid_sfr.vka>= facies.vka_min),'facies'] = facies.facies
    # add color for facies plots

# drop routing segments before calculating distances
gdf_sfr = grid_sfr[~grid_sfr.iseg.isin(drop_iseg)]
gdf_sfr = gdf_sfr.set_index(['iseg','ireach'])[['rchlen','strtop','strhc1', 'vka', 'facies', 'geometry']]
gdf_sfr['Total distance (m)'] = gdf_sfr['rchlen'].cumsum()
gdf_sfr['rch_order'] = np.arange(1,len(gdf_sfr)+1) # reach order for easier referencing
pd_sfr = pd.DataFrame(gdf_sfr.drop(columns=['geometry']))
