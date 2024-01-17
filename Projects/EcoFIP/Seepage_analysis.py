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


# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

bas_dir = join(gwfm_dir, 'BAS6')
proj_dir = join(gwfm_dir,'EcoFIP')
plt_dir = join(proj_dir,'figures/')


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
# model_nam = 'foothill_vani10'
model_nam = 'strhc1_scale'
# model_nam = 'sfr_uzf'

model_ws = loadpth+model_nam


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
vka_quants = pd.read_csv(join(model_ws, 'vka_quants.csv'))

for p in vka_quants.index:
    facies = vka_quants.loc[p]
    grid_sfr.loc[(grid_sfr.vka< facies.vka_max)&(grid_sfr.vka>= facies.vka_min),'facies'] = facies.facies
    # add color for facies plots

# drop routing segments before calculating distances
pd_sfr = grid_sfr[~grid_sfr.iseg.isin(drop_iseg)]
pd_sfr = pd_sfr.set_index(['iseg','ireach'])[['rchlen','strtop','strhc1', 'vka', 'facies']]
pd_sfr['Total distance (m)'] = pd_sfr['rchlen'].cumsum()


# %%
from mf_utility import clean_sfr_df
from flopy_utilities import sfr_load_hds

# %%
sfrdf = clean_sfr_df(model_ws, dt_ref, pd_sfr)
# drop routing segments
sfrdf = sfrdf[~sfrdf.segment.isin(drop_iseg)]
# calculate the effective rate of seepage
sfrdf['Qaquifer_rate'] = sfrdf.Qaquifer/(sfrdf.rchlen*sfrdf.width)

# %% [markdown]
# For EcoFIP Tier 2 they would like to have the monthly average seepage for each reach (cell) scaled by the wetted area to provide them a seepage rate. This could be a place where a regression might be helpful in testing how consistent rates and how to best scale seepage by stream stage.

# %%
grp_cols = ['segment','reach']
sfr_mon = sfrdf.groupby(grp_cols).resample('MS').mean(numeric_only=True).drop(columns=grp_cols)
sfr_mon = sfr_mon.reset_index()
sfr_mon['month'] = sfr_mon.dt.dt.month

# %%
r=100
# test the monthly at a reach to look at relationship of flow and seepage
mon_chk = sfr_mon[sfr_mon['Total distance (m)']==pd_sfr['Total distance (m)'].iloc[r]]
fig,ax = plt.subplots(2,1, layout='constrained')
mon_chk.boxplot(by='month',column='Qaquifer', ax=ax[0])
mon_chk.boxplot(by='month',column='Qaquifer_rate', ax=ax[1])



# %%
cols = ['Qaquifer','Qaquifer_rate']
cv = (mon_chk.groupby('month').std()[cols]/mon_chk.groupby('month').mean())[cols]

# %% [markdown]
# The box plot shows there is significant variability due to varying hydrologic conditions, especially streamflow. Now consider if this diminishes after scaling by wetted area.
# - when plotted on a rate scale (m/day) it's clearer that the variablity among water years and seasons is much smaller since it is within an order of magnitude. the coefficient of variation went down as well

# %%

# %% [markdown]
# The flow seepage relationships are not the most helpful because it assumes a relationship of depth to flow that I have baked in. It's safer to give them the relationship of depth to seepage rather than by month. The relationship is almost linear to begin with but a log-log scaling puts more importance on the higher seepage rates.
# - iterate over all the reaches and fit a regression line. Ultimately we may need to give them the daily rates and let them decide how to relate it to their RAS model output.
# - switching to the rate especially created a piece-wise linear relationship which is the Darcy relationship

# %% [markdown]
# ## Output to cbec
# Share a spreadsheet with the segment, reach, row, column, stream top, and polygon so they can spatially reference.
# For a realization share the data in long-format for

# %%
fig,ax = plt.subplots()
mon_chk.plot(x='depth',y='Qaquifer_rate', kind='scatter',ax=ax)
# ax.set_xscale('log')
# ax.set_yscale('log')

# %%
