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
from os.path import join, basename,dirname, exists, expanduser
import sys
import glob
import pandas as pd
import numpy as np
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns

# standard geospatial python utilities
# import pyproj # for converting proj4string
# import shapely
import geopandas as gpd
# import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator


# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
    
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study
proj_dir = gwfm_dir + '/Oneto_Denier/'
dat_dir = proj_dir+'Stream_level_data/'

fig_dir = proj_dir+'/Streambed_seepage/figures/'
hob_dir = join(gwfm_dir, 'HOB_data')
sfr_dir = gwfm_dir+'/SFR_data/'




# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

add_path(doc_dir+'/GitHub/flopy')
import flopy 
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)
from mf_utility import get_dates, get_layer_from_elev, clean_wb
from map_cln import gdf_bnds, plt_cln

# from importlib import reload
# import mf_utility
# reload(mf_utility)

# %%

loadpth =  'C:/WRDAPP/GWFlowModel/Cosumnes/Stream_seepage'
model_ws = join(loadpth, 'repeat_hydrology')

    
# model_ws = join(loadpth,'parallel_oneto_denier','realization000')
load_only = ['DIS','BAS6','UPW','SFR','OC', ]
m = flopy.modflow.Modflow.load('MF.nam', model_ws= model_ws, 
                                exe_name='mf-owhm.exe', version='mfnwt',
                              load_only=load_only,
                              )


# %% [markdown]
# The simulation took 4 hr 24 min to run the 5 cycles. It seems that the cycle set up failed because it looks like each year is repeated in a row rather than the cycle of years.

# %%
out_dir = 'restoration_out'
ncycle=5

# %%
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')


# %%
hdobj = flopy.utils.HeadFile(join(model_ws,out_dir,'MF.hds'))


# %%
wb, out_cols, in_cols = clean_wb(join(model_ws, out_dir), dt_ref)
wb_cols = np.append(out_cols, in_cols)
fig,ax= plt.subplots(3,1, sharex=True)
wb.plot(y='PERCENT_ERROR', ax=ax[0])
wb.plot(y=out_cols, ax=ax[1], legend=True)
wb.plot(y=in_cols, ax=ax[2], legend=True)


# %% [markdown]
# The forward project fails to show long-term trends because the GHB outflow shoots up in the (2015 repeat) which drains out all of the water from the system that was built up. Realistically pumping is being curtailed long-term and recharge would be offsetting so we should see less of this extreme outflow.
#

# %%
fig, ax = plt.subplots(3,1, sharex=True)
for n, v in enumerate(['dSTORAGE_sum', 'GHB_NET', 'ET_OUT']):
    wb[v].resample('AS-Oct').mean().plot(ax=ax[n])
    ax[n].set_ylabel(v)
# wb.GHB_NET.resample('AS-Oct').mean().plot(ax=ax[1])
# wb.ET_OUT.resample('AS-Oct').mean().plot(ax=ax[2])

# %%
wb.dSTORAGE_sum.plot()



# %%
import mf_utility
from importlib import reload
reload(mf_utility)
from mf_utility import clean_sfr_df

# %%
# double check the no_reconnection has updated fully
base_ws = join(loadpth,'oneto_denier_upscale4x_2014_2020')
grid_sfr = pd.read_csv(join(base_ws,'grid_sfr.csv'),index_col=0)
grid_sfr = grid_sfr[grid_sfr.strhc1!=0]
# grid_sfr['vka'] = vka[grid_sfr.k, grid_sfr.i, grid_sfr.j]

pd_sfr = grid_sfr.set_index(['iseg','ireach'])[['rchlen','strtop', 'facies', 'strthick', 'slope']]
pd_sfr['Total distance (m)'] = pd_sfr['rchlen'].cumsum()
# if 'Logger Location' in XSg.columns:
#     drop_iseg = XSg[~XSg['Logger Location'].isna()].iseg.values
#     # remove stream segments for routing purposes only
#     grid_sfr = grid_sfr[~grid_sfr.iseg.isin(drop_iseg)]
sfrdf =  clean_sfr_df(join(model_ws,out_dir), dt_ref, pd_sfr, name='MF')
# gradient is stage - Ha/str thick, and strthick=1
sfrdf['h_aquifer'] = -(sfrdf.gradient*sfrdf.strthick - sfrdf.stage)

# %%
sfrdf_last = sfrdf[sfrdf['Total distance (m)']==sfrdf['Total distance (m)'].max()]

# %%
sfrdf_last
plt_strt = pd.to_datetime('2014-10-1')
plt_end = pd.to_datetime('2020-9-30')
fig,ax = plt.subplots()

for n in np.arange(0,ncycle):
    yr_add = pd.DateOffset(years = n*6)
    df = sfrdf_last[plt_strt+yr_add:plt_end+yr_add].copy()
    df = df.sort_values('Qout', ascending=False)
    df['prob'] = np.arange(0, len(df))/(len(df)+1)
    df.plot(x='prob',y='Qout', ax=ax, label=n)
    ax.set_yscale('log')

plt.legend()

# %%


