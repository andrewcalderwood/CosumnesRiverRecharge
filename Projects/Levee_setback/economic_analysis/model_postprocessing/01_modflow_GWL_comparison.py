# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Goal: Review GWL at the end of simulation between scenarios.
# Create contours of baseline and heat map difference with scenarios

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename
import sys
import glob
from importlib import reload
import shutil

import pandas as pd
import numpy as np
from scipy.stats import hmean, gmean

# import calendar
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import geopandas as gpd

# # mapping utilities
# import contextily as ctx
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as fm

# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
gwfm_dir


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

from flopy_utilities import zone_clean, reach_data_gdf

from mf_utility import get_dates, get_layer_from_elev, clean_wb, clean_hob
from map_cln import gdf_bnds, plt_cln
from report_cln import nse, base_round

from map_obs_plt import get_top_active_layer, plot_head_simple


# %%
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

sfr_dir = gwfm_dir+'/SFR_data/'
# grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_p['easting'] = grid_p.geometry.centroid.x
grid_p['northing'] = grid_p.geometry.centroid.y

dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')
lak_grid_clip = gpd.read_file(gwfm_dir+'/Levee_setback/lak_grid_clip/lak_grid_clip.shp')

# %%
run_dir = 'C:/WRDAPP/GWFlowModel'
run_dir = 'F://WRDAPP/GWFlowModel'
# run_dir = 'D://WRDAPP/GWFlowModel'

# loadpth = run_dir +'/Cosumnes/levee_setback/streamflow/'
# # model_nam = 'setback_streamflow'
# model_nam = 'historical_streamflow'

loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology_reconnection'
# model_nam = 'input_write_2014_2020'
model_nam = 'input_write_2014_2022'

base_model_ws = loadpth+model_nam

# %%
load_only = ['DIS','BAS6','SFR']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm', version='mfnwt',
                              load_only = load_only,
                              )

nrow,ncol,nlay,delr,delc = (m.dis.nrow, m.dis.ncol, m.dis.nlay, m.dis.delr, m.dis.delc)

# if 'LPF' in m.get_package_list():
#     gel_nam = 'LPF'
# else:
#     gel_nam = 'UPW'
# gel = m.__getattr__(gel_nam)

# %%
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')


# %%
dt_ref['totim'] = (dt_ref.dt-dt_ref.dt.iloc[0]).dt.days+1


# %%
sfr = m.sfr
# load sfr data 
grid_sfr = reach_data_gdf(sfr, grid_p)

# %% [markdown]
# # Iterate over different economic runs
#
# Goal: 
# - check change in storage is consistent, i.e., overall mass balance is consistent
# - check that pumping dynamics are consistent on seasonal and annual basis, timing will be different due to daily vs timed irrigation periods
# - recharge may be the most affected due to larger irrigation events and also because the SWB is different between native and non-native and that irrigated lands have a simpler SWB run in the winter. Currently soil moisture is reset to field capacity (I think, maybe wilting point) instead of passing the actual soil moisture from previous run end.
#
#

# %%
loadpth = run_dir +'/Cosumnes/Economic/'
model_nam = 'input_write_2014_2022'
scen = 'R200'
# it is probably better to create a slightly different file name then to copy these over for a set scenario
econ_model_ws = join(loadpth, model_nam+'_'+scen, 'crop_modflow')

all_run_dates = pd.read_csv(join(econ_model_ws, 'all_run_dates.csv'))



# %%
# load_only = ['DIS']
# m_per = flopy.modflow.Modflow.load('MF.nam', model_ws= join(econ_model_ws,d), 
#                                 exe_name='mf-owhm', version='mfnwt',
#                               load_only = load_only,
#                               )
# it looks like the start date time isn't correct in the economic model ws so can't really use this
# strt_date_per, end_date_per, dt_ref_per2 = get_dates(m_per.dis, ref='strt')


# %%
# wb_per_all = pd.DataFrame()
for n, d in enumerate(all_run_dates.date[:-1]):
    # create locally referenced dt_ref file to avoid reloading dis
    dt_ref_per = dt_ref[(dt_ref.dt<=all_run_dates.date[n+1])&(dt_ref.dt>=d)].copy()
    dt_ref_per.kstpkper = list(zip(np.zeros(len(dt_ref_per),dtype=int), np.arange(len(dt_ref_per))))

    econ_ws_yr =join(econ_model_ws, d)
    # wb_per, out_cols_per, in_cols_per  = clean_wb(econ_ws_yr, dt_ref_per)
    # wb_per_all = pd.concat((wb_per_all, wb_per))
# econ_ws_yr

# %%
# pretty slow to load (<1 min though)
hdobj = flopy.utils.HeadFile(base_model_ws+'/MF.hds')
spd_stp = hdobj.get_kstpkper()
times = hdobj.get_times()
# cbc = base_model_ws+'/MF.cbc'

# %%
def plt_hob_map(head, y, s, hob=True, nd_chk=None, rch=False, contour=False, hk=False, step=15):
    fig,ax=plt.subplots(figsize=(8, 8))
    m_domain.plot(ax=ax,color='None')
    mapview = flopy.plot.PlotMapView(model=m,ax=ax)
    if hk:
        csa = mapview.plot_array(gel.hk.array[1,:,:], norm=mpl.colors.LogNorm())
        cb = plt.colorbar(csa, shrink=0.5,ax=ax)
        cb.set_label('Horiz. Cond. (m/d)')

#     csa = mapview.plot_array(-pump_rate.mean(axis=0), vmin=-np.quantile(pump_rate, .95))
    if rch:
        csa = mapview.plot_array(rech.mean(axis=0)- pump_rate.mean(axis=0), vmin=-np.quantile(pump_rate, .95))
        cb = plt.colorbar(csa, shrink=0.5,ax=ax)
        cb.set_label('Recharge (m/d)')
    if hob:
        # code where year and season were specified, this will get simplified
        hob_gpd_plt = hob_gpd[(hob_gpd.index.year==y)&(hob_gpd.season==s)]

    # stns[stns.botm_elev > stns.screen_elev].plot(color='red',marker='x',ax=ax)
    grid_sfr.plot(ax=ax,color='black')
    if contour:
        head_ma = np.ma.masked_where(head==-999.99, head)  
        # pull water table and layer below
        head_ma = get_top_active_layer(head_ma)
        hmin, hmax = base_round(head_ma.min(), step), base_round(head_ma.max(), step)
        levels = np.arange(hmin, hmax, step)
        contour_set = mapview.contour_array(head_ma,
                                    masked_values=[-999.99], levels=levels, ax=ax)
        hcb = plt.colorbar(contour_set, shrink = 0.5,ax=ax)
        hcb.set_label('Head (m)')
        ax.clabel(contour_set, contour_set.levels[0::], inline=True, fontsize=8)
#     foothills.plot(ax=ax, alpha=0.5, edgecolor='black', color='grey')
    if nd_chk != None:
        hob_gpd_plt = hob_gpd_plt[hob_gpd_plt.node.isin(nd_chk)]
    if hob:
        # hob_gpd.plot('error',scheme='EqualInterval', k= 6, ax=ax,legend=True,cmap='magma')
        hob_gpd_plt.plot('error',markersize='abs_error',scheme='Quantiles', k = 6, ax=ax,
                          legend=True,cmap='bwr_r',legend_kwds={'loc':(1.1,0.9),'title':'Error (Sim - Obs)'})
        hob_gpd_plt.apply(lambda x: ax.annotate(str(x.node), xy=(x.geometry.x, x.geometry.y), ha='right'),axis=1);
    
        gdf_bnds(hob_gpd_plt, ax=ax, buf=2E3)
        return(hob_gpd_plt)
    plt_cln(ax=ax)
    return None

    # ax.legend(loc=(1,0.5))


# %%
scen = 'R204'
# it is probably better to create a slightly different file name then to copy these over for a set scenario
econ_model_ws = join(loadpth, model_nam+'_'+scen, 'crop_modflow')


# %%
# all_run_dates.iloc[[-2]].date
# dt_ref_per.dt
# dt_ref_per.dt

# %%

# obs_ts_df_all = pd.DataFrame()
for n, d in enumerate(all_run_dates.date[:-1]):
# for n, d in enumerate(all_run_dates.iloc[[-2]].date):
    # create locally referenced dt_ref file to avoid reloading dis
    dt_ref_per = dt_ref[(dt_ref.dt<=all_run_dates.date[n+1])&(dt_ref.dt>=d)].copy()
    dt_ref_per.kstpkper = list(zip(np.zeros(len(dt_ref_per),dtype=int), np.arange(len(dt_ref_per))))
    dt_ref_per['totim'] = (dt_ref_per.dt-dt_ref_per.dt.iloc[0]).dt.days+1
    
    econ_ws_yr =join(econ_model_ws, d)
    hdobj_yr = flopy.utils.HeadFile(econ_ws_yr+'/MF.hds')


# %%
# extract heads from last timestep of each year, results show last doesn't exist
# so use second to last
use_per = dt_ref_per.iloc[-2]
hds = hdobj_yr.get_data(kstpkper= use_per.kstpkper)


# %%
plt_hob_map(hds, use_per.loc['dt'].year, 'fall', hob=False, nd_chk=None, rch=False, contour=True, hk=False, step=5)

# %%
