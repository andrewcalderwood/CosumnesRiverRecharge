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
proj_dir = join(gwfm_dir,'Mapping')
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

# from mf_utility import get_layer_from_elev
# from map_cln import gdf_bnds, plt_cln
import map_obs_plt as mop
from map_obs_plt import plt_bc_hk, plot_head_simple, plot_dtw_simple
from map_cln import gdf_bnds, plt_cln
from report_cln import base_round

from flopy_utilities import zone_clean, reach_data_gdf
from mf_utility import get_dates, clean_hob
# from mf_utility import clean_wb # not yet because use GW_IN instead


# %%
import regional_utilities
reload(regional_utilities)
from regional_utilities import clean_wb

# %%
sfr_dir = gwfm_dir+'/SFR_data/'
# grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_p['easting'] = grid_p.geometry.centroid.x
grid_p['northing'] = grid_p.geometry.centroid.y

dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')
lak_grid_clip = gpd.read_file(gwfm_dir+'/Levee_setback/lak_grid_clip/lak_grid_clip.shp')


# %%
run_dir = 'C://WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'
base_model_ws = join(loadpth, model_nam)
# model_nam = 'historical_geology_cal'
# model_nam = 'historical_extended'
model_nam = 'input_write_2000_2022'

model_ws = loadpth+model_nam


# %%
load_only = ['DIS','BAS6','UPW','OC','SFR','LAK',
            'RCH',
             'WEL'
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
# gel = flopy.modflow.ModflowUpw.load(join(model_ws, 'MF.upw'), model=m)
# sfr = flopy.modflow.ModflowSfr2.load(join(model_ws, 'MF.sfr'), model=m)

# %%
sfr = m.sfr
vka = gel.vka.array
# load sfr data 
grid_sfr = reach_data_gdf(sfr, grid_p)
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
all_obs = pd.read_csv(model_ws+'/input_data/all_obs_grid_prepared.csv',index_col=0, parse_dates=['date'])
all_obs.index = all_obs.index.rename('date')
all_obs = all_obs.reset_index()


# %% [markdown]
# # Head plots and contours

# %%
hdobj = flopy.utils.HeadFile(model_ws+'/MF.hds')
spd_stp = hdobj.get_kstpkper()
times = hdobj.get_times()
cbc = model_ws+'/MF.cbc'

# strt_date = pd.to_datetime(m.dis.start_datetime)
# dates = strt_date+(np.asarray(times)-1).astype('timedelta64[D]')

# dt_ref = pd.DataFrame(dates, columns=['dt'])
# dt_ref['kstpkper'] = spd_stp
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')
# round of the steady state period
dt_ref['dt'] = dt_ref.dt.dt.round('D')
dt_ref = dt_ref[~dt_ref.steady]

# %%
rech = m.rch.rech.array[:,0,:,:]
rech_avg = pd.Series(rech.mean(axis=(1,2)))[1:]
# rech_avg.index=dt_ref.dt[:-1]

# %%
# convert pumping to array
pump = np.zeros((m.dis.nper,m.dis.nrow,m.dis.ncol))
for n in np.arange(0,m.dis.nper):
    wel_n = m.wel.stress_period_data[n]
    pump[n, wel_n.i, wel_n.j] += wel_n.flux*-1
pump_rate = pump/(m.dis.delr[0]*m.dis.delc[0])

# %%
rch_total = rech.sum(axis=(1,2))
pump_total = pump_rate.sum(axis=(1,2))

plt.plot(rch_total.cumsum(), label='rch')
plt.plot(pump_total.cumsum(), label='pump')
plt.legend()

# %%
# plt.imshow(rech.mean(axis=0), vmax=0.001)
# # plt.imshow(pump_rate.mean(axis=0), vmax=0.005)

# plt.colorbar(shrink=0.6)
# rech.mean(axis=0).sum()

# %% [markdown]
# ## Water Budget check

# %%
wb, out_cols, in_cols = clean_wb(model_ws, dt_ref)

# %%
wb_ann = wb.resample('AS-Oct').sum(numeric_only=True)
fig,ax = plt.subplots( sharex=True)
plt.axhline(0, color='black')
af_scale = 1/(0.3048**3)/43560/1000

wb_ann[out_cols].multiply(-1).multiply(af_scale).plot( kind='bar', ax=ax, stacked=True)
wb_ann[in_cols].multiply(af_scale).plot( kind='bar', ax=ax, stacked=True)
plt.legend(loc=(1.05,0.3))
plt.ylabel('TAF/year')

# %%
fig,ax = plt.subplots(3,1, figsize=(6.5,4), sharex=True)
wb.plot(y='PERCENT_ERROR', ax=ax[0])
wb[in_cols].multiply(1E-6).plot(y=in_cols, ax=ax[1], legend=True)
wb[out_cols].multiply(1E-6).plot(y=out_cols, ax=ax[2], legend=True)

ax[1].set_ylabel('Inflow\n($10^6 m^3/day$)')
ax[2].set_ylabel('Outflow\n($10^6 m^3/day$)')

# ax[1].set_ylim(0,1E3)
# ax[2].set_ylim(0,1E5)

# %%

fig,ax= plt.subplots(2,1, figsize=(6.5, 4), sharex=True, dpi=300)
# in_labels = ['Foothills GW Inflow','Boundary GW Inflow','GW Recharge', 'Stream Leakge to GW']
# wb[in_cols].multiply(1E-6).plot(y=in_cols, ax=ax[0], legend=True, label=in_labels)
# out_labels = ['Foothills GW Outflow','GW Pumping','Boundary GW Outflow','GW Leakage to Streams']
# wb[out_cols].multiply(1E-6).plot(y=out_cols, ax=ax[1], legend=True, label=out_labels)

in_labels = ['GW Recharge', 'Stream Leakge to GW','GW Inflow']
wb[in_cols].multiply(1E-6).plot(y=['RCH_IN','SFR_IN', 'GW_IN'], ax=ax[0], legend=True, 
                                label=in_labels, color=['tab:green','tab:blue','tab:brown'])
out_labels = ['GW Pumping','GW Outflow', 'GW ET']
wb[out_cols].multiply(1E-6).plot(y=['WEL_OUT','GW_OUT','ET_OUT'], ax=ax[1], legend=True, 
                                 label=out_labels, color=['black','tab:brown', 'green'])

ax[0].set_ylabel('Inflow\n($10^6 m^3/day$)')
ax[1].set_ylabel('Outflow\n($10^6 m^3/day$)')
plt.xlabel('Date')
# plt.savefig(join(plt_dir, 'total_water_budget_time_series.png'),  bbox_inches='tight')


# %% [markdown]
# # Sim vs Obs Head
#

# %%
# model_ws = loadpth+'historical_simple_geology'
# model_ws = loadpth+'historical_simple_geology_reconnection'


# %%
hobout = clean_hob(model_ws, dt_ref, split_c = '.', obs_id_spd=False).rename(columns={'Sensor':'obs_site'})


# %%
# temporary add on to drop bad msmts here
hobout = hobout[hobout.obs_nam.isin(all_obs.obs_nam)]

fig, ax = plt.subplots(1,1,figsize=(5,5))

# get boundary values for plotting a 1:1
hobmax = hobout.loc[:,['sim_val','obs_val']].max().min()
hobmin = hobout.loc[:,['sim_val','obs_val']].min().max()

hob_lin = np.array([hobmin, hobmax])

# plot observed vs simulated values
hobout.plot.scatter(x='obs_val', y='sim_val',marker='.',ax=ax)
ax.plot(hob_lin,hob_lin,'red')
# plot buffer lines with 30 ft, 10 ft limit
hob_lin_adj = 10*0.3048
ax.plot([hobmin-hob_lin_adj/2, hobmax-hob_lin_adj/2], [hobmin+hob_lin_adj/2, hobmax+hob_lin_adj/2], 'gray')
ax.plot([hobmin+hob_lin_adj/2, hobmax+hob_lin_adj/2], [hobmin-hob_lin_adj/2, hobmax-hob_lin_adj/2], 'gray')

ax.set_xlabel('Observed Values (m)')
ax.set_ylabel('Simulated Values (m)')
# plt.xlabel('Observed Values (m)')
# plt.ylabel('Simulated Values (m)')

# lim2 = hobout.loc[:,['obs_val']].max().min()
# lim1 = hobout.loc[:,['obs_val']].min().max()
# ax.set_ylim(lim1,lim2)

fig_nam = plt_dir+'GSP_WaterBudget/sim_vs_obs_heads'

# plt.savefig(fig_nam+'.png',dpi=600,bbox_inches='tight')
# plt.savefig(fig_nam+'.svg',dpi=600,bbox_inches='tight')

# %%
def mak_hob_gpd(hobout, all_obs):
    # join more indepth obs data to output simulated heads
    obs_data = hobout.join(all_obs.set_index('obs_nam'),on=['obs_nam'], how='inner')
    obs_data = obs_data.dropna(subset=['node'])
#     obs_data.loc[:,['row','column','node']] = obs_data.loc[:,['row','column','node']].astype(int)
    obs_data[['row','column','node']] = obs_data[['row','column','node']].astype(int)
    # add hk to plot unit
    obs_data['hk'] = gel.hk.array[obs_data.layer-1, obs_data.row-1, obs_data.column-1]
    # obs_data.index = obs_data.index
    obs_grid = obs_data.join(grid_p.set_index(['row','column']).loc[:,['easting','northing']], 
                             on=['row','column'])
    # # convert back to geospatial
    hob_gpd = gpd.GeoDataFrame(obs_grid, geometry = gpd.points_from_xy(obs_grid.easting, obs_grid.northing),
                              crs = grid_p.crs)
    hob_gpd['error'] = hob_gpd.sim_val - hob_gpd.obs_val 
    hob_gpd['abs_error'] = hob_gpd.error.abs()
    # add error statistics
    hob_gpd['Statistic'] = 0.01
    hob_gpd['StatFlag'] = 'SD'
    # locations with significant difference between RPE GSE and the DEM should have additional uncertainty included
    hob_gpd['Statistic'] += np.round(np.abs(hob_gpd.dem_wlm_gse),4)
    hob_gpd['Weight'] = 1/(hob_gpd.Statistic**2)
    
    if 'date' in hob_gpd.columns:
        hob_gpd = hob_gpd.set_index('date')
        hob_gpd.index = pd.to_datetime(hob_gpd.index)    
        #     groupby values by season
        hob_gpd.loc[(hob_gpd.index.month > 2)&(hob_gpd.index.month < 6),'season'] = 'spring'
        hob_gpd.loc[(hob_gpd.index.month > 8)&(hob_gpd.index.month < 12),'season'] = 'fall'

    # simplify to stations for reference
    stns = hob_gpd.drop_duplicates('site_code', keep='last').reset_index().drop(columns=['date','gwe'])
    stns['botm_elev'] = m.dis.botm[stns.layer-1, stns.row-1, stns.column-1]
    stns.crs = hob_gpd.crs
    
    return(hob_gpd, stns)
    # set date
    


# %%
hob_gpd, stns = mak_hob_gpd(hobout, all_obs)

hob_seasonal = hob_gpd.groupby(['node','season']).mean(numeric_only=True)
hob_seasonal = gpd.GeoDataFrame(hob_seasonal, geometry = gpd.points_from_xy(hob_seasonal.easting, hob_seasonal.northing))
hob_seasonal = hob_seasonal.reset_index()
    


# %%


soswr = (np.sum(np.abs(hob_gpd.sim_val-hob_gpd.obs_val)*hob_gpd.Weight))
print('Sum of absolute difference of OBS and SIM: %.2e' %soswr)

from sklearn.metrics import mean_squared_error
from report_cln import nse

rmse = mean_squared_error(hob_gpd.obs_val, hob_gpd.sim_val, squared=False)
print('Root mean square error is %.2f m' %rmse)
nse_out = nse(hob_gpd.obs_val, hob_gpd.sim_val)
print('NSE is %.2f' %nse_out)

# %%
from report_cln import base_round
from map_cln import plt_cln


# %%
def get_top_active_layer(head_ma):
    """ Sample the top active value for a 3d array for each row,column"""
    if head_ma.mask.any():
        head_loc = pd.DataFrame(np.transpose(np.where(~head_ma.mask)), columns=['k','i','j'])
        head_loc = head_loc.groupby(['i','j']).min().reset_index()
        # top active value for each row,column
        head_top = np.full((head_ma.shape[1], head_ma.shape[2]), np.nan)
        head_top[head_loc.i, head_loc.j] = head_ma[head_loc.k, head_loc.i, head_loc.j]
        head_top = np.ma.masked_invalid(head_top)
    else:
        # if nothing is masked then the first layer is the maximum
        head_top = head_ma[0,:,:]
    return head_top


# %%
def plt_hob_map(y, s, hob=True, nd_chk=None, rch=False, contour=False, hk=False, step=15):
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

    hob_gpd_plt = hob_gpd[(hob_gpd.index.year==y)&(hob_gpd.season==s)]

    # stns[stns.botm_elev > stns.screen_elev].plot(color='red',marker='x',ax=ax)
    grid_sfr.plot(ax=ax,color='black')
    if contour:
        head = hdobj.get_data((0,int(hob_gpd_plt.spd.mean())))
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
# plot plain contours for reference
hob_gpd_plt = plt_hob_map(2019, 'fall', hob=True, rch=False, contour=True, hk=False, step=5)
# hob_gpd_plt = plt_hob_map(2016, 'fall', hob=True, rch=False, contour=True, hk=False, step=5)

# %%

# %%
# nd_chk = [15343, 16733, 11448, 8437, 15314, 14626] +[3103, 5642, 6112, 10746, 6458]
# # nd_chk = [2926, 8437, 12944, 13407]
# nd_chk = [15314, 15343, 13407, 12944, 14626]
# nd_chk = [6458, 8437, 9580, 11448, 15314, 20055]
# # nd_chk = [6458, 8437, 9580, 10884, 11448]
nd_chk = [20055, 16614, 22825] # southeast boundary
# nd_chk = [10161, 10165, 10383, 11078, 11084] # Oneto-Denier
hob_gpd_chk = plt_hob_map(2016, 'fall', nd_chk=nd_chk, rch=False, contour=True)


# %%
# df_chk = hob_gpd_chk.groupby('node').mean(numeric_only=True)
# df_chk['rech'] = rech.sum(axis=0)[df_chk.row.astype(int)-1, df_chk.column.astype(int)-1]
# df_chk[['sim_val','obs_val','avg_screen_depth', 'hk', 'rech', 'layer','abs_error']]

# %%
hobout = clean_hob(model_ws, dt_ref, split_c = '.').rename(columns={'Sensor':'obs_site'}).drop(columns=['spd'])

hob_gpd, stns = mak_hob_gpd(hobout, all_obs)
# simplify dataset to time series 
hob_long = hob_gpd.melt(value_vars=['sim_val','obs_val'], id_vars=['node'], ignore_index=False)

# # a few wells have duplicates in a node (same with site_code), temp fix
# issue was actually the NA values
# hob_long = hob_long.reset_index().drop_duplicates(['date','node','variable']).set_index('date')


# %%
# # ## load simulated output from one alternate scenario for comparison
# alt_ws = join(loadpth, 'historical_simple_geology_reconnection')
# hobout_alt = clean_hob(alt_ws, dt_ref, split_c = '.').rename(columns={'Sensor':'obs_site'}).drop(columns=['spd'])
# hobout_alt = hobout_alt.rename(columns={'sim_val':'sim_alt'}).drop(columns=['obs_val'])
# hob_gpd_alt = hob_gpd[['sim_val','obs_val', 'obs_nam','node']].reset_index().merge(hobout_alt).set_index('date')
# hob_long_alt = hob_gpd_alt.melt(value_vars=['sim_val','obs_val', 'sim_alt'], id_vars=['node'], ignore_index=False)


# %% [markdown]
# When presenting data for clients we should use large axis scale to represent what we are interested in which is regional dynamics

# %%
# chk_lng = hob_long.groupby('node').count().variable
# chk_lng = chk_lng[chk_lng>50].index
# sns.relplot(hob_long[hob_long.node.isin(chk_lng)], x='date',y='value', 

# sns.relplot(hob_long_alt.dropna(subset='value'), x='date',y='value', 
sns.relplot(hob_long.dropna(subset='value'), x='date',y='value', 
            hue='variable',  col='node',
#             col_wrap=10, # for powerpoint
            col_wrap=4,
           facet_kws={'sharex':True, 'sharey':False}
           )


# %%
# nd_chk = 21451
# val = -30
# hob_gpd[hob_gpd.node==nd_chk][hob_gpd[hob_gpd.node==nd_chk].obs_val<val]

# %%
# plt.imshow(rech.mean(axis=0))

# %%
# for n in [10161, 10383, 11078, 11084]:
#     site = stns[stns.node==10161].iloc[0]
#     # always inactive hob
#     # plt.plot(hdobj.get_ts((site.layer-1, site.row-1, site.column-1)))
#     k = site.layer-1
#     n_ib = m.bas6.ibound.array[k, site.row-1, site.column-1]
#     print('ibound', n_ib, 'layer', k+1)

# %%
# all NA values for the Oneto-Denier monitoring wells after increasing floodplain recharge

# hob_gpd[hob_gpd.node==10161]

# %% [markdown]
# ### review observation details

# %%
node = [16733, 15314, 15343]
# node=[6085, 6458, 20055]
stn_chk = stns[stns.node.isin(node)].copy()
cols = ['site_code','node','agency','wlm_gse', 'dem_elev', 'layer', 'hk', 'row','column']
cols = ['node','wlm_gse','top_prf_int','bot_prf_int', 'interp_depth','screen_elev' ]
stn_chk[cols]
# # stn_chk.columns

# %% [markdown]
# # general contour check

# %%
ghb_dir = join(gwfm_dir, 'GHB_data')
year = 2010 # 2016
filename = glob.glob(ghb_dir+'/final_WSEL_arrays/fall'+str(year)+'_kriged_WSEL.tsv')[0]
# convert from ft to meters
hd_strt = np.loadtxt(filename)*0.3048
# extent = (minx,maxx,miny,maxy)

# %%
fig,ax=plt.subplots(figsize=(8, 8))
m_domain.plot(ax=ax,color='None')
mapview = flopy.plot.PlotMapView(model=m,ax=ax)

step=2
levels = np.arange(base_round(hd_strt.min(),step), base_round(hd_strt.max(),step), step)

cs = mapview.contour_array(hd_strt, levels=levels)
# CS = plt.contour(gwl_arr, levels= breaks, extent = extent, colors='blue')

# need to do this after calling zoom to avoid negative impacts
ax.clabel(cs, cs.levels, inline=True, fmt="%2.0f", fontsize=10) #fmt
plt_cln(ax=ax)

# %%
y=2016
s='fall'
filename = glob.glob(ghb_dir+'/final_WSEL_arrays/'+s+str(y)+'_kriged_WSEL.tsv')[0]
# convert from ft to meters
hd_strt = np.loadtxt(filename)*0.3048

row = 50
hob_gpd_plt = hob_gpd[(hob_gpd.index.year==y)&(hob_gpd.season==s)]
head = hdobj.get_data((0,int(hob_gpd_plt.spd.mean())))
head_ma = np.ma.masked_where(head==-999.99, head)  
hmin, hmax = base_round(head_ma.min(), step), base_round(head_ma.max(), step)
levels = np.arange(hmin, hmax, step)
# calculate heads in the water table
# head_wt = np.max(head_ma[:-3], axis=0)
head_wt = np.max(head_ma[0:-2], axis=0)

plt.plot(head_wt[row,:],label='Simulation - WT')
# plt.plot(head_ma[-2][row,:],label='Simulation Laguna')
# plt.plot(head_ma[-1][row,:],label='Simulation Mehrten')
plt.plot(hd_strt[row,:], label='Contour')
plt.plot(dem_data[row,:],label='Ground surface',color='black',linestyle='--')

plt.legend()
print('simulated min %.2f' %head_ma[0][row,:].min(),'and observed %.2f' %hd_strt[row,:].min())


# %% [markdown]
# ## cross-section

# %%
ghb_dir = join(gwfm_dir, 'GHB_data')
year = strt_date.year # 2016
filename = glob.glob(ghb_dir+'/final_WSEL_arrays/spring'+str(year)+'_kriged_WSEL.tsv')[0]
# convert from ft to meters
hd_strt = np.loadtxt(filename)*0.3048

# %%
# model_nam = 'historical_simple_geology_reconnection'
# model_nam = 'foothill_vani10'

# model_ws = loadpth+model_nam


# %%
ibound = m.bas6.ibound.array
hk = m.upw.hk.array
hk_ma = np.ma.masked_where( ~ibound.astype(bool), hk)

# %%
lakarr = m.lak.lakarr.array[0]
lakarr = np.where(lakarr==0, np.nan, lakarr)

# %%
reach_data = pd.DataFrame(m.sfr.reach_data)
# reach_data.set_index(['k','i','j']).loc[chk]
# neither cell with convergence issues had a stream directly in it
# reach_data[(reach_data.j>40)&(reach_data.j<55)]
# reach_data[(reach_data.j>0)&(reach_data.j<20)]

# %%
# sfr array for plotting
sfr_arr = np.full(hk[0].shape,np.nan)
sfr_arr[reach_data.i, reach_data.j] = reach_data.strhc1

# %%
# # chk = (10,43, 41)
# chk_all = [(10,41, 43), (10,44, 46),(16,47,14)]
# # chk = (10,46, 44)
# # chk = (16,14,47)
# chk = chk_all[2]
# lakarr[chk], hk[chk]
# fig,ax = plt.subplots(dpi=600)
# # plt.plot(chk_all)
# # plt.imshow(hk[10])
# plt.imshow(lakarr[7])
# ax.imshow(sfr_arr, norm = mpl.colors.LogNorm())
# for n in chk_all:
#     ax.scatter(n[2], n[1], marker='x')

# %%
hdobj = flopy.utils.HeadFile(model_ws+'/MF.hds')
fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=300) 
# plt.aspect=10
plt.aspect=5 # better when zooming in to -100 to 100

rownum = 50
sp1 = 100
sp2 = 1400
# sp1 = 58
# sp2 = 100

mcs = flopy.plot.PlotCrossSection(model=m, line={'Row' : rownum})
# mcs = flopy.plot.PlotCrossSection(model=m, line={'Column' : 50})

linecollection = mcs.plot_grid(linewidth = 0.3)
ax.add_collection(linecollection)

# need to mask to show land surface
mcs.plot_array(a=hk_ma, norm = mpl.colors.LogNorm())

wt = mcs.plot_surface(a=hd_strt[:,:], color='black')

head = hdobj.get_data(kstpkper = spd_stp[sp1])
wt = mcs.plot_surface(a=head[:,:,:], color='blue')

head_new = hdobj.get_data(kstpkper = spd_stp[sp2])
wt = mcs.plot_surface(a=head_new[:,:,:],color='red')
plt.ylim(-100, 100)

plt.xlabel('Distance from southwestern edge (m)')
plt.ylabel('Elevation (m)')

# %% [markdown]
# A simple 3 layer model runs in 25 minutes while a refined 30ish layer model with 4m layers from 40 to -60 takes 4.5 hours. But most of the water table fluctuations is between 0 m and -20 m. 
# - for calibration runs we could do refined geology from 10 m to -30 m which cut out 30 m of complex saturated flow.
# - I should just try the calibration with the 4.5 hour model first and see if it converges within a couple of days.

# %% [markdown]
# # Plot stream water budget

# %%
model_ws

# %%
from mf_utility import clean_sfr_df
from flopy_utilities import sfr_load_hds

# %%
sfrdf = clean_sfr_df(model_ws, dt_ref, pd_sfr)
hdobj = flopy.utils.HeadFile(model_ws+'/MF.hds')
# alt_ws = join(loadpth, 'fine_tprogs')
# sfrdf = clean_sfr_df(alt_ws, dt_ref, pd_sfr)
# hdobj = flopy.utils.HeadFile(alt_ws+'/MF.hds')

# drop routing segments
sfrdf = sfrdf[~sfrdf.segment.isin(drop_iseg)]

# shouldn't put in function directly as this depends on local units
sfrdf['Qout_cfs'] = sfrdf.Qout/(86400*0.3048**3)

# %%
# plt_dates = ['2016-1-1']
plt_dates = pd.date_range('2017-1-1','2017-4-1')
sfr_heads, avg_heads = sfr_load_hds(hdobj, grid_sfr[~grid_sfr.iseg.isin(drop_iseg)], plt_dates, dt_ref)

# %%
var = ['Qrech','Qbase', 'Qout','stage']
label=['Stream\nLosses ($m^3/d$)', 'Stream\nBaseflow ($m^3/d$)','Streamflow\n($m^3/d$)', 'Elevation (m)']

var = ['Qaquifer', 'Qout','stage']
label=['Stream\nSeepage ($m^3/d$)','Streamflow\n($m^3/d$)', 'Elevation (m)']
fig, ax = plt.subplots(len(var)+1, 1, sharex=True, figsize=(6.5,8),dpi=300)
pd_sfr.plot(x='Total distance (m)', y='strhc1', ax=ax[0], legend=False)
ax[0].set_ylabel('Stream \n$K_{vert}$ (m/d)')
# pd_sfr.plot(x='Total distance (m)', y='vka', ax=ax[0])
ax[0].set_yscale('log')
ax[-1].plot(pd_sfr['Total distance (m)'], sfr_heads[0], label='GWE')
for n, v in enumerate(var):
    sfrdf.loc[plt_dates].groupby('Total distance (m)').mean(numeric_only=True).plot(y=v, ax=ax[n+1], legend=False)
    # sfrdf.groupby('Total distance (m)').mean(numeric_only=True).reset_index().plot(x='Total distance (m)', y=v, ax=ax[n+1], legend=False)
    ax[n+1].set_ylabel(label[n])
ax[-1].legend()
# sfrdf.loc[plt_date].plot(x='Total distance (m)', y='Qout')

# %%
sfr_sum = sfrdf.groupby(['segment','reach','row','column','layer']).mean(numeric_only=True).reset_index()
# sfr_seg_sum = sfrdf.groupby(['segment']).mean().reset_index()

# %%
sfrdf_sum = sfrdf.resample('D').mean(numeric_only=True)

# plt_date = ['2017-1-1']
var = ['Qout', 'Qrech', 'Qbase']
fig, ax = plt.subplots(len(var), 1, sharex=True, figsize=(6.5,3),dpi=300)
# pd_sfr.plot(x='Total distance (m)', y='strhc1', ax=ax[0])
# pd_sfr.plot(x='Total distance (m)', y='vka', ax=ax[0])
# ax[0].set_yscale('log')

for n, v in enumerate(var):
    # sfrdf.loc[plt_date].plot(x='Total distance (m)', y=v, ax=ax[n+1], legend=False)
    sfrdf_sum.plot( y=v, ax=ax[n], legend=False)
    ax[n].set_ylabel(v)

# %%
# flow is based on old rating curve
mcc_flow = pd.read_csv(join(sfr_dir,'MCC_flow_obs_all.csv'), parse_dates=['DATE TIME'])
mcc_flow = mcc_flow.set_index('DATE TIME')
# mcc_flow.plot(x='DATE TIME',y='flow_cfs')
# I had cleaned the stage data before translating it to flows because it isn't ready to use as is
# mcc_stage = pd.read_csv(join(sfr_dir,'McConnell_stage_2000_01_01_to_2020_12_31.csv'), 
#                         parse_dates=['DATE TIME'], na_values='---')


# %%
sensors = pd.read_csv(join(gwfm_dir, 'Mapping', 'allsensor_latlong.csv'))
sensors = gpd.GeoDataFrame(sensors, geometry=gpd.points_from_xy(sensors.Longitude, sensors.Latitude),
                           crs='epsg:4326').to_crs(grid_sfr.crs)
grid_MCC = gpd.sjoin_nearest(grid_sfr, sensors[sensors.Site_id=='MCC'], how='right')
# subset relevant columns to plot flows at MCC
grid_MCC = grid_MCC[['iseg','ireach','strtop','slope','k','i','j','Elev_m_MSL']].iloc[0]
# grid_MCC.columns

# %%
# sensors

# %%
sfrdf_MCC = sfrdf[(sfrdf.segment==grid_MCC.iseg)&(sfrdf.reach==grid_MCC.ireach)]

mcc_flow_plt = mcc_flow[mcc_flow.index.isin(sfrdf_MCC.index)]
mcc_flow_plt = mcc_flow_plt.reindex(sfrdf_MCC.index)

fig,ax = plt.subplots()
sfrdf_MCC.Qout_cfs.plot(ax=ax, label='Simulated')
mcc_flow_plt.plot(y='flow_cfs', ax=ax, label='Observed')
plt.legend()
ax.set_yscale('log')

# %%
local_sites = pd.read_csv(join(sfr_dir,'flow_obs', 'Stream_monitoring_reference.csv'))

local_sites = gpd.GeoDataFrame(local_sites, geometry=gpd.points_from_xy(local_sites.Longitude, local_sites.Latitude),
                           crs='epsg:4326').to_crs(grid_sfr.crs)
grid_local_sites = gpd.sjoin_nearest(grid_sfr, local_sites, how='right')
# # subset relevant columns to plot flows at MCC
grid_local_sites = grid_local_sites[['iseg','ireach','strtop','slope','k','i','j','site','Site_name']]
# local_sites

# %%
# pull out simulated flow for the stream seg/reach with obs
sfrdf_sites = sfrdf.reset_index().merge(grid_local_sites).set_index('dt').copy()

# %%
# sim_flow.plot(y='Qout_cfs')
# obs_flow.index
# local_flow = local_flow.groupby('site').resample('D').mean(numeric_only=True)

# local_flow

# %%
flow_plt = obs_flow[obs_flow.index.isin(sim_flow.index)]
# flow_plt
# obs_flow.index, sim_flow.index

# %%
local_flow = pd.read_csv(join(sfr_dir,'flow_obs', 'stream_flow_cfs_approx.csv'), 
                         parse_dates=['dt'], index_col='dt')
local_flow = local_flow.groupby('site').resample('D').mean(numeric_only=True).reset_index('site')
local_flow['flow_cmd'] =  local_flow.flow_cfs*86400*0.3048**3

site = 'ACR_189'
def plt_gage_flow(ax, site, local_flow, sfrdf_sites, name=None):
    obs_flow = local_flow[local_flow.site==site]
    sim_flow = sfrdf_sites[sfrdf_sites.site==site]
    
    flow_plt = obs_flow[obs_flow.index.isin(sim_flow.index)]
    flow_plt = flow_plt.reindex(sim_flow.index)
    
    sim_flow.plot(y='Qout_cfs', ax=ax, label='Simulated')
    flow_plt.plot(y='flow_cfs', ax=ax, label='Observed')
    ax.set_xlim(flow_plt.dropna().index.min(), flow_plt.dropna().index.max())
    plt.legend()
    if name is not None:
        ax.set_title(name)
    ax.set_yscale('log')
# rooney_flow.plot(x='dt',y='flow_cfs')
fig,ax = plt.subplots(len(local_sites),1, sharex=True, dpi=300)
for n, s in enumerate(local_sites.site):
    plt_gage_flow(ax[n], s, local_flow, sfrdf_sites,
                  name=s+' - '+local_sites.Site_name.iloc[n])

ax[-1].set_xlabel('Date')
fig.supylabel('Flow (cfs)')


# %% [markdown]
# It makes sense that flows are below the observed because we don't have Deer Creek in the system.
# - even with Deer Creek added it still under predicts
# - 
# We want to look at flows at Rooney and Mahon as well

# %% [markdown]
# # Zone budget check

# %%
def zone_clean(cbc,zon, kstpkper):
    zb = flopy.utils.ZoneBudget(cbc, zon, kstpkper)
    zb_df = zb.get_dataframes()
    # ungroup by timestep
    zb_df = zb_df.reset_index()
    names = zb_df.name.unique()
    zb_df = zb_df.pivot(index = 'totim', columns = 'name',values = 'ZONE_1')
    
    # columns to make negative
    to_cols = zb_df.columns[zb_df.columns.str.contains('TO_')]
    # multiply by -1 to have pulled out of water balance on plot
    zb_df.loc[:, to_cols] *= -1
    # correct for storage change
    # to storage is gw increase (positive)
    stor_cols = zb_df.columns[zb_df.columns.str.contains('STORAGE')]
#     zb_df.loc[:, stor_cols] *= -1
    zb_df['dSTORAGE'] = (zb_df.TO_STORAGE + zb_df.FROM_STORAGE) * -1
    zb_df = zb_df.drop(columns=stor_cols)
    zb_df = zb_df.reset_index()
    strt_date = pd.to_datetime(m.dis.start_datetime)
    zb_df.totim = strt_date+(zb_df.totim*24).astype('timedelta64[h]')
    zb_df = zb_df.set_index('totim')
    # convert 4 hr time steps to daily basis
    zb_df = zb_df.resample('D').mean()
    # summarize to monthly sum
    zb_mon = zb_df.resample('MS').sum()
    zb_mon['PERCENT_ERROR'] = zb_mon['IN-OUT']/np.mean((zb_mon.TOTAL_IN, zb_mon.TOTAL_OUT), axis=0)
    return(zb_df, zb_mon)


# %%
zon_poly = gpd.read_file(join(gwfm_dir, 'Mapping','zonebudget', 'hob_waterbudget_zones.shp'))

# %%
zon_cells = gpd.sjoin(grid_p,zon_poly,how='right',predicate='within')

# filter zone budget 
# could also specify layers
zon_arr = np.zeros((grid_p.row.max(),grid_p.column.max()),dtype=int)
# troubleshooting by assigning different zones to look which drives the most inflow
zon_arr[60:, :] = 3
zon_arr[:60:, :] = 4
zon_arr[:, 175:] = 2
zon_arr[:, :130] = 0

zon_arr[zon_cells.row-1,zon_cells.column-1]=1
plt.imshow(zon_arr)
plt.colorbar(shrink=0.6)


# %%
# could start by just extracting the recharge/pumping inputs and sfr from sfr output
# the cell by cell flows can be check after
zon_rech = rech[:, zon_arr.astype(bool)].sum(axis=1)
zon_wel = pump_rate[:, zon_arr.astype(bool)].sum(axis=1)
# zon_wel = 

# %%
plt.plot(zon_rech, label='Rech')
plt.plot(zon_wel, label='Well')
plt.legend()

# %%
plt_strt = '2019-3-1'
plt_end = '2019-3-31'


# %%
# runs too slowly to properly load all time steps
# running one month is reasonable
all_d, all_mon = zone_clean(cbc, zon_arr, dt_ref.set_index('dt').loc[plt_strt:plt_end].kstpkper)

# %%
zon_cols = all_d.columns[all_d.sum()!=0]
zon_cols = zon_cols[pd.Series(zon_cols).str.contains('TO_|FROM_')]
all_plt = all_d[zon_cols]#.drop(columns = ['IN-OUT','PERCENT_DISCREPANCY','TOTAL_IN','TOTAL_OUT'])
fig,ax = plt.subplots()
all_plt.plot(ax=ax, kind='bar', stacked=True)

ax.plot(all_d['dSTORAGE'].values, color='black')
# all_d.sum()
