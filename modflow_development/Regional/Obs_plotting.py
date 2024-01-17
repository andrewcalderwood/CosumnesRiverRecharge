# ---
# jupyter:
#   jupytext:
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
import time

import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio

from sklearn.metrics import r2_score, mean_squared_error

# standard python plotting utilities
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt


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
run_dir = 'C://WRDAPP/GWFlowModel'
=run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'
base_model_ws = join(loadpth, model_nam)
# model_nam = 'foothill_vani10'
model_nam = 'strhc1_scale'
model_nam = 'parallel_realizations/realization005'
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
# gel = flopy.modflow.ModflowUpw.load(join(model_ws, 'MF.upw'), model=m)

# %%
sfr_dir = gwfm_dir+'/SFR_data/'
# grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')
# grid_sfr['Kz'] = m.sfr.reach_data.strhc1
m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_p['easting'] = grid_p.geometry.centroid.x
grid_p['northing'] = grid_p.geometry.centroid.y

lak_grid_clip = gpd.read_file(gwfm_dir+'/Levee_setback/lak_grid_clip/lak_grid_clip.shp')


# %%

vka = gel.vka.array
# load sfr data 
grid_sfr = reach_data_gdf(sfr, grid_p)
grid_sfr[['row','column']] = grid_sfr[['i','j']] +1 # convert to 1 based to match with SFR output
drop_iseg = grid_sfr[grid_sfr.strhc1==0].iseg.values
grid_sfr['vka'] = vka[grid_sfr.k, grid_sfr.i, grid_sfr.j]
vka_quants = pd.read_csv(join(model_ws, 'vka_quants.csv'))

for p in vka_quants.index:
    facies = vka_quants.loc[p]
    grid_sfr.loc[(sfr_vka< facies.vka_max)&(sfr_vka>= facies.vka_min),'facies'] = facies.facies
    # add color for facies plots

# drop routing segments before calculating distances
pd_sfr = grid_sfr[~grid_sfr.iseg.isin(drop_iseg)]
pd_sfr = pd_sfr.set_index(['iseg','ireach'])[['rchlen','strtop','strhc1', 'vka', 'facies']]
pd_sfr['Total distance (m)'] = pd_sfr['rchlen'].cumsum()


# %%
all_obs = pd.read_csv(model_ws+'/input_data/all_obs_grid_prepared.csv',index_col=0)
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
# plt.imshow(rech.mean(axis=0), vmax=0.001)
# plt.imshow(pump_rate.mean(axis=0), vmax=0.005)

# plt.colorbar(shrink=0.6)

# %% [markdown]
# ## Water Budget check

# %%
def clean_wb(model_ws, dt_ref):
    # load summary water budget
    wb = pd.read_csv(model_ws+'/flow_budget.txt', delimiter=r'\s+')

    wb['kstpkper'] = list(zip(wb.STP-1,wb.PER-1))
    wb = wb.merge(dt_ref, on='kstpkper').set_index('dt')

    # calculate change in storage
    wb['dSTORAGE'] = wb.STORAGE_OUT - wb.STORAGE_IN
    # calculate total gw flow, sum GHB, CHD
    wb['GW_OUT'] = wb.GHB_OUT + wb.CHD_OUT
    wb['GW_IN'] = wb.GHB_IN + wb.CHD_IN
    wb = wb.loc[:,~wb.columns.str.contains('GHB|CHD')]
    
    wb_cols = wb.columns[wb.columns.str.contains('_IN|_OUT')]
    wb_cols = wb_cols[~wb_cols.str.contains('STORAGE|IN_OUT')]
    wb_out_cols= wb_cols[wb_cols.str.contains('_OUT')]
    wb_in_cols = wb_cols[wb_cols.str.contains('_IN')]
    # only include columns with values used
    wb_out_cols = wb_out_cols[np.sum(wb[wb_out_cols]>0, axis=0).astype(bool)]
    wb_in_cols = wb_in_cols[np.sum(wb[wb_in_cols]>0, axis=0).astype(bool)]

    return(wb, wb_out_cols, wb_in_cols)


# %%
wb, out_cols, in_cols = clean_wb(model_ws, dt_ref)

# %%
wb_ann = wb.resample('AS-Oct').sum(numeric_only=True)
fig,ax = plt.subplots( sharex=True)
plt.axhline(0, color='black')
wb_ann[out_cols].multiply(-1).plot( kind='bar', ax=ax, stacked=True)
wb_ann[in_cols].plot( kind='bar', ax=ax, stacked=True)
plt.legend(loc=(1.05,0.3))

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

# # general contour check

# %%
ghb_dir = join(gwfm_dir, 'GHB_data')
year = 2019 # 2016
filename = glob.glob(ghb_dir+'/final_WSEL_arrays/fall'+str(year)+'_kriged_WSEL.tsv')[0]
# convert from ft to meters
hd_strt = np.loadtxt(filename)*0.3048
extent = (minx,maxx,miny,maxy)

# %%
# minx, miny, maxx, maxy = m_domain.geometry.unary_union.bounds
# extent = (minx, maxx, miny, maxy)

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
y=2019
s='fall'
row = 50
hob_gpd_plt = hob_gpd[(hob_gpd.index.year==y)&(hob_gpd.season==s)]
head = hdobj.get_data((0,int(hob_gpd_plt.spd.mean())))
head_ma = np.ma.masked_where(head==-999.99, head)  
hmin, hmax = base_round(head_ma.min(), step), base_round(head_ma.max(), step)
levels = np.arange(hmin, hmax, step)
plt.plot(head_ma[0][row,:],label='Simulation')
# plt.plot(head_ma[-2][row,:],label='Simulation Laguna')
# plt.plot(head_ma[-1][row,:],label='Simulation Mehrten')
plt.plot(hd_strt[row,:], label='Contour')
plt.legend()
print('simulated min %.2f' %head_ma[0][row,:].min(),'and observed %.2f' %hd_strt[row,:].min())


# %% [markdown]
# # Sim vs Obs Head
#

# %%
# model_ws = loadpth+'historical_simple_geology'
# model_ws = loadpth+'historical_simple_geology_reconnection'


# %%
def load_hob(model_ws):
    hobout = pd.read_csv(model_ws+'/MF.hob.out',delimiter=r'\s+', header = 0,names = ['sim_val','obs_val','obs_nam'],
                         dtype = {'sim_val':float,'obs_val':float,'obs_nam':object},
                        na_values=[-9999.])
    # if only one obs exists correct naming convention
    one_obs = ~hobout.obs_nam.str.contains('.0')
    hobout.loc[one_obs,'obs_nam'] = hobout.loc[one_obs,'obs_nam']+'.'+str(1).zfill(5)
    return(hobout)
    
hobout = load_hob(model_ws)

# %%
fig, ax = plt.subplots(1,1,figsize=(5,5))

# get boundary values for plotting a 1:1
hobmax = hobout.loc[:,['sim_val','obs_val']].max().min()
hobmin = hobout.loc[:,['sim_val','obs_val']].min().max()

# plot observed vs simulated values
hobout.plot.scatter(x='obs_val', y='sim_val',marker='.',ax=ax)
ax.plot([hobmin,hobmax],[hobmin,hobmax],'red')
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
def mak_hob_gpd(hobout, model_ws):

    all_obs = pd.read_csv(model_ws+'/input_data/all_obs_grid_prepared.csv',index_col=0)
    all_obs.index = all_obs.index.rename('date')
    all_obs = all_obs.reset_index()
    # join more indepth obs data to output simulated heads
    obs_data = hobout.join(all_obs.set_index('obs_nam'),on=['obs_nam'], how='right')
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
    hob_gpd['error'] = hob_gpd.obs_val - hob_gpd.sim_val
    hob_gpd['abs_error'] = hob_gpd.error.abs()
    
    if 'date' in hob_gpd.columns:
        hob_gpd = hob_gpd.set_index('date')
        hob_gpd.index = pd.to_datetime(hob_gpd.index)    
        #     groupby values by season
        hob_gpd.loc[(hob_gpd.index.month > 2)&(hob_gpd.index.month < 6),'season'] = 'spring'
        hob_gpd.loc[(hob_gpd.index.month > 8)&(hob_gpd.index.month < 12),'season'] = 'fall'
    
    return(hob_gpd)
    # set date
    


# %%
hob_gpd = mak_hob_gpd(hobout, base_model_ws)


hob_seasonal = hob_gpd.groupby(['node','season']).mean(numeric_only=True)
hob_seasonal = gpd.GeoDataFrame(hob_seasonal, geometry = gpd.points_from_xy(hob_seasonal.easting, hob_seasonal.northing))
hob_seasonal = hob_seasonal.reset_index()
    


# %%
stns = hob_gpd.drop_duplicates('site_code', keep='last').reset_index().drop(columns=['date','gwe'])
stns['botm_elev'] = m.dis.botm[stns.layer-1, stns.row-1, stns.column-1]
stns.crs = hob_gpd.crs

# %%
hob_gpd['Statistic'] = 0.01
hob_gpd['StatFlag'] = 'SD'
# locations with significant difference between RPE GSE and the DEM should have additional uncertainty included
hob_gpd['Statistic'] += np.round(np.abs(hob_gpd.dem_wlm_gse),4)
hob_gpd['Weight'] = 1/(hob_gpd.Statistic**2)

soswr = (np.sum(np.abs(hob_gpd.sim_val-hob_gpd.obs_val)*hob_gpd.Weight))
print('Sum of absolute difference of OBS and SIM: %.2e' %soswr)


# %%
from report_cln import base_round
from map_cln import plt_cln


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
                          legend=True,cmap='bwr',legend_kwds={'loc':(1.1,0.9),'title':'Error (Obs - Sim)'})
        hob_gpd_plt.apply(lambda x: ax.annotate(str(x.node), xy=(x.geometry.x, x.geometry.y), ha='right'),axis=1);
    
        gdf_bnds(hob_gpd_plt, ax=ax, buf=2E3)
        return(hob_gpd_plt)
    plt_cln(ax=ax)
    return None

    # ax.legend(loc=(1,0.5))


# %%
# plot plain contours for reference
# hob_gpd_plt = plt_hob_map(2019, 'fall', hob=False, rch=False, contour=True, hk=False, step=5)
hob_gpd_plt = plt_hob_map(2019, 'fall', hob=True, rch=False, contour=False, hk=False, step=5)

# %%
# nd_chk = [15343, 16733, 11448, 8437, 15314, 14626] +[3103, 5642, 6112, 10746, 6458]
# # nd_chk = [2926, 8437, 12944, 13407]
# nd_chk = [15314, 15343, 13407, 12944, 14626]
# nd_chk = [6458, 8437, 9580, 11448, 15314, 20055]
# # nd_chk = [6458, 8437, 9580, 10884, 11448]
# hob_gpd_chk = plt_hob_map(2019, 'fall', nd_chk=nd_chk, rch=True, contour=False)


# %%
# df_chk = hob_gpd_chk.groupby('node').mean(numeric_only=True)
# df_chk['rech'] = rech.sum(axis=0)[df_chk.row.astype(int)-1, df_chk.column.astype(int)-1]
# df_chk[['sim_val','obs_val','avg_screen_depth', 'hk', 'rech', 'layer','abs_error']]

# %%
hobout = load_hob(model_ws)

hob_gpd = mak_hob_gpd(hobout)
# find sites with long time series of OBS
hobs_long = (hob_gpd.groupby('site_code').count()>=int(m.dis.nper/365)*2)
hobs_long = hobs_long.index[hobs_long.obs_val].values

# hob_gpd.site_code.unique().shape, voi.shape
# hob_long = hob_gpd[hob_gpd.site_code.isin(hobs_long)]
# hob_long = hob_ts_chk.melt(value_vars=['sim_val','obs_val', 'sim_new'], id_vars=['node'], ignore_index=False)
hob_long = hob_gpd.melt(value_vars=['sim_val','obs_val'], id_vars=['node'], ignore_index=False)

# # a few wells have duplicates in a node (same with site_code), temp fix
# issue was actually the NA values
hob_long = hob_long.reset_index().drop_duplicates(['date','node','variable']).set_index('date')
# # hob_long

# %%
# chk_lng = hob_long.groupby('node').count().variable
# chk_lng = chk_lng[chk_lng>50].index
# sns.relplot(hob_long[hob_long.node.isin(chk_lng)], x='date',y='value', 
sns.relplot(hob_long.dropna(subset='value'), x='date',y='value', 
            hue='variable',  col='node',
#             col_wrap=10, # for powerpoint
            col_wrap=4,
           facet_kws={'sharex':True, 'sharey':False}
           )


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
hdobj = flopy.utils.HeadFile(model_ws+'/MF.hds')
fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=300) # was 40, 20 before

plt.aspect=10

head = hdobj.get_data(kstpkper = spd_stp[0])
head_new = hdobj.get_data(kstpkper = spd_stp[1400])

rownum = 50
# rownum = 0

mcs = flopy.plot.PlotCrossSection(model=m, line={'Row' : rownum})
# colnum = 150
# mcs = flopy.plot.PlotCrossSection(model=m, line={'Column' : colnum})

linecollection = mcs.plot_grid(linewidth = 0.3)
ax.add_collection(linecollection)

mcs.plot_array(a=gel.hk.array, norm = mpl.colors.LogNorm())

wt = mcs.plot_surface(a=hd_strt[:,:], color='black')
wt = mcs.plot_surface(a=head[:,:,:], color='blue')
wt = mcs.plot_surface(a=head_new[:,:,:],color='red')

plt.xlabel('Distance from southwestern edge (m)')
plt.ylabel('Elevation (m)')


# %%
head_new[-1,50,173], head_new[-1,50,200], head_new[-1,50,229]

# %%
head_new[-1,50,175::5]

# %%
hob_gpd_chk[hob_gpd_chk.node==16733].column

# %% [markdown]
# # Plot stream water budget

# %%

from mf_utility import clean_sfr_df
from flopy_utilities import sfr_load_hds

# %%
sfrdf = clean_sfr_df(model_ws, dt_ref, pd_sfr)
# drop routing segments
sfrdf = sfrdf[~sfrdf.segment.isin(drop_iseg)]

# %%
# plt_dates = ['2016-1-1']
plt_dates = pd.date_range('2017-1-1','2017-4-1')
sfr_heads, avg_heads = sfr_load_hds(hdobj, grid_sfr[~grid_sfr.iseg.isin(drop_iseg)], plt_dates, dt_ref)

# %%

var = ['Qrech','Qbase', 'Qout','stage']
label=['Stream\nLosses ($m^3/d$)', 'Stream\nBaseflow ($m^3/d$)','Streamflow\n($m^3/d$)', 'Elevation (m)']
fig, ax = plt.subplots(len(var)+1, 1, sharex=True, figsize=(6.5,8),dpi=300)
pd_sfr.plot(x='Total distance (m)', y='strhc1', ax=ax[0], legend=False)
ax[0].set_ylabel('Stream \n$K_{vert}$ (m/d)')
# pd_sfr.plot(x='Total distance (m)', y='vka', ax=ax[0])
ax[0].set_yscale('log')
ax[4].plot(pd_sfr['Total distance (m)'], sfr_heads[0], label='GWE')
for n, v in enumerate(var):
    # sfrdf.loc[plt_dates].groupby('Total distance (m)').mean(numeric_only=True).plot(y=v, ax=ax[n+1], legend=False)
    sfrdf.groupby('Total distance (m)').mean(numeric_only=True).reset_index().plot(x='Total distance (m)', y=v, ax=ax[n+1], legend=False)
    ax[n+1].set_ylabel(label[n])
ax[4].legend()
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


# %% [markdown]
#

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
mcc_flow_plt = mcc_flow[mcc_flow.index.isin(sfrdf_MCC.index)]
mcc_flow_plt = mcc_flow_plt.reindex(sfrdf_MCC.index)

# %%
sfrdf_MCC = sfrdf[(sfrdf.segment==grid_MCC.iseg)&(sfrdf.reach==grid_MCC.ireach)]
fig,ax = plt.subplots()
sfrdf_MCC.Qout.plot(ax=ax, label='Simulated')
mcc_flow_plt.plot(y='flow_cmd', ax=ax, label='Observed')
plt.legend()
ax.set_yscale('log')


# %% [markdown]
# It makes sense that flows are below the observed because we don't have Deer Creek in the system.

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

# filter zone budget for Blodgett Dam to just within 5 cells or so of the Dam
zon_arr = np.zeros((grid_p.row.max(),grid_p.column.max()),dtype=int)
zon_arr[zon_cells.row-1,zon_cells.column-1]=1


# %%
all_d, all_mon = zone_clean(cbc, zon_arr, dt_ref.kstpkper)

