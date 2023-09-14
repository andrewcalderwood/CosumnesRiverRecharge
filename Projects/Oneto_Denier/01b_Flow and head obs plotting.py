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

py_dir = doc_dir +'GitHub/CosumnesRiverRecharge/python_utilities/'



# %%
x=2


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy

py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

from mf_utility import get_dates, get_layer_from_elev, clean_wb
from map_cln import gdf_bnds, plt_cln

# %%
ext_dir = 'F:/WRDAPP'
c_dir = 'C:/WRDAPP'
if os.path.exists(ext_dir):
    loadpth = ext_dir 
elif os.path.exists(c_dir):
    loadpth = c_dir 
loadpth +=  '/GWFlowModel/Cosumnes/Stream_seepage'


# model_nam = 'oneto_denier_homogeneous_2014_2018'
upscale = 'upscale4x_'
model_nam = 'oneto_denier_'+upscale+'2014_2018'
# model_nam = 'oneto_denier_'+upscale+'2014_2020'

model_ws = join(loadpth,model_nam)

# model_ws = join(loadpth,'parallel_oneto_denier','realization000')
load_only = ['DIS','UPW','SFR','OC', 'EVT', 'BAS6']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= model_ws, 
                                exe_name='mf-owhm.exe', version='mfnwt',
                              load_only=load_only,
                              )



# %%
print('Quantiles: ',[0,0.5,0.6,0.75,1])
print('HK :',np.quantile(m.upw.hk.array,[0,0.5,0.6,0.75,1]))
print('VKA :',np.quantile(m.upw.vka.array,[0,0.5,0.6,0.75,1]))

# %%
# makes it easier to load if I save a set of outputs with alternate names
m_ver = '' # default no alternate output/input names
# m_ver = '_vka10'


# %%
model_grp = 'inset_oneto_denier'
grid_dir = join(gwfm_dir, 'DIS_data/streambed_seepage/grid')
grid_fn = join(grid_dir, model_grp,'rm_only_grid.shp')
grid_p = gpd.read_file(grid_fn)
grid_p.crs='epsg:32610'
m_domain = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [grid_p.unary_union], crs=grid_p.crs)

# %%
XSg = pd.read_csv(join(model_ws,'04_XSg_filled.csv'))
XSg = gpd.GeoDataFrame(XSg, geometry = gpd.points_from_xy(XSg.Easting, XSg.Northing), crs='epsg:32610')

drop_iseg = XSg[~XSg['Logger Location'].isna()].iseg.values
# overwrite SFR segment/reach input relevant to seepage
# sensor_dict = pd.read_csv(join(model_ws, 'sensor_xs_dict.csv'), index_col=0)
# XS_params = sensor_dict.join(params.set_index('Sensor'), on='Sensor')

# %%
sfr = m.sfr
sfrdf = pd.DataFrame(sfr.reach_data)
grid_sfr = grid_p.set_index(['row','column']).loc[list(zip(sfrdf.i+1,sfrdf.j+1))].reset_index(drop=True)
grid_sfr = pd.concat((grid_sfr,sfrdf),axis=1)

# characterize streambed into different hydrofacies
tprogs_quants = np.array([0.590, 0.155, 0.197, 0.058]).cumsum()
# use facies of vka just below stream
vka_sfr = m.upw.vka.array[sfrdf.k, sfrdf.i, sfrdf.j]
vka_quants = pd.DataFrame(np.quantile(vka_sfr, tprogs_quants))
vka_quants.index=['mud','sandy mud','sand','gravel']
grid_sfr['facies'] = 'mud'
for n in np.arange(0,len(vka_quants)-1):
    grid_sfr.loc[vka_sfr > vka_quants.iloc[n].values[0],'facies'] = vka_quants.index[n+1]
# remove stream segments for routing purposes only
grid_sfr = grid_sfr[~grid_sfr.iseg.isin(drop_iseg)]

# %%
params = pd.read_csv(model_ws+'/ZonePropertiesInitial.csv', index_col='Zone')
# convert from m/s to m/d
params['K_m_d'] = params.K_m_s * 86400 
vka = m.upw.vka.array
tprogs_vals = np.arange(1,5)
tprogs_hist = np.flip([0.590, 0.155, 0.197, 0.058])
tprogs_quants = 1-np.append([0], np.cumsum(tprogs_hist)/np.sum(tprogs_hist))
vka_quants = pd.DataFrame(tprogs_quants[1:], columns=['quant'], index=tprogs_vals)
# dataframe summarizing dominant facies based on quantiles
vka_quants['vka_min'] = np.quantile(vka, tprogs_quants[1:])
vka_quants['vka_max'] = np.quantile(vka, tprogs_quants[:-1])
vka_quants['facies'] = params.loc[tprogs_vals].Lithology.values

# %%
sfrdf = pd.DataFrame(sfr.reach_data)
grid_sfr = grid_p.set_index(['row','column']).loc[list(zip(sfrdf.i+1,sfrdf.j+1))].reset_index(drop=True)
grid_sfr = pd.concat((grid_sfr,sfrdf),axis=1)
# group sfrdf by vka quantiles
sfr_vka = vka[grid_sfr.k, grid_sfr.i, grid_sfr.j]
for p in vka_quants.index:
    facies = vka_quants.loc[p]
    grid_sfr.loc[(sfr_vka< facies.vka_max)&(sfr_vka>= facies.vka_min),'facies'] = facies.facies
#     # add color for facies plots
# grid_sfr = grid_sfr.join(gel_color.set_index('geology')[['color']], on='facies')

# %%
# grid_sfr.plot('facies', legend=True)

# %%
lak_shp = join(gwfm_dir,'LAK_data/floodplain_delineation')
lak_extent = gpd.read_file(join(lak_shp,'LCRFR_ModelDom_2017/LCRFR_2DArea_2015.shp' )).to_crs('epsg:32610')


# %% [markdown]
# ## Sensor data and XS data

# %%
rm_grid = pd.read_csv(join(proj_dir, 'mw_hob_cleaned.csv'))
rm_grid = gpd.GeoDataFrame(rm_grid, geometry = gpd.points_from_xy(rm_grid.Longitude,rm_grid.Latitude), 
                           crs='epsg:4326').to_crs(grid_p.crs)
# get model layer for heads
hob_row = rm_grid.row.values-1
hob_col = rm_grid.column.values-1

# %%
gwl_long = pd.read_csv(join(model_ws,'gwl_long.csv'), parse_dates=['dt'])

# %%
# XS are every 100 m
xs_all = pd.read_csv(dat_dir+'XS_point_elevations.csv',index_col=0)
xs_all = gpd.GeoDataFrame(xs_all,geometry = gpd.points_from_xy(xs_all.Easting,xs_all.Northing), crs='epsg:32610')


# %%

# correspond XS to sensors
rm_elev = gpd.sjoin_nearest(XSg, rm_grid, how='right',lsuffix='xs', rsuffix='rm', distance_col='dist_m')
#MW_11, MW_CP1 had doubles with sjoin_nearest due to XS duplicates from Oneto_Denier
rm_elev = rm_elev.drop_duplicates(['Sensor'])
rm_elev = rm_elev.sort_values('iseg')

# %%
# elevation check for Graham to show the difference between wells and the river
fig,ax=plt.subplots()
# rm_elev
rm_elev.plot(x='Sensor',y='z_m_min', ax=ax, kind='scatter', color='red',label='River Min')
rm_elev.plot(x='Sensor',y='z_m_min_cln',ax=ax, kind='scatter', color='black', label='River Min adj')
rm_elev.plot(x='Sensor',y='MPE (meters)',ax=ax, kind='scatter', color='brown', label='Well MPE')
plt.ylabel('Elevation (m AMSL)')
plt.xticks(rotation=90);
# rm_elev.

# %% [markdown]
# ## Model output - time variant

# %%
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')


# %%
chk_ws = join(loadpth,'parallel_oneto_denier_upscale4x_2014_2018','realization011')
chk_ws


# %%

def clean_wb(flow_name, dt_ref):
    # load summary water budget
    wb = pd.read_csv(flow_name, delimiter=r'\s+')

    wb['kstpkper'] = list(zip(wb.STP-1,wb.PER-1))
    wb = wb.merge(dt_ref, on='kstpkper').set_index('dt')

    # calculate change in storage
    wb['dSTORAGE'] = wb.STORAGE_OUT - wb.STORAGE_IN
    wb['dSTORAGE_sum'] = wb.dSTORAGE.cumsum()
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
wb, wb_out_cols, wb_in_cols = clean_wb(model_ws+'/flow_budget.txt', dt_ref)

# %%
wb[wb_out_cols].mean(), wb[wb_in_cols].mean()

# %%

# %%
# wb_chk, wb_out_cols, wb_in_cols = clean_wb(model_ws+'/flow_budget_old.txt', dt_ref)

# %%
fig,ax= plt.subplots(3,1, sharex=True)
wb.plot(y='PERCENT_ERROR', ax=ax[0])
wb.plot(y=wb_out_cols, ax=ax[1], legend=True)
wb.plot(y=wb_in_cols, ax=ax[2], legend=True)


# %%
# out_var = 'GW_OUT'
# in_var = 'GW_IN'
def wb_chk_plt(var):
    out_var = var+'_OUT'
    in_var = var+'_IN'
    out_chk = (wb[out_var] - wb_chk[out_var])/((wb[out_var]+ wb_chk[out_var])/2)
    in_chk = (wb[in_var]- wb_chk[in_var])/((wb[in_var]+ wb_chk[in_var])/2)
    in_chk.plot(label='in'), out_chk.plot(label='out')
    plt.legend()
# wb_chk_plt('SFR')
# wb_chk_plt('ET')
# wb_chk_plt('WEL')

# %% [markdown]
# The ET water budget does seem to be impacted quite a bit during the wet season because more water is available. SFR is affected on a few random dates to the extreme, but overall not too much. Since there is an impact on ET then it seems like it would make sense. Making the rainfall for GDE areas had the RMSE go from 1.66 to and NSE from 0.66 to .67 so not a major change, but slight improvement.

# %%
# wb_chk.SFR_OUT.plot()
# wb.SFR_OUT.plot(alpha=0.6)


# %% [markdown]
# Adjusting the GHB to use heads at 500 m did reduce the baseflow magnitude and number of days, it would be interesting to see how this impacts NSE. Also the NSe went from 0.60 to 0.66 and RMSE went from 1.8 to 1.6 with the increased outflow which supports the need to re-run all realizations with the udpated GHB.

# %%
et_local = m.evt.evtr.array[:,0]
ext_dp = m.evt.exdp.array[0][0]
ievt = m.evt.ievt.array[0][0]
surf = m.evt.surf.array[0][0]


# %%
# wb.iloc[0]

# %%
# hdobj = flopy.utils.HeadFile(model_ws+'/MF.hds')
# plt.contour(hdobj.get_data((0,0))[-1])
# plt.colorbar()

# %% [markdown]
# # Plot Groundwater Observed vs Simulated
# We need to validate that the local model of Oneto-Denier is adequately representing stream-aquifer interactions so it can be used to quantify storage changes and particle age.
#
# - field data shows fall lows around -2 to 2 m with peaks to 8m with winter
# - simulated data shows fall lows around 3-5 meters with peaks to 8 m
# - the simulated data as usual isn't showing as extreme reactions as the real data is showing which could be partly due to an issue with the initial steady state levels, the peaks match which means the stream stage is peaking close to reality. I may need to consider adding evapotranspiration as the GDEs are likely pulling water
# - adding EVT and RCH together brought the steady state conditions to match observed conditions, but way over-estimated the winter peaks

# %%
# hdobj = flopy.utils.HeadFile(model_ws+'/MF.hds')
# # extract time series of heads for each desired location
# mw_hds = hdobj.get_ts(list(zip(rm_grid['lay'], hob_row, hob_col)))
# mw_hds = pd.DataFrame(mw_hds, columns=['time']+rm_grid.Sensor.tolist())
# # convert to hourly to maintain more precision in DT
# mw_hds['dt'] = strt_date+(mw_hds.time.values*24 ).astype('timedelta64[h]')
# mw_gwl = mw_hds.drop(columns=['time'])
# # long format for id join with observed dat
# mw_long = mw_gwl.melt(id_vars='dt', var_name='Well',value_name='sim')
# mw_long = mw_long[mw_long.sim != -1e30]

# %%
# mw_chk = mw_long.join(gwl_long.set_index(['Well','dt']), on=['Well','dt'], how='inner')
# mw_chk = mw_chk.melt(id_vars=['dt', 'Well'],value_vars=['sim','obs'], value_name='gwe', var_name='type')


# %%
# sns.relplot(mw_chk,x='dt',y='gwe',col='Well', hue='type', col_wrap=4)

# %%
from sklearn.metrics import r2_score, mean_squared_error

def nse(targets,predictions):
    return 1-(np.sum((targets-predictions)**2)/np.sum((targets-np.mean(predictions))**2))

def clean_hob(model_ws, name='MF.hob.out'):
    hobout = pd.read_csv(join(model_ws,name),delimiter=r'\s+', header = 0,names = ['sim_val','obs_val','obs_nam'],
                         dtype = {'sim_val':float,'obs_val':float,'obs_nam':object})
    hobout[['Sensor', 'spd']] = hobout.obs_nam.str.split('p',n=2, expand=True)
    hobout['kstpkper'] = list(zip(np.full(len(hobout),0), hobout.spd.astype(int)))
    hobout = hobout.join(dt_ref.set_index('kstpkper'), on='kstpkper')
    hobout.loc[hobout.sim_val.isin([-1e30, -999.99, -9999]), 'sim_val'] = np.nan
    hobout = hobout.dropna(subset='sim_val')
    hobout['error'] = hobout.obs_val - hobout.sim_val
    hobout['sq_error'] = hobout.error**2
    
    return(hobout)



# %%
# hobout = clean_hob(chk_ws)
hobout = clean_hob(model_ws, 'MF.hob_chk.out')

# removing oneto ag because of large depth offset
hobout = hobout[hobout.Sensor != 'MW_OA']

# %%
# technically the RMSE is MSE unless I square root it

# summary stats by well
mw_stats = hobout[['Sensor','sq_error']].groupby('Sensor').sum()
mw_stats['r2'] = 0
for s in hobout.Sensor.unique():
    df_s = hobout[hobout.Sensor==s]
    mw_stats.loc[s,'r2'] = r2_score(df_s.obs_val, df_s.sim_val)
    mw_stats.loc[s,'RMSE'] = mean_squared_error(df_s.obs_val, df_s.sim_val, squared=False) # false returns RMSE instead of MSE
    mw_stats.loc[s,'NSE'] = nse(df_s.obs_val, df_s.sim_val)

t=0
sum_stats = pd.DataFrame(columns=['r2','RMSE','NSE'])
# summary statistics
sum_stats.loc[t,'r2'] = r2_score(hobout.obs_val, hobout.sim_val)
sum_stats.loc[t,'RMSE'] = np.sqrt(hobout.sq_error.sum()/len(hobout))
sum_stats.loc[t,'NSE'] = nse(hobout.obs_val, hobout.sim_val)

sum_stats

# %%
hob_long = hobout.melt(id_vars=['dt', 'Sensor'],value_vars=['sim_val','obs_val'], value_name='gwe', var_name='type')
# hob_long
# hob_long = hobout.melt(id_vars=['dt', 'Sensor'],value_vars=['sim_val','obs_val','sim_4x'], value_name='gwe', var_name='type')


# %%
# hob_long, x='dt',y='
g = sns.relplot(hob_long, x='dt',y='gwe',col='Sensor',hue = 'type',  col_wrap=4)

axes = g.axes.flatten()
mw = hob_long.Sensor.unique()

for n in np.arange(0,len(axes)):
    mw_dat = rm_elev[rm_elev.Sensor ==mw[n]]
    axes[n].axhline(mw_dat['MPE (meters)'].values[0], ls='--', linewidth=3, color='brown')
    axes[n].axhline(mw_dat['z_m_min_cln'].values[0]-1, ls='--', linewidth=3, color='blue')
    et_bot = (surf-ext_dp)[mw_dat.row_rm.iloc[0], mw_dat.column_rm.iloc[0]]
    axes[n].axhline(et_bot, ls='--', linewidth=3, color='green')
#     axes[n].axhline(mw_dat['bot_screen_m'].values[0]-1, ls='--', linewidth=3, color='black')

# %%
# aggregate error for spatial plotting by month average?

hob_diff = hob_long.pivot_table(index=['dt','Sensor'],values='gwe',columns='type')
hob_diff['h_diff'] = hob_diff.sim_val - hob_diff.obs_val

hob_diff_mon = hob_diff.reset_index().set_index('dt').groupby('Sensor').resample('MS').mean()
hob_diff_mon = hob_diff_mon[['h_diff']].reset_index()
# hob_diff_mon

# %%
# t_plt = '2017-11-01'
# t_plt = '2018-06-01'
t_plt = '2018-02-01'

diff_plt = hob_diff_mon[hob_diff_mon.dt == t_plt]
# diff_plt['sign'] = np.sign(diff_plt.h_diff)
diff_plt = rm_grid.join(diff_plt.set_index('Sensor'),on='Sensor')

fig,ax=plt.subplots(figsize=(8, 8))
# m_domain.plot(ax=ax,color='None')
mapview = flopy.plot.PlotMapView(model=m,ax=ax)

gdf_bnd = gdf_bnds(rm_grid,buf=400,ax=ax)

spd = dt_ref[dt_ref.dt==t_plt].kstpkper.values[0]
avg_sfr_lay = int(np.round(grid_sfr.k.mean()))
head = hdobj.get_data((0,spd[-1]))[avg_sfr_lay] #m.dis.top.array - 
# head = np.ma.masked_where(head==-1e30, head)
head[head==-1e30] = np.nan
# m_domain.plot(ax=ax_n,color='none')
im = mapview.contour_array(head, masked_values=[-999.99, -1e30], ax=ax)
# gdf_bnd = gdf_bnds(rm_grid,buf=400,ax=ax)

plt.colorbar(im, ax=ax, shrink = 0.4)
plt.clabel(im)

grid_sfr.plot(color='blue', ax=ax)
lak_extent.plot(color='none', ax=ax)
diff_plt.plot('h_diff', scheme='Quantiles', k = 6, ax=ax,
                  legend=True,cmap='bwr', legend_kwds={'loc':'lower right' ,'title':'Error (Sim - Obs)'})

rm_grid.apply(lambda x: ax.annotate(x.Sensor.replace('MW_',''), xy=x.geometry.coords[0], ha='center', fontsize=6,
                                    xytext = (5,10), textcoords='offset pixels',
#                                     arrowprops = {'shrink':1},
                                    bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2)
                                                        ),axis=1);
ax.set_title(t_plt)

ctx.add_basemap(ax=ax, source = ctx.providers.Esri.WorldImagery, attribution=False, attribution_size=6,
                crs = 'epsg:26910', alpha=0.6)

# %%
mcc_d = pd.read_csv(sfr_dir+'MCC_flow_obs_all.csv', parse_dates = ['DATE TIME'], index_col='DATE TIME')
mcc_d = mcc_d[(mcc_d.index>strt_date)&(mcc_d.index<end_date)]



# %%
# load local stream stage data

# %% [markdown]
# ## Lake plotting

# %%
# elevSteps, volArray, saArray
bathtxt = np.loadtxt(m.model_ws+'/MF.bath', delimiter = '\t')
bath = pd.DataFrame(bathtxt, columns=['elev','vol','area'])
# bath.plot(x='elev',y='vol')

# %%
gage_cols = ['time','stage','volume','conc','inflows','outflows','conductance','error']

def read_gage(gagenam):
    gage = pd.read_csv(gagenam,skiprows=1, delimiter = r'\s+', engine='python')
    cols = gage.columns[1:-1]
    gage = gage.dropna(axis=1)
    gage.columns = cols
    strt_date = pd.to_datetime(m.dis.start_datetime)
    gage['dt'] = strt_date+(gage.Time*24).astype('timedelta64[h]')
    gage = gage.set_index('dt')
    gage['dVolume'] = gage.Volume.diff()
    gage['Total_In'] = gage[['Precip.','Runoff','GW-Inflw','SW-Inflw']].sum(axis=1)
    gage['Total_Out'] = gage[['Evap.','Withdrawal','GW-Outflw','SW-Outflw']].sum(axis=1)
    gage['In-Out'] = gage.Total_In - gage.Total_Out
#     gage['name'] = run
    return(gage)



# %%
lak_out = read_gage(join(model_ws, 'MF_lak.go'))

# %%
# lak_out[lak_out['Percent-Err']>10]

# %%
# model troubleshooting
# fig,ax = plt.subplots(4,1, sharex=True)
# lak_out.plot(y=['Total_In','Total_Out'], ax=ax[0])
# # plt.yscale('log')
# lak_out.plot(y=['In-Out','dVolume'],ax=ax[1])
# # (lak_out['In-Out']-lak_out.dVolume
# lak_out.plot(y='Volume',ax=ax[2])
# lak_out.plot(y=['Percent-Err'],ax=ax[-1])


# %%
lak_out['2015-1-1':'2015-10-1']['Stage(H)'].min()
# min lake stage for lake out is 2.95
# outflow is zero in summer, stage is never 0, volume is zero in summer
# gw inflow is zero in 2015, 2016 summers
# np.sign(lak_out['SW-Outflw']).plot()

# %%
fig,ax = plt.subplots(3,1, sharex=True)
lak_out.plot(y=['Stage(H)'], ax = ax[0]) #dry all the time
# lak_out.columns # GW-Inflw, GW-Outflw, SW-Inflw, SW-Outflw
lak_out.plot(y=['GW-Inflw', 'GW-Outflw'], ax=ax[2]) # there is gw inflow
lak_out.plot(y=['SW-Inflw', 'SW-Outflw'], ax=ax[1]) # there is sw inflow


# %%
fig,ax = plt.subplots(5,1, sharex=True, layout='constrained')
for n, wb_n in enumerate(wb_out_cols):
    wb.plot(y=wb_n, ax=ax[n], legend=False)
    ax[n].set_ylabel(wb_out_cols[n].split('_')[0])

# %% [markdown]
# ## SFR Plotting

# %%
# grid_sfr = pd.DataFrame().from_records(m.sfr.reach_data).rename(columns={'i':'row','j':'column'})
# grid_sfr[['row','column']] += 1 # convert to 1 based to match with SFR output
pd_sfr = grid_sfr.set_index(['iseg','ireach'])[['rchlen','strtop', 'facies']]
pd_sfr['Total distance (m)'] = pd_sfr['rchlen'].cumsum()

def clean_sfr_df(model_ws):
    sfrout = flopy.utils.SfrFile(join(model_ws, m.name+m_ver+'.sfr.out'))
    sfrdf = sfrout.get_dataframe()
    sfrdf = sfrdf.join(dt_ref.set_index('kstpkper'), on='kstpkper').set_index('dt')
    # convert from sub-daily to daily using mean, lose kstpkper
    sfrdf = sfrdf.groupby('segment').resample('D').mean(numeric_only=True)
    sfrdf = sfrdf.reset_index('segment', drop=True)
    sfrdf[['row','column']]-=1 # convert to python
    cmd2cfs = 1/((0.3048**3)*86400) # cubic meters per day to cfs
    sfrdf['month'] = sfrdf.index.month
    sfrdf['WY'] = sfrdf.index.year
    sfrdf.loc[sfrdf.month>=10, 'WY'] +=1
    # add column to track days with flow
    sfrdf['flowing'] = 1
    sfrdf.loc[sfrdf.Qout <= 0, 'flowing'] = 0
#     sfrdf = pd_sfr.join(sfrdf.set_index(['row','column']),on=['row','column'],how='inner',lsuffix='_all')
    sfrdf = sfrdf.join(pd_sfr ,on=['segment','reach'],how='inner',lsuffix='_all')
    # dependent on number of time steps
    sfrdf['Qin_cfs'] = sfrdf.Qin * cmd2cfs
    sfrdf['Qout_cfs'] = sfrdf.Qout * cmd2cfs
    sfrdf['Qaquifer_cfs'] = sfrdf.Qaquifer * cmd2cfs
    
    # create different column for stream losing vs gaining seeapge
    sfrdf['Qrech'] = np.where(sfrdf.Qaquifer>0, sfrdf.Qaquifer,0)
    sfrdf['Qbase'] = np.where(sfrdf.Qaquifer<0, sfrdf.Qaquifer*-1,0 )
    # booleans for plotting
    sfrdf['gaining'] = (sfrdf.gradient == 0)
    sfrdf['losing'] = (sfrdf.gradient >= 0)
    sfrdf['connected'] = (sfrdf.gradient < 1)
    return(sfrdf)



# %%
sfrdf =  clean_sfr_df(model_ws)

# %%
# find last day of flow

# start simple with just year by segment ,'month','facies'
sns.relplot(sfrdf.groupby(['WY','segment']).sum(numeric_only=True), x='segment',y='flowing', hue='WY')

# %%
# fig, ax = plt.subplots(1,3, figsize=(16,4))
# for s in [1,30,55]:
#     sfrdf[sfrdf.segment==s].plot(y='Qout_cfs',ax=ax[0], label=s)
#     sfrdf[sfrdf.segment==s].plot(y='depth',ax=ax[1], label=s)
#     sfrdf[sfrdf.segment==s].plot(y='Qaquifer_cfs',ax=ax[2], label=s)


# %%
# sfrdf.groupby('segment').mean()['gradient'].plot()

# %%

x = int(grid_sfr.facies.isin(['gravel','sand']).sum())
y = sfrdf[sfrdf.facies.isin(['gravel','sand'])].resample('AS').mean(numeric_only=True)
y_plt = y[['Qbase','connected']]
y_plt['num_coarse'] = x
y_plt['year'] = y_plt.index.year
y_sym = ['.','x','*','s', 'D']

fig, ax = plt.subplots(1,2,figsize=(8,3))
ax[0].set_ylabel('Mean Annual Baseflow ($m^3/day$)')
ax[0].set_xlabel('Number of Coarse Stream Segments')    
ax[1].set_ylabel('Fraction of Time Connected')
ax[1].set_xlabel('Number of Coarse Stream Segments')
    
for n, yr in enumerate(y_plt.year):
    # ax[0].scatter(np.repeat(x, len(y1)), y1.Qbase.values)
    y_plt[y_plt.year==yr].plot.scatter('num_coarse','Qbase', ax=ax[0],  marker=y_sym[n])

    # ax[1].scatter(np.repeat(x, len(y1)), y1.connected.values)
    y_plt[y_plt.year==yr].plot.scatter('num_coarse','connected', ax=ax[1], marker=y_sym[n])


plt.legend(y_plt.year)

# %%
# grid_sfr[['iseg','ireach','facies']]
sfr_facies_sum = sfrdf.groupby(['dt','facies']).sum(numeric_only=True)
seep_facies_sum = sfr_facies_sum[['Qrech','Qbase']].melt(ignore_index=False)



# %% [markdown]
# The gaining/losing work I did with Stephen Maples shows periods of connection and disconnection and if we assume that the magnitude of stream stage is higher in the wet years then likely the groundwater system is losing in those years as well. What is likely true about floodplains is that the system is predominantly losing except during flood periods when there are more complex local scale gaining/losing conditions. 
#
# What this work shows is that gravels/sands are more active during these extreme wet periods.

# %%
# som eissue with sharex is hiding mud probably issue of dt type, difference between mud and sandy mud
fig,ax = plt.subplots(2,2, figsize=(12,8), sharex=True, sharey=True)#

df_rech= seep_facies_sum[seep_facies_sum.variable=='Qrech'].reset_index('facies')
for n, f in enumerate(['mud','sandy mud','sand','gravel']):
    ax_n = ax[int(n/2), n%2]
    df_plt = df_rech[df_rech.facies==f]
    df_plt.index = pd.to_datetime(df_plt.index)
    df_plt.plot(y='value', ax=ax_n, legend=False)
    ax_n.set_title(f)
    ax_n.set_yscale('log')

# %%
# sfr_facies_sum[sfr_facies_sum['Qbase']>0]
# sfr_facies_sum


# %%
g = sns.relplot(sfr_facies_sum.melt(value_vars=['Qbase','Qrech'], ignore_index=False),
                x='dt',y='value',hue='variable', kind='line',
            col = 'facies', col_wrap=2)
g.set(yscale='log')

# %% [markdown]
# ### Plot stream discretization

# %%

# spd_hd = dt_ref[dt_ref.dt == '2020-05-21'].kstpkper.values[0]
# head = hdobj.get_data(spd_hd)[0][0]

for t in spd_stp[0::90]: # every 7 days 
#     spd_hd = dt_ref[dt_ref.dt == t].kstpkper.values[0]
    head = hdobj.get_data(t)[grid_sfr.k, grid_sfr.i, grid_sfr.j]
    head = head[head!=-1e30]
    plt.plot(head, color='lightgray')
plt.plot(head,label = 'GWE',  color='lightgray')
plt.plot(m.dis.top.array[grid_sfr.i, grid_sfr.j], label='Model Top', ls='--',color='black')
plt.plot(m.sfr.reach_data.strtop, label= 'Stream Top', ls=':',color='black')
plt.plot(m.sfr.reach_data.strtop-m.sfr.reach_data.strthick, label= 'Stream Bottom', ls=':',color='black')

plt.legend()

# %% [markdown]
# ## Grid wide head distribution

# %%
avg_sfr_lay = int(np.round(grid_sfr.k.mean()))

nx = 2
ny = 4
fig,ax = plt.subplots(ny,nx, figsize=(12,12),sharex=True, sharey=True)

# fig.tight_layout()
for n,t in enumerate(spd_stp[::180][1:]):
    head = hdobj.get_data(t)[avg_sfr_lay] #m.dis.top.array - 
        
    ax_n = ax[int(n / nx), n % nx]
    mapview = flopy.plot.PlotMapView(model=m,ax=ax_n)
    m_domain.plot(ax=ax_n,color='none')
    im = mapview.contour_array(head, masked_values=[-999.99], ax=ax_n)
    grid_sfr.plot(ax=ax_n)
# vmin, vmax from visual inspection but could be added with code
#     im = ax_n.contour(head[avg_sfr_lay])

    ax_n.set_aspect(1)
    plt.colorbar(im, ax=ax_n, shrink = 0.4)
# fig.subplots_adjust(wspace=0.2, hspace=-.5)
fig.tight_layout()

# %%
