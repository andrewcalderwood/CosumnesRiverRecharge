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
import sys
from os.path import basename, dirname, join, exists, expanduser
import glob
import pandas as pd
import numpy as np
import time
from scipy.stats import gmean

# standard geospatial python utilities
import geopandas as gpd
from osgeo import gdal
import rasterio

# import flopy
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

# %%
# statistics functions
from sklearn.metrics import r2_score, mean_squared_error
from sklearn import datasets, linear_model
from scipy.stats import pearsonr, spearmanr, kendalltau


# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study
proj_dir = gwfm_dir + '/Oneto_Denier/'
dat_dir = proj_dir+'Stream_level_data/'

sfr_dir = gwfm_dir+'/SFR_data/'

# %%
out_dir = join(proj_dir, 'output')
fig_dir = join(proj_dir, 'figures')


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

upscale = 4 
upscale_txt = 'upscale'+str(upscale)+'x_'
# model_nam = 'inset_oneto_denier'
model_nam = 'oneto_denier_'+upscale_txt+'2014_2018'

base_model_ws = join(loadpth,model_nam)

# all_model_ws = join(loadpth, 'parallel_oneto_denier')
all_model_ws = join(loadpth, 'parallel_'+model_nam)

# may want to skip loading rch, evt and wel which take up a lot of memory with stress period data
load_only = ['DIS','UPW','SFR','OC']
m = flopy.modflow.Modflow.load('MF.nam', model_ws= base_model_ws, 
                                exe_name='mf-owhm.exe', version='mfnwt')


# %%
homogeneous_ws = join(loadpth, 'oneto_denier_homogeneous_2014_2018')


# %%
m_ver = ''

# %%
grid_dir = join(gwfm_dir, 'DIS_data/streambed_seepage/grid')
grid_fn = join(grid_dir, 'inset_oneto_denier','rm_only_grid.shp')
grid_p = gpd.read_file(grid_fn)
grid_p.crs='epsg:32610'
m_domain = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [grid_p.unary_union], crs=grid_p.crs)

# %%
XSg = pd.read_csv(join(base_model_ws,'04_XSg_filled.csv'))
XSg = gpd.GeoDataFrame(XSg, geometry = gpd.points_from_xy(XSg.Easting, XSg.Northing), crs='epsg:32610')

drop_iseg = XSg[~XSg['Logger Location'].isna()].iseg.values
# overwrite SFR segment/reach input relevant to seepage
# sensor_dict = pd.read_csv(join(model_ws, 'sensor_xs_dict.csv'), index_col=0)
# XS_params = sensor_dict.join(params.set_index('Sensor'), on='Sensor')

# %%
grid_sfr_all = pd.DataFrame()
for r in np.arange(0,100): #100
    folder = 'realization'+ str(r).zfill(3)
    # update model workspace so outputs to right directory
    model_ws = join(all_model_ws, folder)
    grid_sfr = pd.read_csv(model_ws+'/grid_sfr.csv',index_col=0)
    grid_sfr = grid_sfr.drop(columns=['node','geometry','node.1'])
    grid_p_sfr = grid_p.set_index(['row','column']).loc[list(zip(grid_sfr.i+1,grid_sfr.j+1))].reset_index(drop=True)
    grid_sfr = pd.concat((grid_p_sfr,grid_sfr),axis=1)
    grid_sfr_all = pd.concat((grid_sfr_all, grid_sfr.assign(realization=r)))
grid_sfr_all = grid_sfr_all[~grid_sfr_all.iseg.isin(drop_iseg)]

# %%
sfrdf = pd.DataFrame(m.sfr.reach_data)
grid_sfr = grid_p.set_index(['row','column']).loc[list(zip(sfrdf.i+1,sfrdf.j+1))].reset_index(drop=True)
grid_sfr = pd.concat((grid_sfr,sfrdf),axis=1)

# characterize streambed into different hydrofacies
tprogs_quants = np.array([0.590, 0.155, 0.197, 0.058]).cumsum()
vka_quants = grid_sfr.strhc1.quantile(tprogs_quants)
vka_quants.index=['mud','sandy mud','sand','gravel']
grid_sfr['facies'] = 'mud'
for n in np.arange(0,len(vka_quants)-1):
    grid_sfr.loc[grid_sfr.strhc1 > vka_quants.iloc[n],'facies'] = vka_quants.index[n+1]

# add color for facies plots
gel_color = pd.read_csv(join(gwfm_dir,'UPW_data', 'mf_geology_color_dict.csv'), comment='#')
gel_color.geology = gel_color.geology.str.lower()
grid_sfr = grid_sfr.join(gel_color.set_index('geology')[['color']], on='facies')
# remove stream segments for routing purposes only
grid_sfr = grid_sfr[~grid_sfr.iseg.isin(drop_iseg)]

# %%
strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')


# %% [markdown]
# # Groundwater elevation review
# Summarize the mean groundwater elevation in the floodplain, below the channel, in the domain
# - it takes a long time to load the hdobj and longer to load the head timeseries which might make this almost impractical.

# %%
def load_r(r, all_model_ws):
    folder = 'realization'+ str(r).zfill(3)
    # update model workspace so outputs to right directory
    model_ws = join(all_model_ws, folder)
    # upw_r = flopy.modflow.ModflowUpw.load(model_ws+'/MF.upw', model=m)  
    hdobj = flopy.utils.HeadFile(model_ws+'/MF.hds')
    return(hdobj)



# %%
lakarr = m.lak.lakarr.array

# %%
# find where lake existed
lak_lay, lak_row, lak_col = np.where(lakarr[0]==1)
lak_kij = pd.DataFrame(np.transpose(np.where(m.lak.lakarr.array[0]==1)), columns=['k','i','j'])
# get first layer below lake cells
lak_kij = (lak_kij.groupby(['i','j']).max()+1).reset_index()
# create tuples for sampling
lak_idx = list(zip(lak_kij.k, lak_kij.i, lak_kij.j))


# %%
def get_lak_head(hdobj, lak_idx):
    """
    Return the spatially averaged head for the maximum head at the input locations (idx)
    hdobj: flopy head object
    idx: list of tuples as (layer, row, column)
    """
    # get heads under the lake
    lak_ts = hdobj.get_ts(lak_idx)
    lak_ts_df = pd.DataFrame(lak_ts, columns=['totim']+lak_idx)
    lak_ts_df = lak_ts_df.set_index('totim')
    lak_ts_df = lak_ts_df.melt(ignore_index=False)
    lak_ts_df[['k','i','j']] = lak_ts_df.variable.tolist()
    lak_ts_df = lak_ts_df.drop(columns='variable') # drop to speed up groupby
    lak_head = lak_ts_df.groupby(['totim','i','j']).max().groupby('totim').mean()
    return lak_head


# %%
rewrite=True
# extract lake heads for the time series
if rewrite:
    lak_head_all = pd.DataFrame()
    for r in np.arange(0,100):
        if r % 10 ==0:
            print(r, end=', ')
        hdobj = load_r(r, all_model_ws)
        lak_head = get_lak_head(hdobj, lak_idx)
        lak_head_all = pd.concat((lak_head_all, lak_head[['value']].rename(columns={'value':str(r)})),axis=1)
        lak_head_all.to_csv(join(out_dir, 'lak_head_timeseries.csv'))
else:
    lak_head_all = pd.read_csv(join(out_dir, 'lak_head_timeseries.csv'))


# %%

# %% [markdown]
# # Obs checking

# %%
def nse(targets,predictions):
    return 1-(np.sum((targets-predictions)**2)/np.sum((targets-np.mean(predictions))**2))

# hob metadata
rm_grid = pd.read_csv(join(proj_dir, 'mw_hob_cleaned.csv'))


# %%
def clean_hob(model_ws):
    hobout = pd.read_csv(join(model_ws,'MF.hob.out'),delimiter=r'\s+', header = 0,names = ['sim_val','obs_val','obs_nam'],
                         dtype = {'sim_val':float,'obs_val':float,'obs_nam':object})
    hobout[['Sensor', 'spd']] = hobout.obs_nam.str.split('p',n=2, expand=True)
    hobout['kstpkper'] = list(zip(np.full(len(hobout),0), hobout.spd.astype(int)))
    hobout.loc[hobout.sim_val.isin([-1e30, -999.99,-9999]), 'sim_val'] = np.nan
    hobout = hobout.dropna(subset='sim_val')
    hobout = hobout.join(dt_ref.set_index('kstpkper'), on='kstpkper')
    hobout['error'] = hobout.obs_val - hobout.sim_val
    hobout['sq_error'] = hobout.error**2
    return(hobout)


# %%
# 
sum_stats = pd.DataFrame(columns=['r2','RMSE','NSE'], dtype=np.float64)
mw_stats = pd.DataFrame(columns=['realization','SOSE','RMSE','NSE'], dtype=np.float64)
hob_err_chk = pd.DataFrame()
for t in np.arange(0,100):
    model_ws = join(all_model_ws, 'realization'+ str(t).zfill(3))
    hobout = clean_hob(model_ws)
    # removing oneto ag because of large depth offset
    hobout = hobout[hobout.Sensor != 'MW_OA']
    hob_3m = hobout.set_index('dt').groupby('Sensor').resample('3MS').mean(numeric_only=True).reset_index('dt')
    hob_err_chk = pd.concat((hob_err_chk, hob_3m.groupby('dt').mean()))
    # summary stats by well
    mw_stats['realization'] = t
    for s in hobout.Sensor.unique():
        df_s = hobout[hobout.Sensor==s]
        mw_stats.loc[s,'SOSE'] = hobout[['Sensor','sq_error']].groupby('Sensor').sum()
        mw_stats.loc[s,'r2'] = r2_score(df_s.obs_val, df_s.sim_val)
        mw_stats.loc[s,'RMSE'] = mean_squared_error(df_s.obs_val, df_s.sim_val, squared=True)
        mw_stats.loc[s,'NSE'] = nse(df_s.obs_val, df_s.sim_val)

    # summary statistics
    sum_stats.loc[t,'r2'] = r2_score(hobout.obs_val, hobout.sim_val)
    sum_stats.loc[t,'RMSE'] = np.sqrt(hobout.sq_error.sum()/len(hobout))
    sum_stats.loc[t,'NSE'] = nse(hobout.obs_val, hobout.sim_val)


# %%
# filter out realizations who haven't finished running yet
stats_done = sum_stats[sum_stats.NSE!=sum_stats.NSE.min()].copy()

# %%
stats_done.to_csv(join(out_dir, 'hob_fit_stats.csv'))


# %%
# review to see if error is generally similar between realizations
# review hydrographs for realization with worst error
fig,ax = plt.subplots(1,2, figsize=(12,4))
stats_done.plot(y='NSE', ax=ax[0])
stats_done.plot(y='RMSE', ax=ax[1])


# %%
# identify the 10 realizations with the best accuracy
# calculate best score, r2 is tiebreak
stats_done['score'] = (stats_done.NSE >= stats_done.NSE.quantile([0.9]).values[0]).astype(float)
stats_done.score += (stats_done.RMSE <= stats_done.RMSE.quantile([0.1]).values[0]).astype(float)
stats_done.score += (stats_done.r2 >= stats_done.r2.quantile([0.9]).values[0]).astype(float)*0.25
# pull 10 best realizations 
best_realizations = stats_done[stats_done.score >= stats_done.score.quantile([0.9]).values[0]]
print('best realizations', best_realizations.index)
best_realizations.to_csv(join(proj_dir,upscale_txt+'top_10_accurate_realizations.csv'))


# %%
# it seems like in the summer the average error is -1 to -2.75 m with it more sever in the drought

hob_err_chk = hob_err_chk.assign(month=hob_err_chk.index.month, year = hob_err_chk.index.year)
hob_err_chk[hob_err_chk.month==7].boxplot(by='year', column='error')

# %%
# check hydrographs with worst error
# t = sum_stats['RMSE'].idxmax()
t = sum_stats['RMSE'].idxmin()
# approximate median location
# t = sum_stats.sort_values('RMSE').iloc[int(len(sum_stats)/2)].name
# t=11
# t = 45
print(t)
print(sum_stats.loc[t])
def mak_hob_long(t):
    hobout = clean_hob(join(all_model_ws, 'realization'+ str(t).zfill(3)))
    # removing oneto ag because of large depth offset
    hobout = hobout[hobout.Sensor != 'MW_OA']
    hob_long = hobout.melt(id_vars=['dt', 'Sensor'],value_vars=['sim_val','obs_val'], value_name='gwe')
    return(hob_long)
hob_long = mak_hob_long(t)
# hob_long

# %%
hob_h = clean_hob(homogeneous_ws)
t = sum_stats['RMSE'].idxmin()
hob_min = clean_hob(join(all_model_ws, 'realization'+ str(t).zfill(3)))
t=sum_stats['RMSE'].idxmax()
hob_max = clean_hob(join(all_model_ws, 'realization'+ str(t).zfill(3)))
t= sum_stats.sort_values('RMSE').iloc[int(len(sum_stats)/2)].name
hob_med = clean_hob(join(all_model_ws, 'realization'+ str(t).zfill(3)))


# %%
cols = ['dt', 'sim_val','Sensor']
hob_comp = pd.concat((
    hob_med[cols].assign(var='Median'),
    hob_max[cols].assign(var='Max'),
    hob_min[cols].assign(var='Min'),
    hob_h[cols].assign(var='Homogeneous'),
    hob_h[['dt','obs_val','Sensor']].rename(columns={'obs_val':'sim_val'}).assign(var='Observations')
    ))


# %%
nx=4
wells = hob_med.Sensor.unique()
ny = int(np.round(len(wells)/nx))
fig,ax = plt.subplots(4,4, sharex=True, sharey=True, figsize=(6.5, 8), dpi=300)
for n, w in enumerate(wells):
    ax_n = ax[int(n/ny), n%ny]
    hob_med[hob_med.Sensor==w].plot(x='dt',y='sim_val', ax=ax_n, legend=False)
    hob_max[hob_max.Sensor==w].plot(x='dt',y='sim_val', ax=ax_n, legend=False)
    hob_min[hob_min.Sensor==w].plot(x='dt',y='sim_val', ax=ax_n, legend=False)
    hob_h[hob_h.Sensor==w].plot(x='dt',y='sim_val', ax=ax_n, legend=False)
    hob_h[hob_h.Sensor==w].plot(x='dt',y='obs_val', ax=ax_n, legend=False, marker='x', linestyle='', color='black',markersize=0.5)
fig.supylabel('Groundwater Elevation (m)')
fig.supxlabel('Date')

fig.tight_layout(h_pad=0.1, w_pad=-0.5)

for n in np.arange(0,nx):
    ax_n = ax[-1,n]
    ax_n.set_xlabel(None)
    ax_n.set_xticks(pd.date_range(strt_date, end_date, freq='AS'), 
                         pd.date_range(strt_date, end_date, freq='AS').year.astype(str).values, rotation=45)
    ax_n.set_xticks(pd.date_range(strt_date, end_date, freq='3MS'), minor=True)

# %%
# in the wrost case the dynamics match but the magnitude is off (levels start much too low)
# import seaborn as sns
# g = sns.relplot(hob_long, x='dt',y='gwe',col='Sensor',hue='variable', col_wrap=4);




# %% [markdown]
# ## Water Budget check

# %%
# manual columns
wb_out_cols  =['WEL_OUT','ET_OUT','GHB_OUT','SFR_OUT','LAK_OUT']
wb_in_cols = ['RCH_IN','GHB_IN','SFR_IN','LAK_IN']

# %%
wb_all = pd.DataFrame()
for t in np.arange(0,100):
    model_ws = join(all_model_ws, 'realization'+ str(t).zfill(3))
    # load summary water budget
    wb = pd.read_csv(model_ws+'/flow_budget.txt', delimiter=r'\s+')
    # wb = pd.read_csv(loadpth+'/oneto_denier_upscale8x_2014_2018'+'/flow_budget.txt', delimiter=r'\s+')
    wb['kstpkper'] = list(zip(wb.STP-1,wb.PER-1))
    wb = wb.merge(dt_ref, on='kstpkper')
    wb['realization'] = t
    wb_all = pd.concat((wb_all, wb))


# %%
# wb_plt.mean(axis=1)
# need to update homogeneous run for layering and
wb_h = pd.read_csv(homogeneous_ws+'/flow_budget.txt', delimiter=r'\s+')
wb_h['kstpkper'] = list(zip(wb_h.STP-1,wb_h.PER-1))
wb_h = wb_h.merge(dt_ref, on='kstpkper')

# %%
fig,ax = plt.subplots(5,1, sharex=True, layout='constrained')
for n, wb_n in enumerate(wb_out_cols):
    wb_plt = wb_all.pivot_table(index='dt',columns='realization',values=wb_n)
    wb_plt.plot(legend=False, color='gray', ax=ax[n]) 
    wb_plt.mean(axis=1).plot(color='red',linestyle='--',ax=ax[n])
    wb_h.plot(x='dt', y=wb_out_cols[n], color='black',linestyle='--',ax=ax[n])
    ax[n].set_ylabel(wb_out_cols[n].split('_')[0])

# %%
wb_chk_plt = wb_all[wb_all.realization==t]
fig,ax= plt.subplots(3,1, sharex=True)
wb_chk_plt.plot(y='PERCENT_ERROR', ax=ax[0])
wb_chk_plt.plot(y=wb_out_cols, ax=ax[1], legend=True)
wb_chk_plt.plot(y=wb_in_cols, ax=ax[2], legend=True)

