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
git_dir = join(doc_dir,'GitHub')
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(join(git_dir,'flopy'))
import flopy

py_dir = join(git_dir,'CosumnesRiverRecharge/python_utilities')
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



# %% [markdown]
# # Stream seepage plots
# Grouping must be done by segment and reach or by total distance (m).  
#

# %%
# slightly different version of clean_sfr_df is used here to help with realization comparison
# import mf_utility
# from importlib import reload
# reload(mf_utility)
# from mf_utility import clean_sfr_df
from OD_utility import clean_sfr_df

# %%
rewrite = False
# rewrite = True

# %%
# from importlib import reload
# import OD_utility
# reload(OD_utility)
# from OD_utility import clean_sfr_df

# %%
sfrdf =  clean_sfr_df(base_model_ws, drop_iseg, dt_ref)
h_sfrdf =  clean_sfr_df(homogeneous_ws, drop_iseg, dt_ref)




# %%
# standardize grouping values
wy_vals = sfrdf.WY.unique()
facies_vals = ['Mud','Sandy Mud','Sand','Gravel']


# %% [markdown]
# # Aggregate stream data
# 1. Aggregate across all reaches  but save all dates
# 2. Aggregate across dates but save all reaches
#

# %% [markdown]
# ## Gradient plots (spatial)
# 1. Seepage averaged across the year (or between dry and wet season) and the rows could be realizations instead which would help indicate consistency across realizations  
# 2. Heat map of columns with stream reaches, rows of dates and the color blue to red for gaining or losing with the seepage averaged across all realizations
#
# When the gradient is greater than 1 we know we have disconnected conditions, I need to represent the count of days where the system is connected.
# - Since I've focused more baseflow/seepage and streamflow/flowing days the connected/gaining gradient plots tend to repeat the information on baseflow/seepage and are usual lower correlations because they are slightly less impacted by the coarse vs fine conductivity.

# %%
if rewrite:
    # aggregate data by facies and sum to review seepage over time
    t0 = time.time()

    sfr_facies_all = pd.DataFrame()
    for t in np.arange(0,100):
        model_ws = join(all_model_ws, 'realization'+ str(t).zfill(3))
        # remove stream segments for routing purposes
        sfrdf =  clean_sfr_df(model_ws, drop_iseg, dt_ref, 'MF')
        # summing by facies makes sense for seepage
        sfr_facies_sum = sfrdf.groupby(['dt','facies']).sum(numeric_only=True)
        sfr_facies_sum['realization'] = t
        # count number of facies if needed later for explaining rates
        sfr_facies_sum['num_facies'] = sfrdf.groupby(['dt','facies']).count().iloc[:,0].values
        sfr_facies_all = pd.concat((sfr_facies_all, sfr_facies_sum))
    # check time
    t1 = time.time()
    print('Time: %.2f min' % ((t1-t0)/60))
    # save output to speed up reuse
    # sfr_facies_all.to_csv(join(out_dir, 'sfrdf_facies_sum.csv'))
    sfr_facies_all.reset_index(level='facies').to_hdf(join(out_dir, 'sfrdf_facies_sum.hdf5'), key='dt', complevel=4)

# %%
sfr_facies_all = pd.read_hdf(join(out_dir, 'sfrdf_facies_sum.hdf5'))
# sfr_facies_all = sfr_facies_all.set_index('facies', append=True)

# %%
# fig,ax = plt.subplots()

# days connected
# grp_col = 'Total distance (m)'
# df_plt = sfrdf.groupby(grp_col).sum(numeric_only=True)[['connected']]
# ax_conn = df_plt.plot(legend=False)
# plt.ylabel('Connected Days')
# plt.xlabel(grp_col)

# for f in sfrdf.facies.unique():
#     ax_conn.fill_between(sfrdf[grp_col], 0, df_plt.max(), where = sfrdf.facies==f,
#                     color=gel_color.loc[gel_color.geology==f,'color'], alpha=0.5)

# %% [markdown]
# ## Notes
# In these I need to clarify the proportion of the streambed (length or area) that is each facies to show that despite being only a small percent of the streambed sand and gravel make up a significant portion of recharge and baseflow. 

# %% [markdown]
# ## Seepage plots (temporal)
# Aggregate by facies to plot cumulative seepage (by time) to help show variability caused by geology.

# %%
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = [
    Line2D([0], [0], color='gray', label='Heterogeneous\nRealizations'),
    Line2D([0], [0], color='black', linestyle='--',label='Homogeneous\nCase'),
    Line2D([0], [0], color='red', linestyle='--', label='Heterogeneous\nMean')
]

ts_lgd = [
    Line2D([0], [0], color='gray', label='100 Heterogeneous\nRealizations'),
    # Line2D([0], [0], color='black', linestyle='--',label='Homogeneous\nCase'),
    # Line2D([0], [0], color='red', linestyle='--', label='Heterogeneous\nMean'),
    Line2D([0], [0], color='black', label='Median'),
    Line2D([0], [0], color='blue', label='Maximum'),
    Line2D([0], [0], color='red', label='Minimum'),
]

# %% [markdown]
# Despite being only a small percentage of the stream reaches, the sand and gravel produce a significant portion of the stream seepage.

# %%
value = 'Qbase'
ylabel  = 'Baseflow ($m^3/d$)'
def plt_dt_facies(value, ylabel):
    # original plotting was 12 by8 but then text is very small
    fig,ax = plt.subplots(2,2, figsize=(8,5.3), sharex=True, sharey=True, layout='constrained')#
    facies = ['Mud','Sandy Mud','Sand','Gravel']
    for t in np.arange(0,100):
        sfr_facies_sum = sfr_facies_all[sfr_facies_all.realization==t]
        for n,f in enumerate(facies):
            ax_n = ax[int(n/2), n%2]
            df_plt = sfr_facies_sum[sfr_facies_sum.facies==f]
            if df_plt.shape[0]>0:
                df_plt.plot(y=value, ax=ax_n, legend=False, color='gray')
    # plot homogeneous case
    h_sfr_facies_sum = h_sfrdf.groupby(['dt','facies']).sum(numeric_only=True)
    # plot mean of heterogeneous
    sfr_facies_mean = sfr_facies_all.groupby(['dt', 'facies']).mean().reset_index('facies')
    # set axis labels
    for n,f in enumerate(facies):
        ax_n = ax[int(n/2), n%2]
        h_sfr_facies_sum.reset_index('facies').plot(y=value, ax=ax_n, legend=False, color='black', linestyle='--')
        sfr_facies_mean[sfr_facies_mean.facies==f].plot(y=value, ax=ax_n, legend=False, color='red', linestyle='--')
        ax_n.set_title(f)
        ax_n.set_yscale('log')
        ax_n.set_ylabel(ylabel)
        ax_n.set_xlabel('Date')
    # add figure legend
    fig.legend(handles=legend_elements, loc='outside upper center', ncol = 3)


# %%
# value = 'Qbase'
# plt_dt_facies(value, ylabel)
# ax.legend(handles=legend_elements, loc='lower right')


# %% [markdown]
# On an individual realization level the plots grouping baseflow by facies are interesting but on a whole they are not because among realizations there is enough variability to make patterns appear similar. 

# %%
# value = 'Qrech'
# ylabel  = 'Stream Seepage ($m^3/d$)'
# plt_dt_facies(value, ylabel)


# %% [markdown]
# ## Streamflow
# No need to aggregate by facies, instead show impact at downstream end in terms of time step and cumulative

# %%

grp_col = 'Total distance (m)'


# %%
if rewrite:
    sfr_last_all = pd.DataFrame()
    for t in np.arange(0,100):
        model_ws = join(all_model_ws, 'realization'+ str(t).zfill(3))
        sfrdf =  clean_sfr_df(model_ws, drop_iseg, dt_ref, 'MF')
        # plot from last segment (shows cumulative effects)
        sfr_last = sfrdf[sfrdf[grp_col] == sfrdf[grp_col].max()].copy()
        sfr_last['realization'] = t
        sfr_last_all = pd.concat((sfr_last_all, sfr_last))

    # save data
    sfr_last_all.to_hdf(join(out_dir, 'sfrdf_last_seg.hdf5'), key='dt', complevel=4)

# %%
sfr_last_all = pd.read_hdf(join(out_dir, 'sfrdf_last_seg.hdf5'))


# %% [markdown]
# The plots of streambed seepage and streamflow at the outlet present the variability among realizations and provide context for the results we will interpret. They explain we have a dominantly losing system where the mangitude of seepage is shifted up or down due to heterogeneity and the streamflow is also shifted up or down but this is most noticeable during lower flows when seepage has a greater proportional impact on streamflow.  
# - Plot the standard deviation among the results for year/season because this will show how far off a difference of heterogeneity can affect interpretation, e.g., if we don't include heterogeneity we might miss the magnitude of flow imapcts.

# %%
# fig,ax = plt.subplots(3,1, figsize=(10,6.3), sharex=True, sharey=False, layout='constrained',dpi=600)

def plt_ts_quants(sfr_facies_all, sfr_last_all, ax):
    # q_ind = [0,24, 74,99]
    q_ind = [0,49, 99]
    seg_mean = sfr_facies_all.groupby('realization').mean(numeric_only=True)

    rech_quant = seg_mean.sort_values('Qrech').iloc[q_ind].index.values
    base_quant = seg_mean.sort_values('Qbase').iloc[q_ind].index.values
    flow_mean = sfr_last_all.groupby('realization').mean(numeric_only=True)
    flow_quant = flow_mean.sort_values('Qin').iloc[q_ind].index.values
    
    colors = ['red', 'black','blue']
    # plot quantiles
    for n in np.arange(0,3):
        t = rech_quant[n]
        sfr_seg = sfr_facies_all[sfr_facies_all.realization==t].groupby('dt').sum(numeric_only=True)
        sfr_seg.plot(y='Qrech', ax=ax[0], legend=False, color=colors[n])
        t = base_quant[n]
        sfr_seg = sfr_facies_all[sfr_facies_all.realization==t].groupby('dt').sum(numeric_only=True)
        sfr_seg.plot(y='Qbase', ax=ax[1], legend=False, color=colors[n])
        t = flow_quant[n]
        sfr_last = sfr_last_all[sfr_last_all.realization==t]
        sfr_last.plot(y='Qin', ax=ax[-1], legend=False, color=colors[n])


# %%

def plt_ts_100(sfr_facies_all, sfr_last_all, ax):
    for t in np.arange(0,100):
        sfr_last = sfr_last_all[sfr_last_all.realization==t]
        sfr_last.plot(y='Qin', ax=ax[-1], legend=False, color='gray')
        sfr_seg = sfr_facies_all[sfr_facies_all.realization==t].groupby('dt').sum(numeric_only=True)
        # sfr_seg.plot(y='Qaquifer', ax=ax[0], legend=False, color='gray')
        sfr_seg.plot(y='Qrech', ax=ax[0], legend=False, color='gray')
        sfr_seg.plot(y='Qbase', ax=ax[1], legend=False, color='gray')
    return None


# plt_ts_quants(sfr_facies_all, sfr_last_all, ax)
# # # plot homogeneous case
# h_sfr_last = h_sfrdf[h_sfrdf[grp_col]==h_sfrdf[grp_col].max()]
# h_sfr_last.plot(y='Qin', ax=ax[-1], legend=False, color='black',linestyle='--')
# # h_sfrdf.groupby('dt').sum(numeric_only=True).plot(y='Qaquifer', ax=ax[0], legend=False, color='black', linestyle='--')
# h_sfrdf.groupby('dt').sum(numeric_only=True).plot(y='Qrech', ax=ax[0], legend=False, color='black', linestyle='--')
# h_sfrdf.groupby('dt').sum(numeric_only=True).plot(y='Qbase', ax=ax[1], legend=False, color='black', linestyle='--')
# # plot mean of heterogeneous
# sfr_last_mean = sfr_last_all.groupby('dt').mean(numeric_only=True)
# sfr_last_mean.plot(y='Qin', ax=ax[-1], legend=False, color='red',linestyle='--')
# sfr_sum_mean = sfr_facies_all.groupby(['dt', 'realization']).sum(numeric_only=True).groupby('dt').mean()
# # sfr_sum_mean.plot(y='Qaquifer', ax=ax[0], legend=False, color='red', linestyle='--')
# sfr_sum_mean.plot(y='Qrech', ax=ax[0], legend=False, color='red', linestyle='--')
# sfr_sum_mean.plot(y='Qbase', ax=ax[1], legend=False, color='red', linestyle='--')

def plt_ts_axes(ax):
    # set axis labels
    ax[-1].set_xlabel('Date')
    # ax[0].set_ylabel('Total Seepage ($m^3/d$)')
    ax[0].set_ylabel('Total\nRecharge ($m^3/d$)')
    ax[1].set_ylabel('Total\nBaseflow ($m^3/d$)')
    ax[-1].set_ylabel('Daily\nStreamflow ($m^3/d$)')
    
    # need log scale or peaks wash out other data
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    ax[-1].set_yscale('log')
    return None

# t0 = time.time()
# fig,ax = plt.subplots(3,1, figsize=(10,6.3), sharex=True, sharey=False, layout='constrained',dpi=600)

# plt_ts_100(sfr_facies_all, sfr_last_all, ax)

# plt_ts_axes(ax)
# fig.legend(handles=legend_elements, loc='outside upper center', ncol = 3)
# # fig.legend(handles=ts_lgd, loc='outside upper center', ncol = 3)

# # fig.tight_layout()
# # plt.savefig(join(fig_dir, 'time_series_seepage_flow.png'), bbox_inches='tight')
# t1 = time.time()
# print('Time: %.2f min' % ((t1-t0)/60))


# %%
fig,ax = plt.subplots(3,1, figsize=(10,6.3), sharex=True, sharey=False, layout='constrained',dpi=600)

plt_ts_100(sfr_facies_all, sfr_last_all, ax)
plt_ts_quants(sfr_facies_all, sfr_last_all, ax)

plt_ts_axes(ax)
fig.legend(handles=ts_lgd, loc='outside upper center', ncol = 4)

# fig.tight_layout()
# plt.savefig(join(fig_dir, 'time_series_seepage_flow_quantiles.png'), bbox_inches='tight')



# %% [markdown]
# ### In-text numbers, KEEP

# %%
# how does baseflow persist into the summer by WY
# sfr_seg.r.groupby(['realization','WY'])
sfr_base = sfr_facies_all[sfr_facies_all.Qbase !=0].copy()
sfr_base['WY'] = sfr_base.index.year
sfr_base.loc[sfr_base.index.month>=10, 'WY'] +=1

sfr_base_max = sfr_base.reset_index().groupby(['realization','WY']).max()
sfr_base_max[['dt']].groupby('WY').mean()
# sfr_base_max[['dt']].groupby('WY').std() # std deviation is pretty constant 16-20 days

# %%
# introduce the concept that there is  large variability in days with flow between realizations
year_flow = sfr_last_all.groupby(['realization','WY']).sum(numeric_only=True)
year_flow.flowing.groupby('WY').median(), year_flow.flowing.groupby('WY').std()

# %%
## more detailed analysis investigating stream drying in 2017 summer
plt_strt = pd.to_datetime('2017-7-1')
plt_end = pd.to_datetime('2017-9-30')
plt_dates = pd.date_range(plt_strt, plt_end)

zero_flow = sfr_last_all.loc[plt_dates].groupby('realization').mean(numeric_only=True)
zero_flow = sfr_last_all.loc[plt_dates].groupby('realization').sum(numeric_only=True)
# zero_flow[zero_flow.Qout>0]
# zero_flow.hist('Qout')
# zero_flow.flowing.hist()
flow_frac = (zero_flow.flowing/len(plt_dates))
flow_frac.quantile([0,.5,1]), flow_frac.std()
# zero_flow.flowing.hist(bins=np.arange(0, len(plt_dates)+10, 10))
print('Realizations with >50%% of days with no flow: %i' %((zero_flow.flowing<0.5*len(plt_dates)).sum()))
print('Realizations with continuous flow: %i' %(zero_flow.flowing==len(plt_dates)).sum())
# len(plt_dates)
# we want to show the min, max, median days without flow

# %% [markdown]
# The standard deviation for streamflow can be larger than the flow itself because in the later summer the flow gets very small. If I plot the std dev as a fraction of the mean then the error would be 100-1000%. For streamflow it might make more sense to plot the variability in days with flow for a year because this has a more concrete impact on ecosystems as under predicting the flow availability by 1-2 months is worse than 1-2 days.
#
# WY2016 is the only year that starts initially with no flow in the channel (where there is no seepage from July 2015- November 2015). It's an interesting year because there is an initial flow event that makes it to the outlet but then upstream drying occurs.
#

# %%
# fig,ax = plt.subplots(2,1, figsize=(8,5.3), sharex=True, sharey=False, layout='constrained',)
sfr_sum_std = sfr_facies_all.groupby(['dt', 'realization']).sum(numeric_only=True).groupby('dt').std()
# (sfr_sum_std['Qaquifer']/sfr_sum_mean['Qaquifer']).plot()
# (sfr_sum_std['Qaquifer']/sfr_sum_mean['Qaquifer']).mean()
# (sfr_sum_std['Qrech']/sfr_sum_mean['Qrech']).plot()
# (sfr_sum_std['Qrech']/sfr_sum_mean['Qrech']).mean()

# (sfr_sum_std['Qaquifer']).plot()
# sfr_last_std =  sfr_last_all.groupby('dt').std(numeric_only=True)
# (sfr_last_std['Qout']/sfr_last_mean['Qout']).plot(ax=ax[1])
# (sfr_last_std['Qout']).plot(ax=ax[1])

# ax[0].set_yscale('log')
# ax[1].set_yscale('log')


# %%
# (sfr_sum_std['Qbase']/sfr_sum_mean['Qbase']).plot()
# (sfr_sum_std['Qbase']/sfr_sum_mean['Qbase']).mean()

# %% [markdown]
# ## Scatter plots
# **Goal**: relate the heterogeneity (e.g., number of coarse bodies connecting to the stream) to the number of days with streamflow.
#
# Need to save data as an aggregated dataframe to run trend analysis

# %%
if rewrite:
    t0 = time.time()

    sfr_3mon_all = pd.DataFrame()
    sfr_yr_sum_all = pd.DataFrame()

    for t in np.arange(0,100):
        model_ws = join(all_model_ws, 'realization'+ str(t).zfill(3))
        # clean sfr output
        sfrdf =  clean_sfr_df(model_ws, drop_iseg, dt_ref, 'MF')
        sfrdf['realization'] = t
        num_coarse = sfrdf['num_coarse'].mean()
        # aggregate to seasonal values, since model starts in october it groups as oct-dec, jan-mar, apr-jun, jul-sep
        sfrdf_mon = sfrdf.resample('3MS').mean(numeric_only=True)
        sfrdf_mon['realization'] = t
        sfrdf_mon['num_coarse'] = num_coarse
        sfr_3mon_all = pd.concat((sfr_3mon_all, sfrdf_mon))
        # aggregate to annual values for each segment, reach
        sfrdf_yr_sum = sfrdf.groupby(['WY','segment', 'reach']).sum(numeric_only=True)
        sfrdf_yr_sum['realization'] = t
        sfrdf_yr_sum['num_coarse'] = num_coarse
        sfr_yr_sum_all = pd.concat((sfr_yr_sum_all, sfrdf_yr_sum))
    # check time
    t1 = time.time()
    print('Time: %.2f min' % ((t1-t0)/60))

    # save output to speed up reuse
    sfr_3mon_all.to_csv(join(out_dir, 'sfrdf_3month_mean.csv'))
    sfr_yr_sum_all.to_csv(join(out_dir, 'sfrdf_annual_sum_by_segment.csv'))

# %%
# load previously made output to speed up reuse
sfr_3mon_all = pd.read_csv(join(out_dir, 'sfrdf_3month_mean.csv'), parse_dates=['dt'], index_col='dt')
# fix issue where this was averaged
sfr_3mon_all.month = sfr_3mon_all.index.month
sfr_yr_sum_all = pd.read_csv(join(out_dir, 'sfrdf_annual_sum_by_segment.csv'), index_col=['WY','segment','reach'])

# %%
grid_sfr['Total distance (m)'] = grid_sfr.rchlen.cumsum()

# fix total distance broken with sum
sfr_dist = grid_sfr[['iseg','ireach','Total distance (m)']].rename(columns={'iseg':'segment','ireach':'reach'})
sfr_yr_sum_all = sfr_yr_sum_all.drop(columns=['Total distance (m)']).reset_index().merge(sfr_dist).set_index(['WY','segment','reach'])


# %%
## troubleshooting changes in output
# old_folder = 'old_segments'
# sfr_3mon_all = pd.read_csv(join(out_dir, old_folder, 'sfrdf_3month_mean.csv'), parse_dates=['dt'], index_col='dt')
# sfr_yr_sum_all = pd.read_csv(join(out_dir, old_folder, 'sfrdf_annual_sum_by_segment.csv'), index_col=['WY','segment'])
# sfr_facies_all = pd.read_hdf(join(out_dir, old_folder,'sfrdf_facies_sum.hdf5'))


# %%
# after resampling if need to re-identify the number of coarse in a stream
coarse_ref = sfr_3mon_all.groupby('realization').mean(numeric_only=True)[['num_coarse']]
# coarse_ref.to_csv(join(proj_dir, 'num_sfr_coarse.csv'))

# %%
# relating coarse segments to days of flow
# median, std = coarse_ref.loc[zero_flow.flowing==len(plt_dates)].median(), coarse_ref.loc[zero_flow.flowing==len(plt_dates)].std()
# print('No flow cessation had %i coarse segments' %median.iloc[0], 'with a std dev of %.2f' %std.iloc[0])
def flow_summary(flw_frac, dry_frac, zero_flow, coarse_ref, plt_dates):
    # look at wetted days
    df = coarse_ref.loc[zero_flow.flowing>=flw_frac*len(plt_dates)]
    if len(df)>0:
        median, std = df.median(), df.std()
        print(len(df),'Realizations with >%i' %(flw_frac*100), '%% flow had %i coarse reaches' %median.iloc[0], 
              'with a std dev of %.2f' %std.iloc[0])    
        # look at drying days
    df = coarse_ref.loc[len(plt_dates)-zero_flow.flowing > dry_frac*len(plt_dates)]
    if len(df)>0:
        median, std = df.median(), df.std()
        print(len(df), 'Realizations with >%i' %(dry_frac*100), '%% dry had %i coarse reaches' %median.iloc[0], 
              'with a std dev of %.2f' %std.iloc[0])
    return None
    
flow_summary(.75, 0.75, zero_flow, coarse_ref, plt_dates)
flow_summary(1, 0.5, zero_flow, coarse_ref, plt_dates)



# %%
from matplotlib import cm

#normalize item number values to colormap
norm = mpl.colors.Normalize(vmin=0, vmax=99)
# plt.cmap='gray'

# %%
def plt_wy_seg(value, ylabel, grp_col, log=False):
    grp_col = 'Total distance (m)'
    # grp_col = ['segment','reach']
    grp_wy = np.append(['WY'], grp_col).tolist()
    fig,ax = plt.subplots(2,2, figsize=(8,8), sharex=True, sharey=True, layout='constrained')
    wy_unique = sfr_yr_sum_all.index.get_level_values('WY').unique()

    for t in np.arange(0,100):
        sfr_yr_sum = sfr_yr_sum_all[sfr_yr_sum_all.realization==t].reset_index('WY').set_index(grp_col)
        for n,f in enumerate(wy_unique):
            ax_n = ax[int(n/2), n%2]
            df_plt = sfr_yr_sum[sfr_yr_sum.WY==f]
            if df_plt.shape[0]>0:
                df_plt.plot(y=value, ax=ax_n, legend=False, 
                            color='gray',
                            # color=cm.gray(norm(t)),alpha=0.7,
                                         )

    # plot homogeneous case
    h_sfr_yr_sum = h_sfrdf.groupby(grp_wy).sum(numeric_only=True).reset_index('WY')
    # plot mean of heterogeneous
    sfr_yr_sum_mean = sfr_yr_sum_all.groupby(grp_wy).mean().reset_index('WY')
    # set axis labels
    for n,f in enumerate(wy_unique):
        ax_n = ax[int(n/2), n%2]
        h_sfr_yr_sum[h_sfr_yr_sum.WY==f].plot(y=value, ax=ax_n, legend=False, color='black', linestyle='--')
        sfr_yr_sum_mean[sfr_yr_sum_mean.WY==f].plot(y=value, ax=ax_n, legend=False, color='red', linestyle='--')
        ax_n.set_title(f)
        ax_n.set_ylabel(ylabel)
        ax_n.set_xlabel('Total distance (m)')
        if log:
            ax_n.set_yscale('log')
#     fig.tight_layout()
#     fig.legend(handles=legend_elements, loc='center', bbox_to_anchor=[0.55, 1.03], ncol = 3)
    fig.legend(handles=legend_elements, loc='outside upper center', ncol = 3)
    return ax


# %%
# h_sfr_yr_sum =  h_sfrdf.groupby(['WY', grp_col]).sum(numeric_only=True).reset_index('WY')
# h_sfr_yr_sum[h_sfr_yr_sum.WY==2016].plot(y=value, legend=False, color='black', linestyle='--')

# %%
# grid_sfr_all.groupby('iseg').mean(numeric_only=True).strhc1.plot(label='Mean')
# grid_sfr_all.groupby('iseg').std(numeric_only=True).strhc1.plot(label='Std Dev')
# plt.legend()
# # plt.yscale('log')
# plt.ylabel('Streambed Conductivity ($m/day$)')
# plt.xlabel('Segment')
# there is not a particularly helpful pattern from the mean/median and std devation of the longitudinal strhc1 except some indication of the patch that is consistenly sand/gravel.

# %% [markdown]
# Plotting the longitudinal profile for the baseflow, connected days, gaining is not very helpful because there is not a clear way to connect the heterogeneity.  
#
# If we plot baseflow with variable gray scale for each realization we start to see how there is a scattering of baseflow for each realization but it isn't enitrely clear.
# - there is still some value in showing the difference from a homogeneous case
#
# A box plot of baseflow by segment is very messy (mostly zeros) means quite a few fliers
#
# If we don't care about presenting the volumes of baseflow so much as the rates than the existence of it so we could do a count of the number of realizations with baseflow?

# %%
# fig,ax = plt.subplots()
# for t in np.arange(0,100):
#     sfr_seg = sfr_facies_all[sfr_facies_all.realization==t].groupby('dt').sum(numeric_only=True)
#     sfr_seg.plot(y='Qbase', ax=ax, legend=False, color=cm.gray(norm(t)),alpha=0.7)
# plt.yscale('log')

# %%
# sfr_base = sfr_yr_sum_all.groupby(['segment', 'reach','realization']).mean().reset_index('realization')
# # sfr_base.plot(y='Qbase')
# value='Qbase'
# fig, ax_n = plt.subplots()
# ax_n.set_yscale('log')
# for t in np.arange(0,100):
#     df_plt = sfr_base[sfr_base.realization==t]
#     df_plt.plot(y=value, ax=ax_n, legend=False, 
#                 # color='gray',
#                 color=cm.gray(norm(t)),alpha=0.7,
#                )
   
# # plot homogeneous case
# # h_sfr_yr_sum = h_sfrdf.groupby(['segment', 'reach']).sum(numeric_only=True).plot(y=value, ax=ax_n, legend=False, color='black', linestyle='--')
# # # plot mean of heterogeneous
# # sfr_yr_sum_mean = sfr_yr_sum_all.groupby(['segment','reach']).mean().plot(y=value, ax=ax_n, legend=False, color='red', linestyle='--')
# plt.xlabel('Segment')
# plt.ylabel('Average Annual Baseflow ($m^3/day$)')

# %%

# value = 'Qbase'
# ylabel = 'Baseflow'
# value = 'connected'
# ylabel = 'Connected Days'
# plt_wy_seg(value, ylabel, log=True)
# ax.legend(handles=legend_elements, loc='lower right')


# %% [markdown]
# The days with flow drops off at segment 31 because half of the flow becomes over bank flow to the floodplain so less water remains in the channel to recharge near the channel. Then flow returns at segment 51.
#
# The number of connected days is higher with the heterogeneity because there is more groundwater available but what is important is that the number of connected days on average doesn't rise above 150-250 on average while days with flow is ideally much closer to 365 (in wet years this happens). So although we see more days of connection with baseflow overall the system loses water sooner than in a homogeneous case.
#
# Also the key conversation between longitudinal days with flow and streamflow at the outlet over time is that by the last segment we see the cumulative effect of the baseflow while intermediate segments may see more or less flow due to high seepage rates and the floodplain.

# %%
# ax[0,0].get_xticks()

# %%
value = 'flowing'
ylabel = 'Days with Flow'
ax = plt_wy_seg(value, ylabel, grp_col)
# plt.savefig(join(fig_dir, 'days_with_flow_by_WY.png'), bbox_inches='tight')

# %%
# sfr_yr_sum_all.groupby(['segment']).std(numeric_only=True).plot(y='flowing')
# sns.relplot(sfr_yr_sum_all.groupby(['WY','segment']).std(numeric_only=True), x = 'segment',y='flowing', col='WY', col_wrap=2)

# %% [markdown]
# All WY have a fairly tight start (low variance) because the flow at the inlet is least impacted yet by differences in geology. In low flow years we don't see a signal from the floodplain connection but rather we see a fairly consistent decline across all realizations in days with flow with the largest variablility toward the outlet. In WY2016 there is a pattern where the peak variability occurs just before the floodplain diversion, which might be explained by the floodplain distributing water for recharge such that it provides a reduction in seepage rates out of the stream, this pattern might be stronger here because 2016 starts with initially dry streamflow conditions. The WY2016 reduction in standard deviation below the floodplain diversion might also be explained by the fact that less water in the channel means less seepage so the channel might be more uniformly dry  
#
# IN WY2017 we see lower variability in the upper reaches because there is so much flow available that even seepage from the coarse facies can't drain it all, but as we get downstream we start to see more variability as the cumulative effect of coarse segments occurs and this is likely impact by the floodplain return flows having more variability. 

# %%
# need to see why water could be leaving the lake in a few scenarios
# g = sns.relplot(sfr_yr_sum_all, x='segment',y='flowing', col='WY', col_wrap=2, hue='realization', #color='gray', 
#            facet_kws={'sharey': False, 'sharex': True})



# %%
# import seaborn as sns

# g = sns.relplot(sfrdf_all, x='num_coarse',y='Qbase', col= 'month', row='WY', color='gray', 
           # facet_kws={'sharey': False, 'sharex': True})



# %% [markdown]
# # Statistics
# 1. Present overview of all parts of the streamflow system that could be reasonably impacted by heterogeneity, streamflow, streambed recharge, stream baseflow, days with flow, stream-aquifer connectivity and maybe gradient. Present on aggregate basis or annual basis with just the correlation values.  
#
# 2. Present the streamflow components with high correlation values on a seasonal and annual basis to show how heterogeneity plays a larger role at different times of the year.
#
# - For streamflow the annual sum or mean should only be taken on the last segment to look at cumulative effects.
# - For streambed recharge and baseflow the annual sum across all segments should be used with either the sum or mean across all days.
# - For stream-aquifer connectivity it would need to be the mean across all segments and the sum across all days. The gaining losing relationship is equal and opposite so only use gaining.
#
# Is it worth looking at how the correlation is on a daily scale?
#
# - Pearson's r absolute value: .5 to 1 is strong, .3 to .5 is moderate, 0 to .3 is weak
#     - italicize r, no leading zero as maximum is 1, two significant digits is common, if reporting significance of the test then also report degrees of freedom (sample size - 2)   
# - Spearman's coefficient is defined as the Pearson coefficient between the rank variables. X, Y are converted to ranks R(X), R(Y) then pearson's r is calculated.
# - Kendall's tau is a test on whether the order of the variables aligns with the values, that is xi > xj and yi > yj to be considered concordant. The statistics then is based on the number of concordant and discordant pairs.
#
# **After reviewing the paper on correlations, the value of Spearman's is that it is appropriate when both variables are skewed and is not affected by extreme values**
# - it's not needed in this case because the value range doesn't have many extremes and there isn't a particularly strong skew although there is a tendency to lower values since most of the time there is less flow in the system.

# %%
# the connected and gaining variables are interesting from a box plot perspective but the correlations tend to be insignificant or NAs
variables = {'Qbase':'Stream\nBaseflow','Qrech':'Stream\nLosses', 'flowing':'No. Days\nwith flow',
             # 'connected':'Connected', 'gaining':'Gaining',
            'Qout':'Streamflow'}
tests = ['Pearson','Spearman','Kendall']
tests = ['Pearson']
var_names = list(variables.keys())

# %%
# ref_out = pd.read_csv(join(proj_dir, 'num_sfr_coarse.csv'), index_col=0)
ref_out = pd.read_csv(join(proj_dir, 'coarse_reference.csv'), index_col=0)
ref_out=ref_out.rename(columns={'num_coarse':'num_sfr'})
# ref_out


# %%
lak_head_all = pd.read_csv(join(out_dir, 'lak_head_timeseries.csv')).set_index('totim')

# %%

# %%
wb_all = pd.DataFrame()

for r in np.arange(0,100):
    model_ws = join(all_model_ws, 'realization'+ str(r).zfill(3))
    # load summary water budget
    wb = pd.read_csv(model_ws+'/flow_budget.txt', delimiter=r'\s+')
    wb['GHB_NET'] = wb['GHB_IN'] - wb['GHB_OUT']
    wb['kstpkper'] = list(zip(wb.STP-1,wb.PER-1))
    wb_all = pd.concat((wb_all, wb.assign(realization=r)))

# %%
wb_avg = wb_all.groupby('realization').mean(numeric_only=True)

# %%
git_proj = join(git_dir, 'CosumnesRiverRecharge','Projects', 'Oneto_Denier')
wb_dict = pd.read_excel(join(git_proj, 'mf_wb_color_dict.xlsx'), sheet_name = 'owhm_wb_dict', comment='#')
sfr_dict = pd.read_excel(join(git_proj, 'mf_wb_color_dict.xlsx'), sheet_name = 'plt_ref', comment='#')
plt_names = pd.concat((sfr_dict, wb_dict)).set_index('flux')
# # pd.read_excel?

# %% [markdown]
# ## Correlations for all time
# I didn't know this when I first did the statistics but pandas has a function .corr() built in that allows you to compare across all columns the pearson, spearman or kendall correlation coefficient matrix. .corrwith() is pairwise correlations if desired in a specific pairing.

# %%
# the data is already summed across the year and should be averaged across reaches
corr_bool = sfr_yr_sum_all.groupby(['realization', 'segment', 'reach']).sum().groupby('realization').mean()
corr_bool = corr_bool[['flowing','connected','gaining']]
# seepage data just needs to be summed again 
corr_seep = sfr_yr_sum_all.groupby('realization').sum(numeric_only=True)[['Qbase','Qrech']]
# the baseflow needs a log transform represent properly
corr_seep.loc[corr_seep.Qbase==0, 'Qbase'] = 1 # convert 0 values to 1 for log transform
corr_seep.Qbase = np.log10(corr_seep.Qbase)

# flow data should be averaged
# i tried log10 transform but was problematic with zeros and replacement value, also had lower correlation
corr_flow = sfr_last_all[['realization','Qout']]
corr_flow = corr_flow.groupby('realization').mean(numeric_only=True)[['Qout']]

# join together the data for correlations
corr_all = corr_seep.join(corr_bool).join(corr_flow)

corr_all = corr_all[var_names]

# %%
# corr_matrix_in.hist() # all variables are relatively normally distributed

# %% [markdown]
# Showing that baseflow is correlated with groundwater elevations in the floodplain or stream helps break down further that it is groundwater levels that drive the baseflow.

# %%
# I should present the correlation matrix in the appendix
# helps show potential relationships
corr_matrix_in = corr_all.copy()
corr_matrix_in['num_sfr_coarse'] = ref_out.num_sfr.values
corr_matrix_in['num_lak_coarse'] = ref_out.num_lak.values
# corr_matrix_in['lak_gwe_mean'] = lak_head_all.mean(axis=0).values
# check for correlation with other water budget components
wb_cols = ['LAK_IN', 'ET_OUT','GHB_NET'] 
# wb_avg[wb_cols]
corr_matrix_in = pd.concat((corr_matrix_in, wb_avg[wb_cols]), axis=1)
corr_matrix_in.columns = plt_names.loc[corr_matrix_in.columns.values,'name']

# %%
corr = corr_matrix_in.corr('pearson')
# Fill diagonal and upper half with NaNs
mask = np.zeros_like(corr, dtype=bool)
mask[np.triu_indices_from(mask)] = True
corr[mask] = np.nan
corr.style.background_gradient(cmap='coolwarm', axis=None, vmin=-1, vmax=1).highlight_null(color='#f1f1f1').format(precision=3)

# %% [markdown]
# The full simulation correlations should be included to help explain cross-correlations and the potential influence of multiple inputs. The key takeaway is that we should consider the correlations due to the number of coarse reaches because of the consistent higher pearson's r across all areas. For water availability we could use streamflow, stream recharge or ET as an observation since they are all tightly correlated.
# - because streamflow and losses and ET are so correlated then it might not be needed to show pearson's r individually for stream flow and losses.

# %%
# check gaussian distribution
# v = 'Qout'
# fig,ax = plt.subplots(1,2, figsize=(9,3))
# (corr_all[v]+1).apply(np.log10).hist(ax=ax[0])
# (corr_all[v]+1).hist(ax=ax[1])


# %%
def calc_corr_stats(corr_all, coarse_ref):
    # take pearson's r
    pr = corr_all.apply(lambda x : pearsonr(coarse_ref, x), axis=0)
    pr.index=['r','p']
    pr['type'] = 'Pearson'
    sr = corr_all.apply(lambda x : spearmanr(coarse_ref, x), axis=0)
    sr.index=['r','p']
    sr['type'] = 'Spearman'
    kt = corr_all.apply(lambda x : kendalltau(coarse_ref, x), axis=0)
    kt.index=['r','p']
    kt['type'] = 'Kendall'
    # join data together
    corr_out = pd.concat((pr, sr, kt))
    return(corr_out)
corr_out = calc_corr_stats(corr_all, coarse_ref.num_coarse)
# corr_out = calc_corr_stats(corr_all, ref_out.num_sfr.values)
# corr_out = calc_corr_stats(corr_all, ref_out.num_lak.values)
# corr_out = corr_out[[var_names+['type']]]

# the applicability of certain variables is questionable since Qrech, flowing, Qout have median
# p values of > 0.05
corr_out.loc['p']


# %% [markdown]
# The correlations are all below 0.15 when looking at the number of coarse lake segments.
#

# %%
fig, ax = plt.subplots()
# prepare for plotting
corr_out = corr_out.loc['r'].set_index('type').transpose()

corr_out.plot(kind='bar', ax=ax, rot=25)
# fix column names
ax.set_xticks(ticks = np.arange(0,corr_all.shape[1]), labels = [variables[v] for v in corr_all.columns.values])
plt.ylabel('Correlation Coefficient')
# ax.xticks

# %% [markdown]
# ## Correlations grouped by water year
# Did I take these numbers from SASb GSP or from Rich?
# |WY |Types|WYsum|
# |:----|:---------------|-----|
# |2015 | Critical|9.23|
# |2016 | Below Normal|17.47|
# |2017 | Wet|37.82|
# |2018 | Below Normal|12.86|
#
# I recalculated these with the Fair Oaks station and found 15.2, 17.2, 37 and 23 inches of rainfall each Water Year

# %% [markdown]
# ## Correlations grouped by season and year
# Need to be careful with p-values because the plots with low correlations tend to have p-values greater than 0.05 and the plots with no correlation occur when the value is constant for all realizations (e.g., no baseflow or flow in any realization in WY2015 month 7).

# %% [markdown]
# Apply a log transform to baseflow to improve normality.

# %%
# yearly
# the data is already summed across the year and should be averaged across segments
corr_bool = sfr_3mon_all.copy().set_index('realization',append=True)[['flowing','connected','gaining']]
# sum across all segments for seepage
corr_seep = sfr_facies_all.groupby('realization').resample('3MS').sum(numeric_only=True)[['Qbase','Qrech']]
# flow data should be averaged
corr_flow = sfr_last_all.groupby(['realization']).resample('3MS').mean(numeric_only=True)[['Qout']]
# join together the data for correlations
corr_all_in = corr_seep.join(corr_bool).join(corr_flow).reset_index('dt')

# %%
corr_all_in['month'] = corr_all_in.dt.dt.month
corr_all_in['wy'] = corr_all_in.dt.dt.year
corr_all_in.loc[corr_all_in.month==10,'wy'] +=1

corr_all_in = corr_all_in[var_names+['dt','month','wy']]


# %%
# goal is to present the data used to calculate the correlations in box plot format
# box plot can't seem to handle the xtick formatting and throws an error even though there should only be 4 xticks, it thinks tehre are 8
m_cols=[10,1,4,7]
# m_cols = corr_all.month.unique()
fig, ax = plt.subplots(len(var_names),4, figsize=(8,5.3), sharex='col', sharey='row', layout='constrained')
for nv, v in enumerate(variables.keys()): #['Qbase','Qout']
    corr_plt = corr_all_in[[v,'wy','month']] 
    ax[nv, 0].set_ylabel(variables[v])

    for nwy, wy in enumerate(corr_plt.wy.unique()):
        ax_n = ax[nv, nwy]
        corr_wy = corr_plt[corr_plt.wy==wy].drop(columns=['wy']).pivot(columns='month', values=v)
        corr_wy = corr_wy.reset_index(drop=True).rename_axis(None,axis=1).loc[:, m_cols]

        # # correct order of tests and months to fix sorting
        # corr_wy.plot(column=m_cols, kind='box', ax=ax_n, legend=False, rot=0)
        # corr_plt[corr_plt.wy==2015].plot(by='month', column=v, kind='box', ax=ax_n)
        ax_n.boxplot(corr_wy.values, labels=m_cols);
        # ax_n.set(xticklabels=m_cols)

for nwy, wy in enumerate(corr_plt.wy.unique()):
    ax_n = ax[0, nwy]
    ax_n.set_title(wy)
# ax[0,0].ticklabel_format(style='plain')

plt.savefig(join(fig_dir, 'Appendix', 'corr_boxplots_season_WY.png'), bbox_inches='tight')

# %%
corr_all = corr_all_in.copy()
# the baseflow needs a log transform represent properly
corr_all.loc[corr_all.Qbase==0, 'Qbase'] = 1 # convert 0 values to 1 for log transform
corr_all.Qbase = np.log10(corr_all.Qbase)




# %%
corr_out = pd.DataFrame()
for n_dt in corr_all.dt.unique():
# corr_all.groupby('WY').apply(lambda x : pearsonr(coarse_ref.num_coarse, x))
    corr_dt = corr_all[corr_all.dt==n_dt].drop(columns=['dt'])
    corr_out = pd.concat((corr_out, calc_corr_stats(corr_dt, coarse_ref.num_coarse).assign(dt=n_dt)))
    # corr_out = pd.concat((corr_out, calc_corr_stats(corr_dt, ref_out.num_lak).assign(dt=n_dt)))
corr_out['month'] = corr_out.dt.dt.month
corr_out['wy'] = corr_out.dt.dt.year
corr_out.loc[corr_out.month==10,'wy'] +=1
corr_out = corr_out.drop(columns=['dt'])   

# %%
# there are NA values prodcued in some cases during the winter, only for flowing variable
# the NAs occur because in 2016 every realization had 100% of days in winter,spring with flow so there is
# no relationship at all

# %%
# the applicability of certain variables is questionable since Qrech, flowing, Qout have median
# p values of > 0.05
# corr_out.loc['p']
# corr_out.loc['p',variables.keys()][corr_out.loc['p', variables.keys()]>0.05]
# corr_out
pval_large = corr_out.loc['p'][(corr_out.loc['p', var_names]>0.05).any(axis=1).values]
pval_large[pval_large[var_names]<0.05] = np.nan
pval_large


# %% [markdown]
# Helen asked why I used more than one correlation coefficient and I didn't have a clear response other than that Steve did it so I'm simplifying to the Pearson's. I also should use color to highlight significance and draw lines at 0.3 and 0.5 to mark the moderate or strong correlations.
# - I could use: {\*:<0.05, \**:<0.01, \***:<0.001}

# %%
def plt_corr_lines(ax_n, ms):
    ax_n.axhline(0, linewidth=0.5, color='black')
    ax_n.axhline(0.3*ms, color='black',linewidth=0.5, linestyle='--')
    ax_n.axhline(0.5*ms, color='black',linewidth=0.5)


# %%
# v='Qbase'
# corr_plt = corr_out.loc['r'][[v,'type','wy','month']] 
# p_plt = corr_out.loc['p'][[v,'type','wy','month']]
# corr_plt[~(p_plt[v]<.05).values]
# ax_n.get_xticks()

# %%
# p_lab = p_wy[tests].loc[corr_out.month.unique()].copy()
# plt.scatter(np.arange(0,4)[p_lab.Pearson<0.05], corr_wy[tests].loc[corr_out.month.unique()][p_lab.Pearson<0.05],marker='*', color='black')
# plt.scatter(np.arange(0,4)[p_lab.Pearson<0.01], corr_wy[tests].loc[corr_out.month.unique()][p_lab.Pearson<0.01]+0.005,marker='*', color='black')
# ax[0,0].get_xticks()

# %%
fig, ax = plt.subplots(len(var_names),4, figsize=(8,5.3), sharex=True, sharey='row', layout='constrained')
for nv, v in enumerate(variables.keys()): #['Qbase','Qout']
    corr_plt = corr_out.loc['r'][[v,'type','wy','month']] 
    p_plt = corr_out.loc['p'][[v,'type','wy','month']] 
    for nwy, wy in enumerate(corr_plt.wy.unique()):
        ax_n = ax[nv, nwy]
        corr_wy = corr_plt[corr_plt.wy==wy].drop(columns=['wy']).pivot_table(index='month', values=v, 
                                                                             columns='type', dropna=False)        
        p_wy = p_plt[p_plt.wy==wy].drop(columns=['wy']).pivot_table(index='month', values=v, 
                                                                             columns='type', dropna=False)
        # correct order of tests and months to fix sorting
        corr_wy[tests].loc[corr_out.month.unique()].plot(kind='bar', ax=ax_n, legend=False, rot=0)
        # identify signficant correlations
        p_lab = p_wy[tests].loc[corr_out.month.unique()].copy()
        ax_n.scatter(np.arange(0,4)[p_lab.Pearson<0.05], corr_wy[tests].loc[corr_out.month.unique()][p_lab.Pearson<0.05],marker='*', color='black')
        # could further mark out 0.01 vs 0.001 but gets more complicated
        # ax_n.scatter(np.arange(0,4)[p_lab.Pearson<0.01], corr_wy[tests].loc[corr_out.month.unique()][p_lab.Pearson<0.01]-0.05,marker='*', color='black')
        # ax_n.scatter(np.arange(0,4)[p_lab.Pearson<0.001], corr_wy[tests].loc[corr_out.month.unique()][p_lab.Pearson<0.001]-0.1,marker='*', color='black')

        ax_n.set_ylabel(variables[v])
        if corr_wy[tests].min().values[0]<0:
            plt_corr_lines(ax_n, -1)
        if corr_wy[tests].max().values[0]>0:
            plt_corr_lines(ax_n, 1)

for nwy, wy in enumerate(corr_plt.wy.unique()):
    ax_n = ax[0, nwy]
    ax_n.set_title(wy)
#     ax_n.set_xticks(ticks = np.arange(0,corr_mon.shape[0]), labels = corr_out.month.unique())


# ax[0,1].legend(loc='best')
plt.savefig(join(fig_dir, 'corr_stats_season_WY.png'), bbox_inches='tight')

# %% [markdown]
# A lot of good discussion could come from the correlation by year and season figure. FOr example:
# - the flowing variable shows much stronger correlations in the dry season in wet and above normal years because the coarser facies are having a proportionally larger impact when there is less flow in the river. And the correlation is signficiant in all seasons in a critical year because there was such limited streamflow in the winter that the seepage due to the streambed had a larger effect. 
# - For the streambed recharge the correlation is equally strong in season in the wet and above normal year because there is sufficient flow available to make sure of the coarser facies for recharge. While the critical and below normal year see a less signficiant correlation in the summer months because there is not enough water in the stream for the number of coarse facies to make an impact since the stream goes dry.
# - the baseflow has a similar pattern as streambed recharge with a greater tendency to see weaker correlations when there is less water available because the groundwater elevations have already lowered such that there is no possibility of baseflow. This is most significant in the critical year when there is no correlation in spring because there isn't enough water for the coarse facies to recharge during storm events to make it available afterward for baseflow. This is where we start to see that the coarse facies provide baseflow benefits only if there is sufficiently large flow available. And we see that the correlations does not happen in the dry season.
#
# The correlations by year and season help identify that the coarse deposits consistently increase streambed leakage and decrease streamflow availability in all years and seasons, with the exception of summer 2015 when the channel went dry. The baseflow correlation shows that benefits do not occur in the dry season, and are limited to winter months in non-critical years with very limited benefits in the peak wet season for a critical year.

# %% [markdown]
# ### Questions on correlations
# How did you find these three correlative statistics?
# I noticed you adjusted for normallity for Pearson's r with the variable and predictor? Was that determined by looking at the distribution of the variable for the 100 sites? In my case I have 100 realizations
# How did you find the actual value of the statistic relevant? I'm considering primarily looking at the comparison between statistics to show where correlation is stronger?
#
# Pearson's r requires normality, Steve applied log10 transformation selectively on variables and predictors to improve normality. Spearman's rho and Kendall's tau. In this case the variable is the number of coarse stream segments, the histogram shows it sort of normally distributed with a slight peak around 30, but the data has lots of highs and lows probably because I'm plotting the number of segments with coarse and not necessarily the number of coarse bodies connecting to the stream.
#
# - Graham suggested to show the linearity of the data in the model appendix. These plots can be saved then imported all at once to a word document.

# %%
coarse_ref.quantile([0,.5,1])*100/len(grid_sfr)


# %%
bins = len(coarse_ref.num_coarse.unique())
coarse_ref.hist('num_coarse',
                # bins = bins
                bins= np.arange(0, int(coarse_ref.max().values[0])+4, 4) # group with bins of size 4 to account for downscaling from 200 m to 100m gri dcell
               )
plt.xlabel('Number of coarse reaches')
plt.ylabel('Number of realizations')
plt.title(None)
plt.savefig(join(fig_dir, 'Methods','coarse_segments_SFR_histogram.png'), dpi=300, bbox_inches='tight')

# %%
# exampel regression
plt_month = sfr_3mon_all[sfr_3mon_all.month==1]
df_plt = plt_month[plt_month.WY==2017]
# ax = df_plt.plot(x='num_coarse',y='Qbase', kind='scatter', color='gray')

# regr = linear_model.LinearRegression()
# regr.fit(df_plt[['num_coarse']].values, df_plt[['Qbase']].values)
# x_range = np.array([[df_plt.num_coarse.min()], [df_plt.num_coarse.max()]])
# ax.plot(x_range, regr.predict(x_range), color='black', linewidth=3)
# r2_val = r2_score(df_plt[['Qbase']], regr.predict(df_plt[['num_coarse']].values))
# ax.annotate('$r^2$: '+ str(np.round(r2_val,3)), (0.85,0.85), xycoords='axes fraction')

# %%
def corr_plt_seep(name, ylabel, log=False):
    """ name is Qbase or Qrech"""
    fig,ax = plt.subplots(4,4, figsize=(8,5.3), sharex=True, sharey=False, layout='constrained')
    # plot of number of connected days across realizations
    # plot across years in rows
    for yn, wy in enumerate(plt_month.WY.unique()):
        df_wy = sfr_3mon_all[sfr_3mon_all.WY==wy]
        # plot across 3 month blocks in columns
        for mn, mon in enumerate(sfr_3mon_all.month.unique()):
            df_plt = df_wy[df_wy.month==mon].copy()
            if log:
                df_plt.loc[df_plt[name]==0,name] = 1
                df_plt[name] = np.log10(df_plt[name])
            ax_n = ax[yn, mn]
            df_plt.plot(x='num_coarse',y=name, kind='scatter', color='gray', ax = ax_n)
            # linear, regression
#             sns.regplot(x='num_coarse', y=name,data=df_plt,ci=False,ax=ax_n, label='Linear Model', scatter=False)
            regr = linear_model.LinearRegression()
            regr.fit(df_plt[['num_coarse']].values, df_plt[[name]].values)
            x_range = np.array([[df_plt.num_coarse.min()], [df_plt.num_coarse.max()]])
            ax_n.plot(x_range, regr.predict(x_range), color='black', linewidth=3)
            r2_val = r2_score(df_plt[[name]], regr.predict(df_plt[['num_coarse']].values))
            ax_n.annotate('$R^2$: '+ str(np.round(r2_val,3)), (0.1,0.8), xycoords='axes fraction')
            # plot clean up
            ax_n.set_ylabel('')
            ax_n.set_xlabel('')

    for yn, wy in enumerate(plt_month.WY.unique().astype(int)):
        ax_n = ax[yn, 0]
        ax_n.set_ylabel(str(wy))
    labels_3mon = ['Oct-Dec','Jan-Mar','Apr-Jun','Jul-Sep']
    for mn, mon in enumerate(sfr_3mon_all.month.unique()):
        ax_n = ax[0, mn]
        ax_n.set_title(labels_3mon[mn])
    fig.text(-0.01, 0.4, 'Mean Seasonal '+ylabel, rotation=90)
    fig.text(0.4, -0.01, 'Number of Coarse Stream Reaches')

    fig.tight_layout()


# %%
labels=['Recharge ($m^3/day$)','Baseflow ($m^3/day$)','Days with flow','Streamflow ($m^3/day$)']
var = ['Qrech','Qbase', 'flowing', 'Qout']
# var,labels = ['Qbase'], ['Baseflow']
for n, param in enumerate(var):
    corr_plt_seep(param, labels[n],log=True)

    plt.savefig(join(fig_dir,'Appendix', 'corr_scatter_'+param+'.png'), bbox_inches='tight')

# %% [markdown]
# # Application figure
# Summary of net change in baseflow and ET for all realizations compared to the homogeneous case.
#
# And importantly, how does this change baseflow in seasons of concern, regarding the cosumnes functional flow regime. Between the homogeneous and hetergeneous realizations we need to sum across all segments because it's not fair to compare individual segments with the heterogeneity. But each day is a valuable comparison because a lack of flow is potential harmful to up-migrating salmon especially in fall.
#
# Although the number of days with flow doesn't have a significant relationship with heterogeneity I can still present the results holistically.
#
#

# %%
# sfr_facies_all
# sum for the seepage
total_flow = sfr_facies_all.groupby('realization').resample('3MS').sum(numeric_only=True)[['Qbase']]
total_flow = total_flow.reset_index('realization')

# homogeneous case
h_flow = h_sfrdf.resample('3MS').sum(numeric_only=True)[['Qbase']]
h_conn = h_sfrdf.resample('3MS').mean(numeric_only=True)[['flowing']]

# %%

# %%
nseg = sfrdf.segment.unique().shape[0]
# average for days of connection (facies are already summed so sum together then divide by number of segments)
# now it is the fraction of the segments with flow
total_conn = sfr_facies_all.groupby(['realization', 'dt']).sum(numeric_only=True)[['flowing']]/nseg
total_conn = total_conn.reset_index('realization')
# can calculate the fraction of a season with flow
total_conn = total_conn.groupby('realization').resample('3MS').mean().reset_index('realization', drop=True)

# %%

# %% [markdown]
# Similar to the plot showing days with flow by segment, the overall number of days with flow doesn't change much when looking at the histogram. From the earlier plots it looked like flow was extended at the downstream longer while reduced at the upstream. To show how downstream conditions change, we should represent the flow by requiring all segments to be flowing, take minimum of all segments which would essentially be the same as taking the most downstream segment. 

# %%
plt_date = '2015-10-1'
total_conn[total_conn.index==plt_date].hist('flowing')
plt.axvline(h_conn[h_conn.index==plt_date].flowing[0],color='red')

#plot(x='realization',y='flowing',kind='scatter')

# %%
total_flow = total_flow.join(h_flow,rsuffix='_h')
total_flow['Qbase_diff'] = total_flow.Qbase - total_flow.Qbase_h

# %%
total_flow_fall = total_flow[total_flow.index.month==10]
fig, ax = plt.subplots(2,2)
for n, y in enumerate(total_flow_fall.index.unique()):
    ax_n = ax[int(n/2), n%2]
    sns.histplot(total_flow_fall[total_flow_fall.index==y], x='Qbase_diff', ax=ax_n,bins=10)

# %% [markdown]
# ## Last segment comparison
# Show plot of how the number of days of flow averaged doesn't change much from the homogeneous case, but the days with flow at the outlet are more noticeable.
# - based on Graham's recommendation we will see variability above and below the homogeneous case with the heterogeneous realizations, but what is likely more interesting is  showing for an individual realization why the flows were different. We still should provide some variability numbers for managers to understand why heterogeneity in reality might shift from the modeling of a homogeneous case.

# %%
# groupby season and year
total_last = sfr_last_all.groupby('realization', group_keys=False).resample('3MS').mean(numeric_only=True)
h_last = h_sfrdf[h_sfrdf[grp_col]==h_sfrdf[grp_col].max()]
h_last = h_last.resample('3MS').mean(numeric_only=True)


# %%
# we want to give something easy to digest so like a yes or no, or the number of realizations each season
# that are better than the homogeneous case
# total_last['month'] = total_last.index.get_level_values('dt').month
# total_last.groupby(['realization','month']).mean()

# %%
# calculate difference in flow from base case
total_last['Qout_diff'] = (total_last[['Qout']] - h_last[['Qout']]).values
total_last['flowing_diff'] = (total_last[['flowing']] - h_last[['flowing']]).values
# get values for fall only
total_last_fall = total_last[total_last.index.month==10]


# %% [markdown]
# Does it make sense to give a simple response such as this many realization provided more water or the volume difference average.

# %%
for n, n_y in enumerate(total_last_fall.index.year.unique()):
    more_flow = (total_last_fall[total_last_fall.index.year==n_y].flowing_diff>=0)
    more_flow = (total_last_fall[total_last_fall.index.year==n_y].Qout_diff>=0)

    print(n_y+1, more_flow.sum())

# %%
plt_date = '2017-10-1'

total_last[total_last.index==plt_date].hist('flowing')
plt.axvline(h_last[h_last.index==plt_date].flowing[0],color='red')


# %%

# %%
