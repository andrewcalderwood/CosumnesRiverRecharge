# ---
# jupyter:
#   jupytext:
#     formats: py:percent
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
proj_dir = join(gwfm_dir, 'Projects','EcoFIP')
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
# method to set default plot parameters, no longer need to specify 300 dpi each time, may not need to specify dimensions either
# sns.set_theme(rc={"figure.dpi": 300})
plt.rcParams.update({"figure.dpi": 300})


# %%
run_dir = 'C://WRDAPP/GWFlowModel'
# run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'
base_model_ws = join(loadpth, model_nam)
# model_nam = 'foothill_vani10'
model_nam = 'strhc1_scale'
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


# %%
# cleaned version of sfr reach data to save for EcoFIP reference
gdf_sfr_out = gdf_sfr.drop(columns=['strhc1','vka','facies']).rename(columns={'Total distance (m)':'dist_m'})
# gdf_sfr_out.to_file(join(proj_dir, 'GIS','sfr_reach_reference.shp'))

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

# %% [markdown]
# # Output across realizations
# It's pretty slow (minutes) to load the 10 realizations and post-process (10s of seconds) them. It might make sense to process indiviudal

# %%
gdf_sfr = gpd.read_file(join(proj_dir, 'GIS','sfr_reach_reference.shp'))
gdf_sfr = gdf_sfr.rename(columns={'iseg':'segment','ireach':'reach','dist_m':'Total distance (m)'})


# %%
df_sfr = pd.DataFrame(gdf_sfr.drop(columns=['geometry']))


# %%
grp_cols = ['segment','reach', 'realization', 'Total distance (m)','rch_order']

# simplify columns for output 
keep_cols = np.append(grp_cols, ['dt','Qin','Qaquifer','Qout', 'Qaquifer_rate', 'width','depth'])


# %%
rewrite=False

# %%
if rewrite:
    sfrdf_all = pd.DataFrame()
    for r in best10.realization.values:
        r_ws = join(loadpth, 'parallel_realizations','realization'+str(r).zfill(3))
        sfrdf = clean_sfr_df(r_ws, dt_ref, pd_sfr)
        sfrdf_all = pd.concat((sfrdf_all, sfrdf.assign(realization=r)))
    
    # drop routing segments
    sfrdf_all = sfrdf_all[~sfrdf_all.segment.isin(drop_iseg)]
    # calculate the effective rate of seepage
    sfrdf_all['Qaquifer_rate'] = sfrdf_all.Qaquifer/(sfrdf_all.rchlen*sfrdf_all.width)
    # about 220 MB
    sfrdf_all.reset_index()[keep_cols].to_hdf(join(out_dir, 'sfrdf_all.hdf5'),  format='table',
                                  key='all', complevel=4, data_columns = grp_cols, mode='w')

# %%
sfrdf_all = pd.read_hdf(join(out_dir, 'sfrdf_all.hdf5'),  key='all', complevel=4)
# add total distance which was missing 
sfrdf_all = sfrdf_all.merge(df_sfr)


# %%

if rewrite:
    sfr_mon_all = sfrdf_all.groupby(grp_cols).resample('MS').mean(numeric_only=True).drop(columns=grp_cols)
    sfr_mon_all = sfr_mon_all.reset_index()
    sfr_mon_all['month'] = sfr_mon_all.dt.dt.month
    
    # produces file of 6 MB without data_columns
    # file of 12 MB with data columns of grp_cols
    sfr_mon_all[keep_cols].to_hdf(join(out_dir, 'sfrdf_mon_all.hdf5'),  format='table',
                                  key='monthly', complevel=4, data_columns = grp_cols, mode='w')

# %%
sfr_mon_all = pd.read_hdf(join(out_dir,'sfrdf_mon_all.hdf5'), key='monthly')
# add reference columns
sfr_mon_all = sfr_mon_all.merge(df_sfr)
sfr_mon_all['month'] = sfr_mon_all.dt.dt.month

# %%
# get average across time to plot spatial view
sfr_avg_all = sfr_mon_all.groupby(grp_cols).mean(numeric_only=True).reset_index()

# %%
sfr_avg = sfr_mon_all.groupby(['segment','reach']).mean(numeric_only=True).reset_index()
sfr_std = sfr_mon_all.groupby(['segment','reach']).std(numeric_only=True).reset_index()
sfr_cv = sfr_std.Qaquifer_rate/sfr_avg.Qaquifer_rate

# %%
from report_cln import magnitude

# %%
# plt.rcdefaults() # reset

# %%
fig,ax = plt.subplots(figsize=(6.5,5),dpi=300)
for r in best10.realization.values:
    sfr_plt = sfr_avg_all.loc[sfr_avg_all.realization==r]
    sfr_plt.plot(x='Total distance (m)',y='Qaquifer_rate', ax=ax, label=r,legend=False)

    
ax.set_yscale('log')
ax.set_ylim(1E-2, 10**(magnitude(sfr_avg_all.Qaquifer_rate.max())+1))
ax.set_ylabel('Stream loss rate (m/day)')

sfr_avg.plot(x='Total distance (m)',y='Qaquifer_rate', ax=ax, legend=False, color='black', linestyle='--',label='Mean')
plt.legend(title='Realization')


# %%
g = sns.lineplot(sfr_avg_all, x='Total distance (m)', y='Qaquifer_rate', errorbar=('ci',95))
g.set(yscale='log', ylim=(1E-2, 1E2), title = '10 Realizations')
g.set( ylabel='Qaquifer (m/day)')


# %% [markdown]
# - 230 column/5 is 46 so the first 46 rows don't have a scaling
# - the shift in CV is around 160 which is still at column 100
# - the first region with no 1/5 strhc scaling starts at reach 220
# - ultimately there is an impact because that region has more high values
# - throughout the rest of the domain there is lower seepage in the upper reaches because of the way the geologic model was sliced which shows more continuous coarse facies outcropping in the lower lying floodplain region
#     - this could also be driven by the conditioning data over-applying fines?
#     
# **It would take 4 hours to run the models again with a constant strhc scalar**

# %%
print('Max column after reach 160 was %.i' %grid_sfr.iloc[160:].column.max())
print('First reach before column 46 was %.i ' %grid_sfr[grid_sfr.column<46].index.min())

# %%
# sfr_cv.plot()

# %%
# plot histogram of leakage rate
r = 5
sfr_mon = sfr_mon_all[sfr_mon_all.realization==r].copy()

# %%
r=200
# test the monthly at a reach to look at relationship of flow and seepage
mon_chk = sfr_mon[sfr_mon['Total distance (m)']==df_sfr['Total distance (m)'].iloc[r]]
fig,ax = plt.subplots( layout='constrained')
# mon_chk.boxplot(by='month',column='Qaquifer', ax=ax[0])
mon_chk.boxplot(by='month',column='Qaquifer_rate', ax=ax)

# %% [markdown]
# Define limits for data by reach to help decide how much range there is. We can use quartiles and whiskers or confidence intervals.

# %%
# check normality if needed
# from scipy.stats import normaltest
# normaltest(mon_chk.Qaquifer_rate)

# %%
# import statsmodels.stats.api as sms
#     for r in np.arange(0, len(df_sfr)):
#         mon_chk = sfr_mon[sfr_mon['Total distance (m)']==df_sfr['Total distance (m)'].iloc[r]].copy()
#         # get the confidence intervals 
#         ci_l, ci_h = sms.DescrStatsW(mon_chk.Qaquifer_rate).tconfint_mean()

# %%
# fig,ax = plt.subplots()
# sfr_mon.boxplot(by='Total distance (m)', column='Qaquifer_rate',ax=ax)
# ax.set_yscale('log')
# ax.set_ylim(1E-2, None)
# ax.set_xticklabels(ax.get_xticklabels(), rotation=90);
# boxplots are messy 

# seaborn is slow to plot but works well
# g = sns.lineplot(sfr_mon_all, x='Total distance (m)', y='Qaquifer_rate',
#                  hue='realization', errorbar=('ci',95))
r = 5
# r = best10.realization.iloc[1]
sfr_mon = sfr_mon_all[sfr_mon_all.realization==r].copy()

g = sns.lineplot(sfr_mon, x='Total distance (m)', y='Qaquifer_rate', errorbar=('ci',95))
# g = sns.lineplot(sfr_mon, x='Total distance (m)', y='Qaquifer_rate', err_style="bars", errorbar=("se", 2))

g.set(yscale='log', ylim=(1E-2, None), title = 'Realization '+str(r))
g.set(xlim = (13E3, 40E3), ylabel='Qaquifer (m/day)')

# sns.set(rc={'figure.figsize':(6.5,6.5)})
# plt.savefig(join(fig_dir, 'stream_loss_with_CI_r'+str(r)+'.png'), bbox_inches='tight')

# %%

# %%
for r in best10.realization.values:
    sfr_mon = sfr_mon_all[sfr_mon_all.realization==r].copy()
    
    g = sns.lineplot(sfr_mon, x='Total distance (m)', y='Qaquifer_rate', errorbar=('ci',95))
    # g = sns.lineplot(sfr_mon, x='Total distance (m)', y='Qaquifer_rate', err_style="bars", errorbar=("se", 2))
    
    g.set(yscale='log', ylim=(1E-2, None), title = 'Realization '+str(r))
    g.set(xlim = (13E3, 40E3), ylabel='Qaquifer (m/day)')
    # sns.set(rc={'figure.figsize':(6.5,6.5)})
    plt.savefig(join(fig_dir, 'stream_loss_with_CI_r'+str(r)+'.png'), bbox_inches='tight')
    plt.close()

# %% [markdown]
# # Linear regression
# Perform linear regression between leakage rate and depth to identify regions where a scale leakage rate might be appropriate.

# %%
# statistics functions
from sklearn.metrics import r2_score, mean_squared_error
from sklearn import datasets, linear_model
# pands has this built in with the DataFrame.corr()
# from scipy.stats import pearsonr, spearmanr, kendalltau 

# %%
from hyd_utility import mark_outlier

# %%
mon_chk_plt = mark_outlier(mon_chk.Qaquifer_rate)
mon_chk_plt = mon_chk.loc[~mon_chk_plt.flier]
x_var = 'depth'
y_var='Qaquifer_rate'


# %%
# ax_n.set_xscale('log')
def plt_reg(df_plt, x_var, y_var, logx=False, logy=False):
    # linear, regression

    # perform the regression
    regr = linear_model.LinearRegression()
    regr.fit(df_plt[[x_var]].values, df_plt[[y_var]].values)
    x_range = np.array([[df_plt.depth.min()], [df_plt.depth.max()]])
    ax_n.plot(x_range, regr.predict(x_range), color='black', linewidth=1)
    r2_val = r2_score(df_plt[[y_var]], regr.predict(df_plt[[x_var]].values))
    ax_n.annotate('$R^2$: '+ str(np.round(r2_val,3)), (0.1,0.8), xycoords='axes fraction')
    # return the coefficient, intercept, and R2 fit
    return([regr.coef_[0][0], regr.intercept_[0], r2_val])

r=250
mon_chk = sfr_mon[sfr_mon['Total distance (m)']==df_sfr['Total distance (m)'].iloc[r]]
mon_chk_plt = mark_outlier(mon_chk.Qaquifer_rate)
mon_chk_plt = mon_chk.loc[~mon_chk_plt.flier]

fig,ax_n = plt.subplots()
mon_chk_plt.plot(x=x_var,y=y_var, kind='scatter',ax=ax_n)
m, x0, r2_val = plt_reg(mon_chk_plt.copy(), x_var, y_var)
plt.title('Reach '+str(r))


# %%
def get_lin_reg(sfr_mon, df_sfr):
    lin_fit = pd.DataFrame(np.zeros((len(df_sfr), 3)), columns=['coef','intercept','R2'])
    for r in np.arange(0, len(df_sfr)):
        mon_chk = sfr_mon[sfr_mon['Total distance (m)']==df_sfr['Total distance (m)'].iloc[r]]
        mon_chk_plt = mark_outlier(mon_chk.Qaquifer_rate)
        mon_chk_plt = mon_chk.loc[~mon_chk_plt.flier]
        m, x0, r2_val = plt_reg(mon_chk_plt, x_var, y_var)
        lin_fit.loc[r,:] = (m, x0, r2_val)
    return(lin_fit)

lin_fit = get_lin_reg(sfr_mon.copy(), df_sfr)
lin_fit.plot(y='R2')


# %%
lin_fit_all = pd.DataFrame()
for r in best10.realization.values:
    sfr_mon = sfr_mon_all[sfr_mon_all.realization==r].copy()
    lin_fit = get_lin_reg(sfr_mon.copy(), df_sfr)
    lin_fit_all = pd.concat((lin_fit_all, lin_fit.assign(realization=r)))
# make rch_order one based
lin_fit_all.index+=1
lin_fit_all.index.name='rch_order'


# %%
lin_fit_all.to_csv(join(out_dir, 'linear_regression_alltime.csv'))


# %% [markdown]
# For the regions with mostly disconnection, the fit is good. Bad fit in the foothills but we don't care as much about that region and bad fit in the lower most reaches which again we are less concerned about.
# - to further refine this, the seepage regression could also be done by month because in the summer much more of the system including the lower river is disconnected.

# %%
fig,ax = plt.subplots(2,1, sharex=True)
for r in best10.realization.values:
    lin_fit = lin_fit_all[lin_fit_all.realization==r].copy()
    lin_fit.plot(y='R2',ax=ax[0],label=r, legend=False)
    lin_fit.plot(y='coef',ax=ax[1],label=r, legend=False)

# plt.ylim(-1,1)
ax[0].legend( loc=(1.01,-0.4), title='Realization')
plt.xlabel('Reach')
ax[0].set_ylabel('$R^2$ score')
ax[1].set_ylabel('Lin. Reg. Coeff.')
ax[1].set_yscale('log')

# %%
# sfrdf_all.Qaquifer_rate

# %%
# need to use daily data maybe to do regression within monthss
# sfrdf = sfrdf_all[sfrdf_all.realization==r].copy()
# lin_fit = get_lin_reg(sfrdf.copy(), df_sfr)


# %%
