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
run_dir = 'F://WRDAPP/GWFlowModel'
loadpth = run_dir +'/Cosumnes/Regional/'

# model_nam = 'historical_simple_geology'
model_nam = 'historical_simple_geology_reconnection'
base_model_ws = join(loadpth, model_nam)
# model_nam = 'foothill_vani10'
# model_nam = 'strhc1_scale'
# model_nam = 'sfr_uzf'
model_nam = 'parallel_realizations/realization005'

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
grid_sfr['facies'] = 'Mud'
for p in vka_quants.index:
    facies = vka_quants.loc[p]
    grid_sfr.loc[(grid_sfr.vka< facies.vka_max)&(grid_sfr.vka>= facies.vka_min),'facies'] = facies.facies
    # add color for facies plots

# drop routing segments before calculating distances
gdf_sfr = grid_sfr[~grid_sfr.iseg.isin(drop_iseg)]
gdf_sfr = gdf_sfr.set_index(['iseg','ireach'])[['rchlen','strtop','strhc1', 'vka', 'facies', 'geometry']]
gdf_sfr['Total distance (m)'] = gdf_sfr['rchlen'].cumsum()
pd_sfr = pd.DataFrame(gdf_sfr.drop(columns=['geometry']))


# %%
# cleaned version of sfr reach data to save for EcoFIP reference
gdf_sfr_out = gdf_sfr.drop(columns=['strhc1','vka','facies']).rename(columns={'Total distance (m)':'dist_m'})
gdf_sfr_out.to_file(join(proj_dir, 'GIS','sfr_reach_reference.shp'))

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
grp_cols = ['segment','reach', 'realization', 'Total distance (m)']

# simplify columns for output 
keep_cols = np.append(grp_cols, ['dt','Qin','Qaquifer','Qout', 'Qaquifer_rate', 'width','depth'])


# %%
sfrdf_all = pd.DataFrame()
for r in best10.realization.values:
    r_ws = join(loadpth, 'parallel_realizations','realization'+str(r).zfill(3))
    sfrdf = clean_sfr_df(r_ws, dt_ref, pd_sfr)
    sfrdf_all = pd.concat((sfrdf_all, sfrdf.assign(realization=r)))

# %%
# drop routing segments
sfrdf_all = sfrdf_all[~sfrdf_all.segment.isin(drop_iseg)]
# calculate the effective rate of seepage
sfrdf_all['Qaquifer_rate'] = sfrdf_all.Qaquifer/(sfrdf_all.rchlen*sfrdf_all.width)

# %%
# about 220 MB
sfrdf_all.reset_index()[keep_cols].to_hdf(join(out_dir, 'sfrdf_all.hdf5'), 
                              key='monthly', complevel=4, data_columns = grp_cols)

# %%
sfr_mon_all = sfrdf_all.groupby(grp_cols).resample('MS').mean(numeric_only=True).drop(columns=grp_cols)
sfr_mon_all = sfr_mon_all.reset_index()
sfr_mon_all['month'] = sfr_mon_all.dt.dt.month

# %%
# produces file of 6 MB without data_columns
# file of 12 MB with data columns of grp_cols
sfr_mon_all[keep_cols].to_hdf(join(out_dir, 'sfrdf_mon_all.hdf5'), 
                              key='monthly', complevel=4, data_columns = grp_cols)

# %%
# sfr_mon_chk = pd.read_hdf(join(out_dir,'sfrdf_mon_all.hdf5'), key='monthly')

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
fig,ax = plt.subplots()
for r in best10.realization.values:
    sfr_plt = sfr_avg_all.loc[sfr_avg_all.realization==r]
    sfr_plt.plot(x='Total distance (m)',y='Qaquifer_rate', ax=ax, legend=False)
    
ax.set_yscale('log')
ax.set_ylim(1E-2, 10**(magnitude(sfr_avg_all.Qaquifer_rate.max())+1))
ax.set_ylabel('Stream loss rate (m/day)')
sfr_avg.plot(x='Total distance (m)',y='Qaquifer_rate', ax=ax, legend=False, color='black', linestyle='--')

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
sfr_cv.plot()
