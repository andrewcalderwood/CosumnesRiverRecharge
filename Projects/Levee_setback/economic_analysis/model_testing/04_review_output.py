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

# %% [markdown]
# Goal: script to review more specific parts of the output as the testing phase

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import glob
import sys
import time
from importlib import reload
import h5py

import pandas as pd
import numpy as np

# standard geospatial python utilities
import shapely
import geopandas as gpd

import matplotlib.pyplot as plt

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# %%
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Economic/'
# model_ws = loadpth+'historical_simple_geology_reconnection'

# model_ws = loadpth+'rep_crop_soilbudget'

model_ws = join(loadpth,'input_write_2014_2020', 'rep_crop_soilbudget')


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
new_dir = dirname(os.getcwd())

add_path(new_dir)
import functions.Basic_soil_budget_monthly as swb


# %%

# %%

def read_crop_arr_h5(crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "r") as f:
        grp = f['array'] 
        # print(grp.keys())
        dset = grp[crop]
        arr = dset[:]
    return(arr)
        
    

# %%
# pc_all.shape

# %%
crop='Alfalfa'
# crop='Pasture'
crop='Corn'
# crop='Grape'
# crop='Misc Grain and Hay'
year=2016
# need separte hdf5 for each year because total is 300MB, group by crop in array
fn = join(model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
pc_all = read_crop_arr_h5(crop, fn)

# # applied water (GW and SW are separate)
fn = join(model_ws, 'field_SWB', "GW_applied_water_WY"+str(year)+".hdf5")
irr_gw = read_crop_arr_h5(crop, fn)

fn = join(model_ws, 'field_SWB', "SW_applied_water_WY"+str(year)+".hdf5")
irr_sw = read_crop_arr_h5(crop, fn)

# %%
# get total irrigation (sum of surface water and groundwater)
irr_tot = irr_gw+irr_sw
# get season total water use
irr_tot_sum = irr_tot.sum(axis=1)
# irrigation totals should be around 42" as in the cost study
(irr_tot_sum/0.3048)*12

# %%
# plt.plot(np.transpose(irr_tot));

# plt.plot(np.transpose(irr_sw));

# %%
# # applied water (GW and SW are separate)
fn = join(model_ws, 'field_SWB', "yield_WY"+str(year)+".hdf5")
Y_A = read_crop_arr_h5(crop, fn)

fn = join(model_ws, 'field_SWB', "profit_WY"+str(year)+".hdf5")
pi = -read_crop_arr_h5(crop, fn)

# %% [markdown]
# The profit reported by the individual site optimization is very high (\\$ 750) compared to the average soil data for distinct DTW profiles (\\$ 450)
# - but the printed out profit shows there is a step so this is just an error? or the irr. eff. overweights the profit since the printed values were before updating the irrigation to account for irrigation efficiency.
# - Yeah it should be the reverse since the printed values have lower profit for the field by field profit
# - also updated code to use minimum of rooting depth and soil depth to potentially avoid excess water availability
#
# Max total profit is 7 tons * 225 \\$/ton so \\$1,575 - \\$78.86 = \\$1496.14
# - with surface water \\$8\acre-inch at 42 in gives \\$336 water cost so
# - I had not finished updating the cost/profit table which might be why the profit was so high

# %%
# fig,ax = plt.subplots(1,2, figsize=(6,3))
# ax[0].plot(Y_A)
# ax[1].plot(pi)
# fig.tight_layout()

# %% [markdown]
# # Save output to table format

# %%
var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)


# %%
crop_in = pd.read_csv(join(model_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'))
print(crop_in.name.unique())
crop_in = crop_in[crop_in.name==pred_dict[crop]]
print(pred_dict[crop])

# %%
(crop_in.pod_bool==1).any()

# %%
dtw_df = pd.read_csv(join(model_ws, 'field_SWB', 'dtw_ft_WY'+str(year)+'.csv'), index_col=0)
# temporary for code testing
# dtw_df = dtw_df.loc[:, crop_in.parcel_id.values.astype(str)]

dtw_df_mean = dtw_df.mean().values
# temporary adjustment to account for representative testing
# dtw_df_mean = np.hstack([dtw_df_mean[::2]]*2)
# dtw_df_mean = np.hstack([dtw_df_mean]*2)
dtw_df_mean.shape

# %%
# summary output from hdf5 into csv
out_summary = pd.DataFrame(np.transpose((irr_gw.sum(axis=1), irr_sw.sum(axis=1), Y_A, pi, dtw_df_mean)), 
             columns=['irr_gw_m','irr_sw_m', 'Y_A','pi', 'dtw_mean_ft'])
out_summary[['irr_gw_in','irr_sw_in']] = out_summary[['irr_gw_m','irr_sw_m']]*12/0.3048
out_summary.to_csv(join(model_ws,'output_summary_'+crop+'.csv'))

# %%
# temporary fix
out_summary.loc[out_summary.irr_gw_in>50,:] = np.nan
# out_summary

# %%

# %%
# temporary fix
# out_summary = out_summary.iloc[:19,:]


plt_cols = ['irr_gw_in','irr_sw_in', 'Y_A', 'pi', 'dtw_mean_ft']
fig, ax  = plt.subplots(len(plt_cols), 1,dpi=300, sharex=True, figsize=(8,6.5))
for n, v in enumerate(plt_cols):
    ax_n = ax[n]
    out_summary[v].plot(ax=ax_n)
    ax_n.set_ylabel(v)
    ax_n.ticklabel_format(style='plain', useOffset=False)
    
fig.suptitle(crop)
ax[-1].set_xticks(out_summary.index.values[::2], out_summary.dtw_mean_ft[::2].round(0), rotation=90)
ax[-1].set_xlabel('Season Mean DTW (ft)')
fig.tight_layout()

# %% [markdown]
# # Comparison of output

# %%
summary_lin = pd.read_csv(join(loadpth+'rep_crop_soilbudget','output_summary_'+crop+'.csv'), index_col=0)
summary_lin['dtw_id'] = np.arange(0, len(summary_lin))

# simplifying to just data without POD
summary_lin = summary_lin.iloc[:19]
# temporary drop extremes
summary_lin = summary_lin[summary_lin.irr_gw_in<50] 

# %%
# summary_lin

# %%
out_summary = pd.read_csv(join(loadpth+'crop_soilbudget','output_summary_'+crop+'.csv'), index_col=0)
# pre sort by dtw so values are easier to view when plotting
out_summary = out_summary.sort_values('dtw_mean_ft')

out_summary['id'] = np.arange(0, len(out_summary))



# %%
# crop_in
# currently no POD

# %%
out_summary_dtw = out_summary[['id','dtw_mean_ft']].sort_values('dtw_mean_ft')
out_lin = pd.merge_asof(out_summary_dtw, summary_lin, on='dtw_mean_ft', direction='nearest')
# sort values for plotting and add back the interpoalted DTW
out_lin = out_lin.rename(columns={'dtw_mean_ft':'dtw_mean_ft_sim'}).merge( summary_lin[['dtw_id','dtw_mean_ft']])
out_lin = out_lin.sort_values('id').reset_index(drop=True)


# %%
# out_lin.loc[out_lin.irr_gw_in >50] = np.nan

# %%
# out_summary.plot(x='dtw_mean_ft',y='dtw_mean_ft')

# %%
plt_cols = ['irr_gw_in','irr_sw_in', 'Y_A', 'pi', 'dtw_mean_ft']
fig, ax  = plt.subplots(len(plt_cols), 1,dpi=300, sharex=True, figsize=(8,6.5))
x_col = 'id'
for n, v in enumerate(plt_cols):
    ax_n = ax[n]
    out_summary.plot(x=x_col, y=v, ax=ax_n, label='Field by field')
    out_lin.plot(x=x_col, y=v, ax=ax_n, label = 'Nearest linear DTW value')
    ax_n.set_ylabel(v)
    
fig.suptitle(crop)

fig.tight_layout()

# %%
