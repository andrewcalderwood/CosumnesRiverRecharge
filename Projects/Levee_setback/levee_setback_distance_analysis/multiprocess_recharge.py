#%%
import os
from os.path import basename, dirname, join
import sys

from multiprocessing import Pool
import multiprocessing
import subprocess

import numpy as np
import numpy.ma as ma
import pandas as pd
from scipy.optimize import minimize_scalar

from time import time

import geopandas as gpd
####################################################################################################
#%% 
## Set up directory referencing
# Package data
usr_dir = os.getcwd()
while basename(usr_dir) != 'Users':
    usr_nam = basename(usr_dir)
    usr_dir = dirname(usr_dir)
usr_dir = join(usr_dir, usr_nam)
git_dir = join(usr_dir,'Documents/GitHub/CosumnesRiverRecharge')
gwfm_dir = join(usr_dir, 'Box/research_cosumnes/GWFlowModel')
print( gwfm_dir)

# set box directory for output figures and data
box_dir = gwfm_dir+'/Levee_setback/levee_setback_distance_analysis/'

# tprogs_id = '' # original tprogs with conditioning data in output tsim
# tprogs_id = '_no_conditioning'
tprogs_id = '_no_cond_c3d'

data_dir = box_dir+ tprogs_id+'/data_output/'
fig_dir = box_dir+tprogs_id+'/figures/'

chan_dir = box_dir+'channel_data/'
gis_dir = chan_dir+'GIS/'

fxn_dir = git_dir+'/python_utilities'
if fxn_dir not in sys.path:
    sys.path.append(fxn_dir)
# sys.path
# import muskingum_recharge as mr

from importlib import reload
# reload(mr)

nrow = 100
ncol = 230
# dem data for cropping above land surface
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_linear.tsv')

# import muskingum_recharge
# reload(muskingum_recharge)
# if I comment out geopandas in the muskingum_recharge script then no issue
from muskingum_recharge import min_Q, mannings, calc_depth_arr, xs_setback

print('Header done')
####################################################################################################
#%% 
## Input data ## 

## channel data##
setbacks = np.arange(0, 3400,200)
# original XS data
xs_all_cln = pd.read_csv(chan_dir+'Elevation_by_XS_number_meters.csv', index_col='dist_from_center_m')
# smoothed XS data used for setback analysis
# xs_all_cln = pd.read_csv(chan_dir+'xs_levee_smooth.csv', index_col='dist_from_center_m')
num_segs = xs_all_cln.shape[1]
# wse_grid.to_file(gis_dir+'wse_grid.shp')

import h5py
# load array identifying row,col to setback ID (1,17)
f = h5py.File(join(chan_dir, 'setback_locs.hdf5'), "r")
# upper = 3, middle = 2, lower=1
local_str_setbacks = f['setbacks']['local'][:]
str_setbacks = f['setbacks']['regional'][:]
f.close()

# load array identifying row,col to XS id (1,28)
xs_arr = np.loadtxt(chan_dir+'XS_num_grid_reference.tsv')

# load flood typology characteristics (based on daily data 1908 - 2014) - median values 
#"cms_pk" for peak discharge, "pk_loc" for time to peak, and "log_no_d" for duration
flood_type = pd.read_csv(join(box_dir, 'whipple_grp6_w97ftmedians.csv'),index_col='Group.1')

# tprogs near surface data
soil_thick=2
# fn = chan_dir+'/tprogs_geomK_'+str(soil_thick)+'m_depth.tsv' # old version with linear dem and old function for sampling tprogs
fn = chan_dir+'/tprogs_geomK_'+str(soil_thick)+'m_depth_dem_mean.tsv'
# units of m/day
soil_K_out = np.loadtxt(fn, delimiter='\t')
soil_K = np.reshape(soil_K_out, (100, nrow, ncol))
# convert soil conductivity from m/d to m/s and apply vertical anisotropy factor
vani = 100
soil_K = (soil_K/vani)/86400 

# high flow array maps
flow_percentile = 6 #95
hf_tot_in =  np.loadtxt(data_dir+'surface_highflow_by_realization_'+str(flow_percentile)+'.tsv',delimiter = '\t')
hf_tot = np.reshape(hf_tot_in, (100, nrow, ncol))

# zonal statistics for elevation
zs = gpd.read_file(gwfm_dir+'/DIS_data/grid_elevation_m_statistics.shp')
# columns with different quantiles 0 to 100% of elevation
q_cols = zs.columns[zs.columns.str.contains('perc')]
df_elevs = zs[q_cols]

# convert quantile dataframe to a 3D array
arr_elev = np.zeros((df_elevs.shape[1], zs.row.max(),zs.column.max()))
for n in np.arange(0,df_elevs.shape[1]):
    arr_elev[n, zs.row-1, zs.column-1] = df_elevs.iloc[:,n]

# the minimum elevation for each subsegment applied outward across the transects
xs_mins_arr = np.loadtxt(chan_dir+'subsegments_xs_mins.tsv', delimiter='\t')
# need to correct segment definition to where the xs_mins subsegment data is
xs_arr[np.isnan(xs_mins_arr)] = np.nan

print('Input data loaded')
####################################################################################################
#%% 
## Pre-process ##

# rating curves for each segment and setback
xs_flow_all = pd.read_csv(join(chan_dir,'all_xs_50pt_rating_curves.csv'))

def depth_match(seg_flow, flow):
    """ Given a XS (nseg, setback) return the expected depth (m) given a flow (cms)"""
    # find flows above and below the input flow
    flow_diff = (seg_flow.flow_cms-flow)
    f_high = flow_diff[flow_diff>0].argsort().index[0]
    f_low = flow_diff[flow_diff<0].argsort().index[-1]
    match_d = seg_flow.loc[[f_low, f_high]].sort_values('flow_cms')
    # linearly interpolate to calculate exact depth
    flow_slope = (match_d.iloc[1].flow_cms-match_d.iloc[0].flow_cms)/(match_d.iloc[1].depth_m-match_d.iloc[0].depth_m)
    out_depth = match_d.iloc[0].depth_m + (flow-match_d.iloc[0].flow_cms)/flow_slope
    return(out_depth)
####################################################################################################
#%% Recharge analysis ##
def arr_to_h5(Q, rch_hf_arr, d_arr, d_xs, cell_frac, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    f = h5py.File(h5_fn, "w")
    grp = f.require_group('array') # makes sure group exists
    grp.attrs['units'] = 'cubic meters/second'
    grp.attrs['description'] = 'Each layer of the array is a day in the event'
    dset = grp.require_dataset('flow', Q.shape, dtype='f', compression="gzip", compression_opts=4)
    dset[:] = Q
    dset = grp.require_dataset('rch_hf', rch_hf_arr.shape, dtype='f', compression="gzip", compression_opts=4)
    dset[:] = rch_hf_arr
    dset = grp.require_dataset('depth', d_arr.shape, dtype='f', compression="gzip", compression_opts=4)
    dset[:] = d_arr
    dset = grp.require_dataset('XS_depth', d_xs.shape, dtype='f', compression="gzip", compression_opts=4)
    dset[:] = d_xs
    dset = grp.require_dataset('cell_frac', cell_frac.shape, dtype='f', compression="gzip", compression_opts=4)
    dset[:] = cell_frac
    f.close()
    

def realization_recharge(t, str_setbacks, region, ft):
    # typical winter baseflow, peak flow, peak location, total time (days)
    # flow of 23 m3/s listed by Whipple as floodplain cutoff
    q_base = 23 # 200*(0.3048**3)
    q_peak = flood_type.loc[ft,'cms_pk']
    # total duration in days 
    T = int(10**flood_type.loc[ft,'log_no_d'])
    p_l = flood_type.loc[ft,'pk_loc']
    tp = int(p_l*T)
    q_rise = np.linspace(q_base, q_peak, tp)
    q_fall = np.linspace(q_peak, q_base, (T-tp+1))
    q_in = np.append(q_rise, q_fall[1:])

    # allocate arrays - num flow steps, num setbacks, num segments
    Q = np.zeros((q_in.shape[0], len(setbacks), xs_all_cln.shape[1]+1)) # discharge for each XS
    d_xs = np.zeros((q_in.shape[0], len(setbacks), xs_all_cln.shape[1]+1)) # depth for each XS

    # set inflow for segment 1 across all setbacks and for all times
    # rate of cubic meters per second
    Q[:,:,0] = np.repeat(q_in.reshape(-1,1), len(setbacks), axis=1)
    # save depth arrays for each setbacks
    d_arr = np.zeros((q_in.shape[0], len(setbacks), nrow, ncol))
    wse_arr = np.zeros((len(setbacks), nrow, ncol))
    cell_frac = np.zeros((q_in.shape[0], len(setbacks), nrow, ncol))
    # save high recharge flows
    rch_hf_arr = np.zeros((q_in.shape[0], len(setbacks), nrow, ncol))

    tic = time()
    n = 0.048 # assume constant roughness for now

    # iterate across streamflows
    for qn in np.arange(0, q_in.shape[0]):
        # iterate across all cross-sections
        for nseg in np.arange(0,xs_all_cln.shape[1]):
            # iterate across all setbacks
            for s,setback in enumerate(setbacks):
                # for a given setback imagine there is an impenetrable levee blocking overbank flow
    #             xs_elevs = xs_all_cln.iloc[:,nseg][3100-setback:3300+setback]
                xs_elevs = xs_setback(xs_all_cln.iloc[:,nseg], setback)
                # boolean of row,col cells that fall within the segment and setback
                # fp_zon = (xs_arr==nseg)&(str_setbacks <= s+1)
                fp_zon = (xs_arr==nseg)&(str_setbacks[s]==1)
                # solve for depth that matches given flow, assume less than 1 cms is too small to calculate
                if Q[qn, s,nseg] >1:
                    # new depth solver uses pre-calculated rating curves
                    seg_flow = xs_flow_all[(xs_flow_all.nseg==nseg)&(xs_flow_all.setback==setback)]
                    depth = depth_match(seg_flow, flow=Q[qn, s, nseg])
                else:
                    depth = 0
                # save depth to an array for export
                d_xs[qn, s, nseg] = depth
                # join depth calculated at cross-section to corresponding model cells and corresponding setback
                # add elevation to minimum to apply segment midpoint as elevation rather than lowest point
                wse_arr[s, fp_zon] = depth + xs_mins_arr[fp_zon]
                # calculate depth of water at different elevation percentiles for segment and setback
                # to account for flood connectivity and avoid negligible, require a minimum of 0.1 m 
                diff = wse_arr[s,:]*fp_zon - 0.1 - arr_elev
                # when depth is negative remove
                diff[diff<0] = 0 #np.NaN
                # only keep cells where water level is above lowest elevation
                diffmax = np.nanmax(diff, axis=0)
                # keep cells where diffmax >0 
                x,y = np.where(diffmax>0)
                # find the highest elevation above which there is water, subtract 1 uses lower percentile
                # where the wse was below the minimum ignore
                bot_q = np.argmin(diff, axis=0)[x,y] - 1
                bot_q[bot_q<0] = 0
                top_q = np.argmin(diff, axis=0)[x,y] 
                # find percentage of interim quantile
                perc_q = (wse_arr[s, x,y] - arr_elev[bot_q, x,y])/(arr_elev[top_q, x,y] +1E-3 - arr_elev[bot_q, x,y])
                # need to account for when top_q == bot_q
                perc_q = np.where(arr_elev[bot_q, x,y]==arr_elev[top_q, x,y], 0, perc_q)
                # adjust for when wse > top_q
                perc_q = np.where(wse_arr[s, x,y]>arr_elev[top_q, x,y], 0, perc_q)
                # percent of cell area covered by flood
                cell_frac[qn, s,x,y] = (bot_q + perc_q)/10
                # depth for each cell is difference between water surface and average flooded ground elevation
                d_arr[qn, s,fp_zon] = np.nanmean(diff, axis=0)[fp_zon] #depth 
                # identify wse above surface elevation 
#                 d_arr[qn,:] = d_arr[qn,:] * cell_frac[qn] # caused overflow error
            # calculate vertical seepage with Darcy's equation assuming a saturated zone thickness similar to the lake bed in modflow
            # hydraulic conductivity is in m/s, hydraulic gradient is unitless, area is 200x200 m^2
            # q_seep = (soil_K[t,:,:])*hf_tot[t,:,:]*(200*200)*((d_arr[qn,:] + soil_thick)/soil_thick)
            q_seep = (soil_K[t,:,:])*(200*200)*((d_arr[qn,:] + soil_thick)/soil_thick)
            rch_hf_arr[qn,:,:,:] += (xs_arr==nseg) * q_seep * cell_frac[qn]
            # identify when the flow is less than the recharge predicted and recharge > 0 
            dry = (Q[qn, :, nseg] < np.nansum(rch_hf_arr[qn,:, xs_arr==nseg], axis=0)) & (np.nansum(rch_hf_arr[qn,:, xs_arr==nseg], axis=0)>0)
            if any(dry):
                # where the cells will end dry, reduce recharge so it sums to the flow into segment
    #             scale = (Q[qn, :, nseg][dry]/np.nansum(rch_hf_arr[qn,:, xs_arr==nseg], axis=(0))[dry])
                scale = (Q[qn, :, nseg]/(np.nansum(rch_hf_arr[qn,:, xs_arr==nseg], axis=(0)))) 
                scale = np.where(scale > 1, 1,scale) # only scale when flow is less than recharge
#             if any(scale<1):
#                 rch_hf_arr[qn, :,xs_arr==nseg][:,dry] *= scale
                rch_hf_arr[qn, :,xs_arr==nseg] *= scale
                print('scale', np.round(np.median(scale),2))
            Q[qn, :, nseg+1] = Q[qn, :, nseg] - np.nansum(rch_hf_arr[qn,:, xs_arr==nseg], axis=(0))
            if any(Q[qn,:,nseg+1]<0):
                print(qn, nseg+1)
                
                
    # base_fn = join(data_dir, 'type'+str(ft), 'r'+str(t).zfill(3)+'_')
    base_fn = join(data_dir, region, 'type'+str(ft), 'r'+str(t).zfill(3)+'_')
    arr_to_h5(Q, rch_hf_arr, d_arr, d_xs, cell_frac, base_fn+'output.hdf5')

    toc = time()
    print((toc-tic)/3600)
    return(Q, rch_hf_arr, d_arr, cell_frac)

###############################################################################
#%% Make short code to loop over local zones

# choose one function to use, regional or local
# def run_rech(t):
#     for zone in [1,2,3]:
#         for ft in [1,2,3]:
#             # 1, 2, 3 are floods long enough to apply to analysis
#             base_fn = join(data_dir, 'local_'+str(zone), 'type'+str(ft))
#             os.makedirs(base_fn, exist_ok=True)
#             realization_recharge(t, np.where(local_str_setbacks==zone, 1, 0), 'local_'+str(zone), ft)

def run_rech(t):
    region = 'regional'
    for ft in [1,2,3]:
        # 1, 2, 3 are floods long enough to apply to analysis
        realization_recharge(t, str_setbacks, 'regional', ft)

###############################################################################
#%% Multiprocess

# realization_recharge(1)

def main():
    # pool = Pool(processes=multiprocessing.cpu_count()-2)  # set the processes max number to the number of cpus
    pool = Pool(processes=16)  # with 100/25 that is 4 sets, 16 is 6.25 sets
    # result = pool.map(realization_recharge, range(100)) # original without adding inputs
    result = pool.map(run_rech, range(100))
    pool.close()
    pool.join()
    print(result)
    print('end')

if __name__ == "__main__":
    tic = time()
    main()
    toc = time()
    print('Total time: %.2f minutes' %((toc-tic)/60))