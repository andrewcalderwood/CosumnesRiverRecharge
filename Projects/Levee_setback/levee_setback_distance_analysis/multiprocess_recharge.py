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


####################################################################################################
#%% 
## Input data ## 

## channel data##
setbacks = np.arange(0, 3400,200)
# smoothed XS data used for setback analysis
xs_levee_smooth = pd.read_csv(chan_dir+'xs_levee_smooth.csv', index_col='dist_from_right_m')
num_segs = xs_levee_smooth.shape[1]
# wse_grid.to_file(gis_dir+'wse_grid.shp')

# load array identifying row,col to setback ID (1,17)
str_setbacks = np.loadtxt(chan_dir+ 'regional_str_setback_id_arr.tsv', delimiter='\t').astype(int)
# str_setbacks = ma.masked_where(str_setbacks==0,str_setbacks)
str_setbacks = np.where(str_setbacks==0,np.NaN, str_setbacks)

# load array identifying row,col to XS id (1,28)
xs_arr = np.loadtxt(chan_dir+'XS_num_grid_reference.tsv')

# load flood typology characteristics (based on daily data 1908 - 2014) - median values 
#"cms_pk" for peak discharge, "pk_loc" for time to peak, and "log_no_d" for duration
flood_type = pd.read_csv(join(box_dir, 'whipple_grp6_w97ftmedians.csv'),index_col='Group.1')


# tprogs near surface data
soil_thick=2
fn = chan_dir+'/tprogs_geomK_'+str(soil_thick)+'m_depth.tsv'
# units of m/day
soil_K_out = np.loadtxt(fn, delimiter='\t')
soil_K = np.reshape(soil_K_out, (100, nrow, ncol))
# convert soil conductivity from m/d to m/s and apply vertical anisotropy factor
vani = 100
soil_K = soil_K/86400/vani 

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
####################################################################################################
#%% 
## Pre-process ##
# find minimum from channel center
xs_mins = xs_levee_smooth.loc[3100:3300].min(axis=0)
xs_mins.index = xs_mins.index.astype(int)
slope = xs_mins.diff().rolling(2, center=True, closed='right').mean().bfill()/2000*-1
adj_xs_mins = np.append(xs_mins[0], (xs_mins[0]-slope.cumsum()*2000))

# for n in flood_type.index:
# 1, 2, 3 are floods long enough to apply to analysis
ft=3
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

# print(q_in)
####################################################################################################
#%% Recharge analysis ##



def realization_recharge(t):
    # allocate arrays - num flow steps, num setbacks, num segments
    Q = np.zeros((q_in.shape[0], len(setbacks), xs_levee_smooth.shape[1]+1))
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
        for nseg in np.arange(0,xs_levee_smooth.shape[1]):
            # iterate across all setbacks
            for s,setback in enumerate(setbacks):
                # for a given setback imagine there is an impenetrable levee blocking overbank flow
    #             xs_elevs = xs_levee_smooth.iloc[:,nseg][3100-setback:3300+setback]
                xs_elevs = xs_setback(xs_levee_smooth.iloc[:,nseg], setback)
                # boolean of row,col cells that fall within the segment and setback
                fp_zon = (xs_arr==nseg)&(str_setbacks <= s+1)
                # solve for depth that matches given flow, assume less than 1 cms is too small to calculate
                if Q[qn, s,nseg] >1:
                    # solve for depth that matches given flow
                    res = minimize_scalar(min_Q, args = (xs_elevs, n, slope.iloc[nseg], Q[qn, s,nseg]), 
                                              bounds=(1E-3, 10), method='Golden', tol=1E-3)
                    # if large flow error, likely stuck in local minimum try setting higher or lower
                    if mannings(res.x, xs_elevs, n, slope.iloc[nseg]) > Q[qn, s,nseg]*1.05:
                        res = minimize_scalar(min_Q, args = (xs_elevs, n, slope.iloc[nseg], Q[qn, s,nseg]), 
                                      bounds=(1E-3,res.x-0.1), method='bounded')
                    if mannings(res.x, xs_elevs, n, slope.iloc[nseg]) < Q[qn, s,nseg]*0.95:
                        res = minimize_scalar(min_Q, args = (xs_elevs, n, slope.iloc[nseg], Q[qn, s,nseg]), 
                                      bounds=(res.x+0.1,10), method='bounded')
                    if res.fun>0.05*Q[qn, s,nseg]: # greater than 5% difference try to fix with bounded solving
                        print(str(nseg),' ', s, '%.2f'%res.fun, 'iter %i'%res.nit, ', ',res.success, 'd %.2f'%res.x)
                    depth = np.copy(res.x)
                else:
                    depth = 0
                # join depth calculated at cross-section to corresponding model cells and corresponding setback
                # add elevation to minimum to apply segment midpoint as elevation rather than lowest point
#                 wse_arr[s, fp_zon] = depth + xs_elevs.min() + slope.iloc[nseg]*1000
                wse_arr[s, fp_zon] = depth + xs_mins_arr[fp_zon]
#                 if any(np.isnan(xs_mins_arr[fp_zon])):
#                     print('Seg',str(nseg),'Setback',str(s), xs_mins_arr[fp_zon])
                # calculate depth of water at different elevation percentiles for segment and setback
                diff = wse_arr[s,:]*fp_zon - arr_elev
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
                # percent of cell area covered by flood
                cell_frac[qn, s,x,y] = (bot_q + perc_q)/10
                # depth for each cell is difference between water surface and average flooded ground elevation
                d_arr[qn, s,fp_zon] = np.nanmean(diff, axis=0)[fp_zon] #depth 
                # identify wse above surface elevation 
#                 d_arr[qn,:] = d_arr[qn,:] * cell_frac[qn] # caused overflow error
            # calculate vertical seepage with Darcy's equation assuming a saturated zone thickness similar to the lake bed in modflow
            # hydraulic conductivity is in m/s, hydraulic gradient is unitless, area is 200x200 m^2
            q_seep = (soil_K[t,:,:])*hf_tot[t,:,:]*(200*200)*((d_arr[qn,:] + soil_thick)/soil_thick)
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
                
                
    base_fn = join(data_dir, 'type'+str(ft), 'r'+str(t).zfill(3)+'_')
    # saving all of the flow at all steps, setbacks is needed to post-process
    Q_out = np.reshape(Q, ((q_in.shape[0]*len(setbacks), xs_levee_smooth.shape[1]+1)))
    np.savetxt(base_fn+'flow.tsv', Q_out)
    # for recharge we want to aggregate across time steps but look at differences across setbacks
    rch_out = np.reshape(np.nansum(rch_hf_arr, axis=0), (len(setbacks)*nrow, ncol))
    np.savetxt(base_fn+'recharge.tsv', rch_out)
    toc = time()
    print((toc-tic)/3600)
    return(Q, rch_hf_arr, d_arr, cell_frac)



###############################################################################
#%% Multiprocess

# realization_recharge(1)

def main():
    pool = Pool(processes=multiprocessing.cpu_count()-2)  # set the processes max number to the number of cpus
    result = pool.map(realization_recharge, range(100))
    pool.close()
    pool.join()
    print(result)
    print('end')

if __name__ == "__main__":
    tic = time()
    main()
    toc = time()
    print('Total time: %.2f minutes' %((toc-tic)/60))