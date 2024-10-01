# %%
# Calculated with an explicity soil water budget including runoff, evapotranspiration and percolation applied in order to avoid iterative solving. The parameters are all based on SSURGO data, the calculated percolation should be scaled by the VKA to avoid excess percolation.
# - verified on 3/4/2024 that the parcel percolation data covers most of the domain including foothills with small gaps likely due to gaps between parcels, no land use ID or urban land category.

# %%
import sys
from os.path import join, exists, dirname, basename, expanduser
import glob

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


import geopandas as gpd

# %%
# 
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/CosumnesRiverRecharge/python_utilities')

# %%
strt_date = pd.to_datetime('2016-10-1')
end_date = pd.to_datetime('2020-9-30')

# %%
years = pd.date_range(strt_date,end_date+pd.DateOffset(years=1),freq='AS-Oct')
yr_ind = (years-strt_date).days
years = years[:-1]
n=6
# years[n]

# %%
# function can't handle data less than a year
from swb_utility import load_swb_data, yr_lim_from_dates

# %%
uzf_dir = join(gwfm_dir, 'UZF_data')

# %%
perc = load_swb_data(strt_date, end_date, 'field_percolation', uzf_dir)



# %%
AW = load_swb_data(strt_date, end_date, 'field_applied_water', uzf_dir)



# %%
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')
nrow,ncol = dem_data.shape

# %%
# raster cropping will be done in outside script so the only part read in will be the final array
ghb_dir = gwfm_dir+'/GHB_data'

strtyear = strt_date.year
endyear = end_date.year+1
kriged_fall = np.zeros((int(endyear-strtyear),nrow,ncol))
kriged_spring = np.zeros((int(endyear-strtyear),nrow,ncol))

# keep track of which place in array matches to year
year_to_int = np.zeros((endyear-strtyear,2))

for t, year in enumerate(np.arange(strtyear,endyear)):
    # load and place spring kriged data in np array, load spring first
    filename = glob.glob(ghb_dir+'/final_WSEL_arrays/spring'+str(year)+'_kriged_WSEL.tsv')[0]
    # convert from feet to meters
    kriged_spring[t,:,:] = np.loadtxt(filename)*0.3048
    # load and place fall kriged data in np array
    filename = glob.glob(ghb_dir+'/final_WSEL_arrays/fall'+str(year)+'_kriged_WSEL.tsv')[0]
    # convert from feet to meters
    kriged_fall[t,:,:] = np.loadtxt(filename)*0.3048

    year_to_int[t,0] = t
    year_to_int[t,1] = year

# %%
# stack fall and spring before resampling
kriged_arr = np.vstack((kriged_spring, kriged_fall))

# Set kriged water table elevations that are above land surface to land surface minus 15 ft (based on historical levels)
# in floodplain elevations can come up to ground surface
# dem_offset = 15*0.3048
dem_offset = 0
kriged_arr = np.where(kriged_arr>dem_data, dem_data- dem_offset, kriged_arr)


# %%
gde_dir = join(uzf_dir,'shp_GDE_TFT')

# locations of GDes and native vegetation
# the pre-calculated rooting depth assumes deeper roots for GDEs than regular natives
GDE_all = gpd.read_file(join(gde_dir,'GDE_and_native_cell.shp'))

# convert rooting depth to array format for modflow input, hydrographs in wells show drawdown to about 10 m
# below ground so should use 10 m for all gde
ext_dp = np.full((nrow,ncol),2)
ext_dp[(GDE_all.row-1).astype(int), (GDE_all.column-1).astype(int)] = GDE_all.rtg_dp

# %%
## where dtw is on minimum >10 m (GSP number is 30 ft) don't set ET 
# use the year with max elev to insure that extreme years are captured
avg_dtw = dem_data - np.max(kriged_arr,axis=0)

# %%
evt_active = np.zeros((nrow,ncol), dtype='bool')
evt_active[(avg_dtw<10)] = True


# %%
## Potential ETo spatial interpolation from CIMIS
fn = glob.glob(join(uzf_dir,'CIMIS','Cosumnes_dailyET_precip*.csv'))
daily_data = pd.DataFrame()
for file in fn:
    new_data = pd.read_csv(file, index_col = ['Date'], parse_dates = True)
    daily_data = pd.concat((daily_data, new_data))
# units of mm
data_in = daily_data[daily_data['Stn Name']=='Fair Oaks']
# clean up data so columns are by location, units of Precip are in mm
rain_in = data_in.pivot_table(index = 'Date', columns = 'Stn Name', values = 'Precip (mm)')
rain_m = rain_in/1000
# clean up data so columns are by location, units of Precip are in mm
ETo_in = data_in.pivot_table(index = 'Date', columns = 'Stn Name', values = 'ETo (mm)')
ETo_m = ETo_in/1000

# create array for every period of rainfall
rain_df = rain_m[strt_date:end_date].resample('D').interpolate('zero')['Fair Oaks']
rain_arr = np.repeat(np.repeat(np.reshape(rain_df.values, (rain_df.shape[0],1,1)), nrow, axis=1),ncol, axis=2)


# %% [markdown]
# The recharge from field percolation needs to be adjusted to account for where groundwater et may steal from the percolation (no double counting) and in spots where it is expected that water will sit for longer (no runoff losses).
# - clipping to days above the 23 cms threshold brings it down considerably (17,500 to 12,500) and increasing to 75 cms brings it down further (10,000)
# - the spatial distribution makes sense as well because locations along the river with riparian plants see precip all the time and other places see higher than average rates ( should assume 75 cms rate to align with inundation with levees from Whipple)

# %%

# %%
inflow = pd.read_csv(join(gwfm_dir, 'SFR_data', 'MB_daily_flow_cfs.csv'), index_col = 'datetime', parse_dates = True)
# inflow = pd.read_csv(sfr_dir+'MB_daily_flow_cfs_2010_2019.csv', index_col = 'datetime', parse_dates = True)
# covnert flow from cubic feet per second to cubic meters per day
cfs2cmd = (86400/(3.28**3))
inflow['flow_cmd'] = inflow.flow_cfs * cfs2cmd
inflow = inflow.loc[strt_date:end_date]

# %%
# for larger floodplain area, only use rainfall when flow is above threshold
fp_days = (inflow.flow_cmd>75*86400).values
# evt_active array for each stress period
evt_active_all = np.repeat(np.reshape(evt_active, (1, nrow,ncol)), len(fp_days),axis=0).astype(bool)
evt_active_all[~fp_days] = False


# %%
# where GDEs are active the rain should be directly applied instead of SWB percolation
# as we don't want to double count ET
adj_perc = perc.copy()
# initially rain was applied only where there was a deep rooting depth
et_rain_bool = (evt_active)&(ext_dp>2)
adj_perc[:, et_rain_bool] = rain_arr[:, et_rain_bool]
adj_perc_min = adj_perc.copy()
# model results show under prediction of levels in the floodplain
# so it may make sense to use rain everywhere there is ET
# only do it for days with floodplain inundation
adj_perc[evt_active_all] = rain_arr[evt_active_all]
adj_perc_fp = adj_perc.copy()
# do it for all time
adj_perc[:, evt_active] = rain_arr[:, evt_active]
adj_perc_max = adj_perc.copy()

# might be too much water for steady state for whole floodplain

adj_perc = adj_perc_fp.copy()

# # remove excess recharge in foothills that would become runoff in reality
# perc_adj_bool = adj_lowK_arr[drain_layer]
# vka_adj = vka[drain_layer,  perc_adj_bool]
# adj_perc[:, perc_adj_bool] = np.where(adj_perc[:, perc_adj_bool] > vka_adj, vka_adj, adj_perc[:, perc_adj_bool])
# cumulative line doesn't change much from the final adj_perc 
# there are only 152 cells/times that are reduced

# %%
af_scale = 200*200/(0.3048**3)/43560/1000/(end_date.year-strt_date.year)
# plt.plot(agETc.sum(axis=(1,2)).cumsum(),label='ag ET')
plt.plot(AW.sum(axis=(1)).cumsum()*af_scale,label='ag AW')
# plt.plot(natETc.sum(axis=(1,2)).cumsum(),label='native ET')

plt.plot(perc.sum(axis=(1,2)).cumsum()*af_scale, label='Perc')
plt.plot(adj_perc_min.sum(axis=(1,2)).cumsum()*af_scale, label='Perc adjusted for deep ET')
plt.plot(adj_perc_fp.sum(axis=(1,2)).cumsum()*af_scale, label='Perc adjusted for all ET and flood days')
plt.plot(adj_perc_max.sum(axis=(1,2)).cumsum()*af_scale, label='Perc adjusted for all ET')
# plt.plot(adj_perc.sum(axis=(1,2)).cumsum()*af_scale, label='Perc adjusted for foothill')
plt.ylabel('TAF/year')
plt.legend()
plt.show()
# 
