# +
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
import seaborn as sns

# standard geospatial python utilities
import geopandas as gpd


# +
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
    
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study
proj_dir = gwfm_dir + '/Oneto_Denier/'
dat_dir = proj_dir+'Stream_level_data/'

fig_dir = proj_dir+'/Streambed_seepage/figures/'


# +
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

add_path(doc_dir+'/GitHub/flopy')
import flopy 
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)
from mf_utility import get_dates, get_layer_from_elev, clean_wb
from map_cln import gdf_bnds, plt_cln
# -

# scenario = '' # baseline, levee removal occurred in 2014
# create identifier for scenario if levee removal didn't occur
scenario = 'no_reconnection'

# +
ext_dir = 'F:/WRDAPP'
c_dir = 'C:/WRDAPP'
if os.path.exists(ext_dir):
    loadpth = ext_dir 
elif os.path.exists(c_dir):
    loadpth = c_dir 
loadpth +=  '/GWFlowModel/Cosumnes/Stream_seepage'

upscale = 'upscale4x_'
# model_nam = 'oneto_denier_'+upscale+'2014_2018'
model_nam = 'oneto_denier_'+upscale+'2014_2020'
model_ws = join(loadpth,model_nam)

if scenario != '':
    model_ws += '_' + scenario
    
# model_ws = join(loadpth,'parallel_oneto_denier','realization000')
load_only = ['DIS','UPW','SFR','OC', "EVT"]
m = flopy.modflow.Modflow.load('MF.nam', model_ws= model_ws, 
                                exe_name='mf-owhm.exe', version='mfnwt',
                              load_only=load_only,
                              )


# -

model_ws0 = join(loadpth,model_nam)


# +
model_grp = 'inset_oneto_denier'
grid_dir = join(gwfm_dir, 'DIS_data/streambed_seepage/grid')
grid_fn = join(grid_dir, model_grp,'rm_only_grid.shp')
grid_p = gpd.read_file(grid_fn)
grid_p.crs='epsg:32610'
m_domain = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [grid_p.unary_union], crs=grid_p.crs)

sfrdf = pd.DataFrame(m.sfr.reach_data)
grid_sfr = grid_p.set_index(['row','column']).loc[list(zip(sfrdf.i+1,sfrdf.j+1))].reset_index(drop=True)
grid_sfr = pd.concat((grid_sfr,sfrdf),axis=1)

# -

# identify the different kind of facies
coarse_cutoff = 70 # sandy mud is 20 m/d, sand is 120 m/d
seep_vka = np.copy(m.upw.vka.array)
coarse = (seep_vka > coarse_cutoff).astype(int)
fine = (seep_vka < coarse_cutoff).astype(int)


zon_color_dict = pd.read_excel('mf_wb_color_dict.xlsx',sheet_name='mf_wb_dict', header=0, index_col='flux',comment='#').color.to_dict()
zon_name_dict = pd.read_excel('mf_wb_color_dict.xlsx',sheet_name='mf_wb_dict', header=0, index_col='flux',comment='#').name.to_dict()


# dictionary to transfer flopy zonebudget names to owhm format for consistency
zb_alt = pd.read_excel('mf_wb_color_dict.xlsx',sheet_name='flopy_to_owhm', header=0, 
                       index_col='flopy',comment='#').owhm.to_dict()


lak_shp = join(gwfm_dir,'LAK_data/floodplain_delineation')
# shapefile rectangle of the area surrounding the Dam within about 5 cells
lak_gdf = gpd.read_file(join(lak_shp,'LCRFR_ModelDom_2017/LCRFR_2DArea_2015.shp' )).to_crs('epsg:32610')
lak_cells = gpd.sjoin(grid_p,lak_gdf,how='right',predicate='within')
# create array of floodplain cells
zon_lak = np.zeros((grid_p.row.max(),grid_p.column.max()),dtype=int)
zon_lak[lak_cells.row-1,lak_cells.column-1]=1

# A distance of 500 m includes the whole reconnected floodplain while only part of the floodplain is included with 400 m.

# create riparian corridor boundary (500 m)
riparian_gdf = gpd.GeoDataFrame(pd.DataFrame([0]), 
                                geometry = [grid_sfr.buffer(500).unary_union], crs=grid_p.crs)
# identify grid cells in the riparian corridor
riparian_grid = gpd.sjoin(grid_p, riparian_gdf).drop(columns=['index_right'])
zon_rip = np.zeros((grid_p.row.max(),grid_p.column.max()),dtype=int)
zon_rip[riparian_grid.row-1, riparian_grid.column-1]=1

fig,ax = plt.subplots()
riparian_gdf.plot(ax=ax, color='None')
lak_gdf.plot(ax=ax, color='None')


# # Zonebudget for understanding riparian zone impacts
# Could compare the ET in the floodplain, within 2-3 cells (300-500 m) of the river, overall?
# - the benefit of a riparian zonebudget would be that rather than interpreting increased streambed recharge as ET, pumping, gw storage, or GW out it is simplified to gw out with limited near bank storage. The larger gw system might mask some localized effects of floodplain recharge unless we subset.
#
# Using only the floodplain leaves the stream out of the water budget (only ET, GW and storage)  
# Riparian corridor and model domain can have all components

def zone_clean(cbc,zon,  kstpkper):
    zb = flopy.utils.ZoneBudget(cbc, zon, kstpkper)
    zb_df = zb.get_dataframes()
    # ungroup by timestep
    zb_df = zb_df.reset_index()
    names = zb_df.name.unique()
    zb_df = zb_df.pivot(index = 'totim', columns = 'name',values = 'ZONE_1')
    
    # columns to make negative
    to_cols = zb_df.columns[zb_df.columns.str.contains('TO_')]
    # get net GHB
    zb_df['GHB_NET'] = zb_df.TO_HEAD_DEP_BOUNDS - zb_df.FROM_HEAD_DEP_BOUNDS
    # to storage is gw increase (positive)
    stor_cols = zb_df.columns[zb_df.columns.str.contains('STORAGE')]
    zb_df['dSTORAGE'] = (zb_df.TO_STORAGE - zb_df.FROM_STORAGE)
    zb_df['dSTORAGE_sum'] = zb_df.dSTORAGE.copy().cumsum()
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


hdobj = flopy.utils.HeadFile(join(model_ws, 'MF.hds'))
spd_stp = hdobj.get_kstpkper()


# +
def make_zonebud(zon, name):
    cbc = join(model_ws, 'MF.cbc')
    zb_df, zb_mon = zone_clean(cbc, zon, spd_stp)
    zb_df.to_csv(join(model_ws, 'MF_zonebud_'+name+'_daily.csv'))
    zb_mon.to_csv(join(model_ws, 'MF_zonebud_'+name+'_monthly.csv'))
    print('No reconnection '+name+' facies zonebudget done')
    cbc = join(model_ws0, 'MF.cbc')
    zb_df0, zb_mon0 = zone_clean(cbc, zon, spd_stp)
    zb_df0.to_csv(join(model_ws0, 'MF_zonebud_'+name+'_daily.csv'))
    zb_mon0.to_csv(join(model_ws0, 'MF_zonebud_'+name+'_monthly.csv'))
    print('Baseline '+name+' zonebudget done')

# 


# -

make_zonebud(zon_lak, 'floodplain')

# +
# cbc = join(model_ws, 'MF.cbc')
# zb_df, zb_mon = zone_clean(cbc, zon_lak, spd_stp)
# zb_df.to_csv(join(model_ws, 'MF_zonebud_floodplain_daily.csv'))
# zb_mon.to_csv(join(model_ws, 'MF_zonebud_floodplain_monthly.csv'))
# print('No reconnection floodplain zonebudget done')
# cbc = join(model_ws0, 'MF.cbc')
# zb_df0, zb_mon0 = zone_clean(cbc, zon_lak, spd_stp)
# zb_df0.to_csv(join(model_ws0, 'MF_zonebud_floodplain_daily.csv'))
# zb_mon0.to_csv(join(model_ws0, 'MF_zonebud_floodplain_monthly.csv'))
# print('Baseline floodplain zonebudget done')
# -


make_zonebud(zon_rip, 'riparian')

# +
# cbc = join(model_ws, 'MF.cbc')
# zb_df, zb_mon = zone_clean(cbc, zon_rip, spd_stp)
# zb_df.to_csv(join(model_ws, 'MF_zonebud_riparian_daily.csv'))
# zb_mon.to_csv(join(model_ws, 'MF_zonebud_riparian_monthly.csv'))
# print('No reconnection riparian zonebudget done')
# cbc = join(model_ws0, 'MF.cbc')
# zb_df0, zb_mon0 = zone_clean(cbc, zon_rip, spd_stp)
# zb_df0.to_csv(join(model_ws0, 'MF_zonebud_riparian_daily.csv'))
# zb_mon0.to_csv(join(model_ws0, 'MF_zonebud_riparian_monthly.csv'))
# print('Baseline riparian zonebudget done')
# -


make_zonebud(fine, 'fine')

make_zonebud(coarse, 'coarse')
