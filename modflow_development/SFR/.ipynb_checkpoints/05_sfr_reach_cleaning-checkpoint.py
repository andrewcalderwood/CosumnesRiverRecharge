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
# Cosumnes SFR second step pre-processing. Cleans up the Cosumnes and Deer Creek stream data into one stream network dataset prepared for the SFR package.
#

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import sys
import glob
from importlib import reload

import pandas as pd
import numpy as np
from scipy.stats import hmean, gmean

# import calendar
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import geopandas as gpd

# # mapping utilities
# import contextily as ctx
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as fm



# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
print(gwfm_dir)
sfr_dir = gwfm_dir+'/SFR_data/'


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')

import flopy 

from importlib import reload
# importlib.reload
# reload(flopy)

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)


# %%

# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')

# %%
# can auto-calculate delr if we assume grid is square
delr = np.sqrt(grid_p.area.iloc[0]).round(0)

# %%
# original method for sfr grid cells
# grid_sfr_in = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')
# slightly updated method for SFR grid cells (flexible for different streams)
grid_sfr_in = gpd.read_file(sfr_dir+'/final_grid_sfr/Cosumnes_River_sfr.shp')
grid_sfr_in['rname'] = 'Cosumnes River'
# mb_seg = 1 # default value of michigan bar segment is 1


# %%
## could insert Deer Creek reaches here, make them the upstream most (2024-3-7)
deer_ck = gpd.read_file(sfr_dir+'/final_grid_sfr/Deer_Creek_sfr.shp')
deer_ck['rname'] = 'Deer Creek'

# deer creek will segment 1 now
deer_ck['iseg'] = 1
deer_ck['reach_order'] = np.arange(1, len(deer_ck)+1)
deer_ck['reach'] = np.arange(1, len(deer_ck)+1)

# define michigan bar segment (1-based)
mb_seg = 2

# %%
fp_logger = pd.read_csv(join(gwfm_dir,'LAK_data','floodplain_logger_metadata.csv'))
fp_logger = gpd.GeoDataFrame(fp_logger, geometry = gpd.points_from_xy(fp_logger.Easting, fp_logger.Northing), crs='epsg:32610')
# find grid cell it is within
fp_grid = gpd.sjoin(fp_logger, grid_p, how='left',predicate='within')

# %%
# using floodplain logger locations to update XSlocs
OD_locs = fp_grid[['Logger Location','geometry']].copy()
OD_locs = OD_locs[OD_locs['Logger Location'].isin(['OD_Excavation','SwaleBreach_1'])]


# %%
XS8pt = pd.read_csv(sfr_dir+'8pointXS.csv')
XSlocs = gpd.read_file(sfr_dir+'8pointXS_locs/8pointXS_locs.shp')
# new shapefile with an extra point for blodgett dam as site 16.5
# XSlocs = gpd.read_file(gwfm_dir+'/Blodgett_Dam/geospatial/8pointXS_locs/8pointXS_locs.shp')
XSlocs.crs = 32610
# this is better done manually (causes weird numbering)
# split segments where Oneto-Denier floodplain starts and ends for easier diversion to LAK
# XSlocs = pd.concat((XSlocs, OD_locs.drop(columns='Logger Location')))
# identify the nearest XS to the OD sites and update the site name to improve referencing uniqueness
OD_XSlocs = gpd.sjoin_nearest(XSlocs, OD_locs, how='right').drop(columns=['index_left','Logger Location'])
OD_XSlocs['Site'] += [0.4, 0.8]
XSlocs = pd.concat((XSlocs, OD_XSlocs))

# it would be better to do the sjoin_nearest after creating an udpating XSg
# XSg  = gpd.sjoin(grid_sfr, XSlocs, how = "inner", predicate= "contains", lsuffix = 'sfr',rsuffix = 'xs')
XSg_in  = gpd.sjoin_nearest(grid_sfr_in, XSlocs, how = "right", lsuffix = 'sfr',rsuffix = 'xs')

# sort by site to make sure any XS added are properly included
XSg_in = XSg_in.sort_values('Site')

# to account for the iseg splitting added by Oneto-Denier forward fill information on XS
# XSg_in['Site'] = XSg_in['Site'].ffill(axis=0)

XSg_in['iseg'] = np.arange(2,len(XSg_in)+2) # add the segment that corresponds to each cross section


# %%
# XSg = XSg_in.copy()
# if blodgett == 'none':
#     # if no blodgett dam scenario then remove the extra cross section
#     XSg = XSg.loc[(XSg.Site!=16.5)]
#     XSg = XSg.loc[(XSg.Site!=16.2)]
# #     XSg = XSg.loc[XSg.Site!=16.2]
# elif blodgett == 'actual':
#     XSg_side = XSg.loc[XSg.Site==16.5]
#     XSg_side.loc[:,'Site'] = 16.4
#     XSg = XSg.append(XSg_side)
# elif blodgett == 'design':
#     # may or may not want to remove the segment before
#     XSg = XSg.loc[(XSg.Site!=16.2)]

# if the scneario is the restructured or designed dam then no change in the segments is necessary


# %%

# if blodgett == 'dam':
#     # designed scenario flow through dam only
#     new_xs = pd.read_csv(gwfm_dir+'/Blodgett_Dam/geospatial/02_designed_XS.csv', skiprows=1)
# elif blodgett =='actual':
#     # current situation, flow around dam and after dam
#     new_xs = pd.read_csv(gwfm_dir+'/Blodgett_Dam/geospatial/03_actual_XS.csv', skiprows=1)
# elif blodgett =='new':
#     # depending scenario, use different input cross sections for 16.5
#     new_xs = pd.read_csv(gwfm_dir+'/Blodgett_Dam/geospatial/01_New_wide_XS.csv')

# if there is a scneario then need to add the new XS
# if blodgett != 'none':
#     XS8pt = pd.concat([XS8pt,new_xs],axis=1)

# %% [markdown]
# ## Oneto-Denier floodplain connection segments
# For the floodplain connection segments, the Site that it is connected to does not matter because they conductivity will be set to zero so the depth computed will not impact calculations. The stream bottom elevation will have an impact on the return flow because it must be below the lake bottom.

# %%
## new version from Oneto-Denier
XSg = XSg_in.copy()

# if scenario != 'no_reconnection':
# identify XS to be copied for diversion reaches
fp_grid_xs = fp_grid[['Logger Location','geometry']].copy()
fp_grid_xs = fp_grid_xs.sjoin_nearest(XSg.reset_index(), how='inner') #.drop(columns=['index_right'])
# od_breach is the sensor location where the breach was made in the levees for flow to leave the river
od_breach = fp_grid_xs[fp_grid_xs['Logger Location']=='OD_Excavation'].copy()
od_breach['Site'] -= 0.2 # adjust Site to set sorting order
# swale is the return to downstream so it should be just upsream of the segment if possible
od_swale = fp_grid_xs[fp_grid_xs['Logger Location']=='SwaleBreach_1'].copy()
od_swale['Site'] -= 0.01 # adjust Site to set sorting order

# need to adjust elevation so transfer segment from floodplain diversion to stream is positive slope
# return takes flow from the breach and returns it to the channel
od_return = od_breach.copy()
od_return['Site'] += 0.1 # adjust Site to set sorting order
od_return.z_min = XSg[XSg.Site==int(np.ceil(od_return.Site).values[0])].z_min.max() + od_return.slope*delr
# need another segment to transfer water between the breach and lake to avoid doubling
od_lake = od_return.copy()
od_lake['Site'] -= 0.05
od_lake.z_min = od_breach.z_min - od_lake.slope*delr

# add reaches for diversion
# XSg = pd.concat((XSg.reset_index(), od_breach, od_return, od_swale)) #, od_swale (old)
# need a segment to control diversion to lake to avoid doubling water
XSg_all = pd.concat((XSg.reset_index(), od_breach, od_return, od_lake, od_swale)) #, od_swale
# else:
#     XSg = XSg.reset_index().copy()

# %%
# this needs to be done once the sfr reaches are in place
# redefine Site/iseg
def clean_reach_routing(XSg):
    XSg['reach_order'] = np.arange(0, len(XSg)) # fix reach order
    
    # iterate over the added segments to fix reach ordering
    for xsn in XSg[~XSg['Logger Location'].isna()].Site:
    # for xsn in XSg[~XSg['Logger Location'].isna()].Site:
        rn = XSg.loc[XSg.Site==xsn,'reach_order'].values[0]
        if XSg.loc[XSg.reach_order==rn-1, 'iseg'].values== XSg.loc[XSg.reach_order==rn, 'iseg'].values:
            XSg.loc[XSg.Site>=xsn, 'iseg'] += 1
    # # now make sure to add 1 to downstream segments only if theyare equal
    for xsn in XSg[~XSg['Logger Location'].isna()].Site:
        rn = XSg.loc[XSg.Site==xsn,'reach_order'].values[0]
        if XSg.loc[XSg.reach_order==rn, 'iseg'].values== XSg.loc[XSg.reach_order==rn+1, 'iseg'].values:
            XSg.loc[XSg.reach_order>rn, 'iseg'] += 1
    # need to fix reach numbering after fixing segment numbers
    # min_reach = XSg.groupby('iseg').min(numeric_only=True)['ireach']
    # for ns in min_reach[min_reach!=1].index:
    #     XSg.loc[XSg.iseg==ns,'ireach'] = np.arange(1,np.sum(XSg.iseg==ns)+1)
    return(XSg)
    
XSg  = XSg_all.sort_values(['Site','reach']).copy()
XSg = clean_reach_routing(XSg)
# need to start all segments one further from Michigan bar (they assume MB starts at 1)
XSg.iseg += mb_seg-1


# %%
def make_xs_sfr(grid_sfr_in, XSg, mb_seg):
    xs_sfr = grid_sfr_in.merge(XSg[['row','column','Logger Location','Site', 'iseg']],how='left')
    
    # specify reach 1 will have iseg from Michigan Bar icalc=4
    # after adding Deer Creek, Michigan Bar is reach 2
    xs_sfr.loc[xs_sfr.reach==1,'iseg'] = mb_seg
    xs_sfr = xs_sfr.sort_values(['reach', 'iseg'])
    # forward fill iseg numbers 
    xs_sfr.iseg = xs_sfr.iseg.ffill()
    # add forward fill of sites for Deer Creek split of segment
    xs_sfr.Site = xs_sfr.Site.ffill()

    # rename old reach numbers to save
    xs_sfr = xs_sfr.rename(columns={'reach':'reach_order'})
    # specify new reach number for each segment
    xs_sfr['reach'] = 1
    for ns, seg in enumerate(xs_sfr.iseg.unique()):
        xs_sfr.loc[xs_sfr.iseg==seg,'reach'] = np.arange(1,(xs_sfr.iseg==seg).sum()+1)
        
    # get total lengths ( should be separate for Deer Creek and Cosumnes)
    xs_sfr['dist_m'] = xs_sfr.length_m.cumsum()
    xs_sfr.dist_m -= xs_sfr.dist_m.iloc[0]
    return(xs_sfr)




# %%
def find_seg_join(xs_sfr, deer_ck):
    ## find where deer creek enters the Cosumnes River
    # use last reach of deer creek (centroid to only select nearest stream cell edge)
    deer_ck_end = deer_ck[deer_ck.reach_order==deer_ck.reach_order.max()].copy()
    deer_ck_end.geometry = deer_ck_end.centroid
    # clean up columns
    deer_ck_end = deer_ck_end[['row','column','z_m', 'z_min','geometry']]
    # find nearest Cosumnes River reaches
    deer_ck_join = gpd.sjoin_nearest(xs_sfr[xs_sfr.iseg!=1], deer_ck_end, how='right')
    # use the maximum segment number to represent the most likely downstream segment
    dc_seg, dc_rch, dc_rch_order = deer_ck_join[['iseg','reach','reach_order']].iloc[0].astype(int)
    return(dc_seg, dc_rch, dc_rch_order)



# %%
def update_rch_seg(xs_sfr, dc_seg, dc_rch):
    ## update reach ordering for deer creek inflow
    # filter to the  downstream stream segments
    # then filter to the downstream reaches of the segment
    dc_rch_adj = (xs_sfr.iseg==dc_seg)&(xs_sfr.reach>=dc_rch)
    dc_seg_adj = (xs_sfr.iseg > dc_seg)|dc_rch_adj
    # add one to the segments that need adjusting
    xs_sfr.loc[dc_seg_adj, 'iseg'] += 1
    # update reach ordering to be one to the next order
    xs_sfr.loc[dc_rch_adj,'reach'] = np.arange(1, dc_rch_adj.sum()+1)
    return(xs_sfr)


# %%
xs_sfr = make_xs_sfr(grid_sfr_in, XSg, mb_seg)
# join deer creek and Cosumnes River sfr datasets (2024-3-7)
xs_sfr = pd.concat((deer_ck, xs_sfr))
# identify where the reaches need an update
dc_seg, dc_rch, dc_rch_order = find_seg_join(xs_sfr, deer_ck)
# update the segment and reach numbering
xs_sfr = update_rch_seg(xs_sfr, dc_seg, dc_rch)
# update dc_seg to account for the fact it now starts one below
dc_seg += 1

# %%
# # modflow assumes a segment goes into reach 1 of a new segment so it's pulling a higher elevation
# fig,ax = plt.subplots()
# deer_ck_join.plot(ax=ax, color='red')
# xs_sfr[xs_sfr.iseg.isin([dc_seg, dc_seg+1])].plot('iseg', ax=ax)

# %%
xs_sfr_out = xs_sfr.drop(columns=['geometry'])
xs_sfr_out.to_csv(join(gwfm_dir, 'SFR_data','xs_sfr.csv'),index=False)

# %%
xs_sfr_out
