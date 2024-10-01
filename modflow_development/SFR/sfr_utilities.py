# +
''' sfr_utilities module to further clean the SFR stream segment data and prepare it for use in the model
first iteration: 2024-8-5
'''

import numpy as np
import pandas as pd
import geopandas as gpd


# -

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

