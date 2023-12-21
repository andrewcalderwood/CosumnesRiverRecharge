"""
flopy_utility module. 
Different functions for working with flopy objects/classes set up with general python functions
First iteration as a Module December 2023
Author: Andrew Calderwood
"""

# import matplotlib.pyplot as plt
# import geopandas as gpd
import numpy as np
import pandas as pd
# import os
# from os.path import join, exists
import flopy
#%% Model development


def spd_2_arr(sp_data, sp_col, dis):
    """Given the stress_period_data from flopy return the data in an array format
    """
    # convert pumping to array
    arr = np.zeros((dis.nper,dis.nrow,dis.ncol))
    for n in np.arange(0,dis.nper):
        data_n = sp_data[n]
        # only index array if there is data for a stress period
        if data_n is not None:
            arr[n, data_n.i, data_n.j] += data_n[sp_col]
    return(arr)

def reach_data_gdf(sfr, grid_p):
    """ Join the reach and routing segment data to model grid geodataframe"""
    sfr_reach = pd.DataFrame(sfr.reach_data)
    grid_sfr = grid_p.set_index(['row','column']).loc[list(zip(sfr_reach.i+1,sfr_reach.j+1))].reset_index(drop=True)
    grid_sfr = pd.concat((grid_sfr,sfr_reach),axis=1)

    # load segment data to identify where diversions are
    sfr_seg = pd.DataFrame(sfr.segment_data[0])
    sfr_seg = sfr_seg[['nseg','outseg','iupseg']].rename(columns={'nseg':'iseg'})
    grid_sfr = grid_sfr.merge(sfr_seg)
    return(grid_sfr)

def get_avg_head(hdobj, idx):
    """
    Return the spatially averaged head for the maximum head at the input locations (idx)
    hdobj: flopy head object
    idx: list of tuples as (layer, row, column)
    """
    # get heads under the lake
    ts = hdobj.get_ts(idx)
    ts_df = pd.DataFrame(ts, columns=['totim']+idx)
    ts_df = ts_df.set_index('totim')
    ts_df = ts_df.melt(ignore_index=False)
    ts_df[['k','i','j']] = ts_df.variable.tolist()
    ts_df = ts_df.drop(columns='variable') # drop to speed up groupby
    head = ts_df.groupby(['totim','i','j']).max().groupby('totim').mean()
    return head

def find_active_ET(m, dt_ref):
	"""
	Identify cells and times when groundwater head is above ET rooting depth
	"""
	et_surf = m.evt.surf.array[0,0]
	ext_dp = m.evt.exdp.array[0,0]
	# bottom elevation of roots
	et_botm =et_surf - ext_dp

	et_row, et_col = np.where(ext_dp>2)
	et_lay = get_layer_from_elev(et_botm[et_row, et_col], m.dis.botm[:, et_row, et_col], m.dis.nlay)
	

	hdobj = flopy.utils.HeadFile(join(m.model_ws, 'MF.hds'))
	et_act = np.zeros((m.dis.nper, m.dis.nrow, m.dis.ncol))

	for t in np.arange(0,m.dis.nper):
	    head = hdobj.get_data(dt_ref.kstpkper.iloc[t])
	    head = np.ma.masked_where(head==-999.99, head)
	    # identify which GDE ET cells would active based on head
	    b = head[et_lay, et_row, et_col] > et_botm[et_row, et_col]
	    et_act[t, et_row[b], et_col[b]] = 1
	return et_act