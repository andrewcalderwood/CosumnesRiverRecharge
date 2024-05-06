# %%
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
# %% Model development


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


def zone_clean(cbc,zon,  kstpkper, strt_date):
    """ assumes totim units in days
    """
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

def sfr_load_hds(hdobj, grid_sfr, plt_dates, dt_ref):
    """ Pull the groundwater elevation from below the stream cells
    Input:
    hdobj: flopy head object
    grid_sfr: pd.dataFrame of reach data
    plt_dates: list of dates to sample
    dt_ref: reference between dates and stress periods
    Output:
    sfr_heads: heads for each row,column with the sfr layer
    avg_heads: head averaged for top 10 layers
    """
    # runs pretty quickly with hdobj.get_data
    sfr_heads = np.zeros((len(plt_dates), len(grid_sfr)))
    avg_heads = np.zeros((len(plt_dates), len(grid_sfr)))
    for n, plt_date in enumerate(plt_dates):
        spd = dt_ref.loc[dt_ref.dt==plt_date, 'kstpkper'].values[0]
    
        head = hdobj.get_data(spd)
        head = np.ma.masked_where(head ==-999.99, head)
        sfr_heads[n,:] = head[grid_sfr.k, grid_sfr.i, grid_sfr.j]
        # pull head for top 10 layers to compare
        avg_heads[n,:] = np.mean(head[:10, grid_sfr.i, grid_sfr.j], axis=0)
    return(sfr_heads, avg_heads)

