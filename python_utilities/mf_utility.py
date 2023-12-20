"""
mf_utility module. 
Different functions for modflow set up with general python functions
First iteration as a Module April 2023
Author: Andrew Calderwood
"""

import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
import pandas as pd
import os
from os.path import join, exists
import flopy
#%% Model development


def get_layer_from_elev(elev, botm_slice, nlay):
    """  Return uppermost model layer (0-based) occupied at least partly by some elevation data
    Parameters
    ----------
    elev: 1D array (n) with elevations matching model elevation units
    botm: 2D array (nlay, n) with layer elevations of model using same x,y locations at elev1D
    """
    elev_lay = np.zeros(len(elev))
    for k in np.arange(0,nlay-1):
        for j in np.arange(0,len(elev)):
            if botm_slice[k,j] > elev[j]:
                elev_lay[j] = k + 1
    return(elev_lay.astype(int))
    

def param_load(model_ws, file_dir, file_name):
    """ Check if a file exists in the model directory if not copy from Box location"""
    if file_name in os.listdir(model_ws):
        print('exists in model workspace')
        params = pd.read_csv(join(model_ws, file_name))
    else:
        print('added to model workspace')
        params = pd.read_csv(join(file_dir,file_name))
        params.to_csv(join(model_ws,file_name), index=False)
    return(params)


#%% Post-processing

def clean_hob(model_ws, dt_ref):
    hobout = pd.read_csv(join(model_ws,'MF.hob.out'),delimiter=r'\s+', header = 0,names = ['sim_val','obs_val','obs_nam'],
                         dtype = {'sim_val':float,'obs_val':float,'obs_nam':object})
    hobout[['Sensor', 'spd']] = hobout.obs_nam.str.split('p',n=2, expand=True)
    hobout['kstpkper'] = list(zip(np.full(len(hobout),0), hobout.spd.astype(int)))
    hobout = hobout.join(dt_ref.set_index('kstpkper'), on='kstpkper')
    hobout.loc[hobout.sim_val.isin([-1e30, -999.99,-9999]), 'sim_val'] = np.nan
    hobout = hobout.dropna(subset='sim_val')
    hobout['error'] = hobout.obs_val - hobout.sim_val
    hobout['sq_error'] = hobout.error**2
    
    return(hobout)


def get_dates(dis, ref='end'):
    """ Given a MODFLOW DIS file return datetimes given the model start datetime
    input:
    dis = flopy.modflow.ModflowDis object
    ref = 'strt' or 'end' to specify how datetime should be kept
    """
    # specify string characters for timedelta and conversion factors from model unit to seconds
    td_dict = {1:'s',2:'m',3:'h',4:'D'}
    td_uni = 'timedelta64['+td_dict[dis.itmuni]+']'
    tc_dict = {1:1,2:60,3:3600,4:86400}
    # create copies for easier referencing
    perlen = dis.perlen.array.copy()
    nstp = dis.nstp.array.copy()
    tsmult = dis.tsmult.array.copy()
    ## process datetime reference ##
    strt_date = pd.to_datetime(dis.start_datetime)
    # the length of each period should correspond to the correct units
    end_date = (strt_date + pd.Series(perlen.sum()).astype(td_uni))[0]
    # with SS period near 0 no longer minus one
    dates_per = strt_date + (perlen.cumsum()).astype(td_uni)
    perlen_stp = np.repeat(perlen, nstp)
    stp_len = []
    for tn in np.arange(0, len(nstp)):
        # calculate the multiplier for each step
        stp_mult = ((tsmult[tn]*np.ones(nstp[tn])).cumprod()/tsmult[tn])
        # fraction of period for each step times the period length to get the length of each step
        stp_len = np.append(stp_len, perlen[tn]*stp_mult/stp_mult.sum())
    # convert from model unit to seconds to maintain best accuracy with timedelta
    # we want to create dates for the end of each timestep when the output is
    td = (stp_len.cumsum()*tc_dict[dis.itmuni]).astype('timedelta64[s]')
    dates_stps = strt_date + td
    # in some cases it is better to represent output by date of start of each period
    if ref=='strt':
        dates_stps = pd.concat((pd.Series(strt_date), pd.Series(dates_stps[:-1])))
    # get ALL stress periods and time steps list, not just those in the output
    kstpkper = []
    for n,stps in enumerate(nstp):
        kstpkper += list(zip(np.arange(0,stps),np.full(stps,n)))

    dt_ref = pd.DataFrame(dates_stps, columns=['dt'])
    dt_ref['kstpkper'] = kstpkper
    return(strt_date, end_date, dt_ref)
    
def clean_wb(model_ws, dt_ref):
    """Give a volumetric water budget output from MF-OWHM,
    return the water budget with a datetime column and 
    the in and out columns with at least one non-zero value"""
    # load summary water budget
    wb = pd.read_csv(model_ws+'/flow_budget.txt', delimiter=r'\s+')

    wb['kstpkper'] = list(zip(wb.STP-1,wb.PER-1))
    wb = wb.merge(dt_ref, on='kstpkper').set_index('dt')

    # calculate change in storage
    wb['dSTORAGE'] = wb.STORAGE_OUT - wb.STORAGE_IN
    # calculate the cumulative storage change
    wb['dSTORAGE_sum'] = wb.dSTORAGE.cumsum()
    # calculate net groundwater flow
    wb['GHB_NET'] = wb.GHB_IN - wb.GHB_OUT

    # calculate total gw flow, sum GHB, CHD (not uniformly helpful, may be confusing)
    # wb['GW_OUT'] = wb.GHB_OUT + wb.CHD_OUT
    # wb['GW_IN'] = wb.GHB_IN + wb.CHD_IN
    # wb = wb.loc[:,~wb.columns.str.contains('GHB|CHD')]
    
    wb_cols = wb.columns[wb.columns.str.contains('_IN|_OUT')]
    wb_cols = wb_cols[~wb_cols.str.contains('STORAGE|IN_OUT')]
    wb_out_cols= wb_cols[wb_cols.str.contains('_OUT')]
    wb_in_cols = wb_cols[wb_cols.str.contains('_IN')]
    # only include columns with values used
    wb_out_cols = wb_out_cols[np.sum(wb[wb_out_cols]>0, axis=0).astype(bool)]
    wb_in_cols = wb_in_cols[np.sum(wb[wb_in_cols]>0, axis=0).astype(bool)]

    return(wb, wb_out_cols, wb_in_cols)
    
    
    
def read_gage(gagenam):
    """Load the Lake Gage file """
    gage = pd.read_csv(gagenam, skiprows=1, delimiter = r'\s+', engine='python')
    cols = gage.columns[1:-1]
    gage = gage.dropna(axis=1)
    gage.columns = cols
    strt_date = pd.to_datetime(dis.start_datetime)
    gage['dt'] = strt_date+(gage.Time*24).astype('timedelta64[h]')
    gage = gage.set_index('dt')
    gage['dVolume'] = gage.Volume.diff()
    gage['Total_In'] = gage[['Precip.','Runoff','GW-Inflw','SW-Inflw']].sum(axis=1)
    gage['Total_Out'] = gage[['Evap.','Withdrawal','GW-Outflw','SW-Outflw']].sum(axis=1)
    gage['In-Out'] = gage.Total_In - gage.Total_Out
    return(gage)
    
    
def clean_sfr_df(model_ws, dt_ref, pd_sfr=None, name='MF'):
    """Load sfr.out file and create new columns """
    sfrout = flopy.utils.SfrFile(join(model_ws, name+'.sfr.out'))
    sfrdf = sfrout.get_dataframe()
    sfrdf = sfrdf.join(dt_ref.set_index('kstpkper'), on='kstpkper').set_index('dt')
    # convert from sub-daily to daily using mean, lose kstpkper
    sfrdf = sfrdf.groupby(['segment','reach']).resample('D').mean(numeric_only=True)
    sfrdf = sfrdf.reset_index(['segment','reach'], drop=True)
    sfrdf[['row','column']]-=1 # convert to python
    # cmd2cfs = 1/((0.3048**3)*86400) # cubic meters per day to cfs
    sfrdf['month'] = sfrdf.index.month
    sfrdf['WY'] = sfrdf.index.year
    sfrdf.loc[sfrdf.month>=10, 'WY'] +=1
    # add column to track days with flow
    sfrdf['flowing'] = 1
    sfrdf.loc[sfrdf.Qout <= 0, 'flowing'] = 0
    if pd_sfr is not None:
        sfrdf = sfrdf.join(pd_sfr ,on=['segment','reach'],how='inner',lsuffix='_all')

    # create different column for stream losing vs gaining seeapge
    sfrdf['Qrech'] = np.where(sfrdf.Qaquifer>0, sfrdf.Qaquifer,0)
    sfrdf['Qbase'] = np.where(sfrdf.Qaquifer<0, sfrdf.Qaquifer*-1,0 )
    # booleans for plotting
    sfrdf['gaining'] = (sfrdf.gradient == 0)
    sfrdf['losing'] = (sfrdf.gradient >= 0)
    sfrdf['connected'] = (sfrdf.gradient < 1)
    return(sfrdf)

