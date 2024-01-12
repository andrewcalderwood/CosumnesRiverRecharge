"""
OD_utility module. 
Different functions for modflow set up with general python functions
First iteration as a Module October 2023
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

from scipy.stats import ttest_rel, linregress


def lin_reg_plt(a, b):
    """ Plot scatter plot and linear regression, print ttest results
    """
    slope, intercept, r_value, p_value, std_err = linregress(a, b)

    plt.scatter(a, b)
    x_range = np.array([[np.min((a,b))], [np.max((a,b))]])
    plt.plot(x_range, slope*x_range + intercept, color='black', linewidth=1)
    plt.annotate('y = '+str(np.round(slope,3))+'x + '+ str(np.round(intercept,2)), (0.1,0.8), xycoords='axes fraction')
    plt.ylabel('No Reconnection')
    plt.xlabel('Baseline')
def run_ttest(a, b, term, freq):
    """ Run ttest and summarize to tables
    """
    t_out = ttest_rel(a, b)

    t_df = pd.DataFrame([t_out.statistic, t_out.pvalue]).transpose()
    t_df.columns=['statistic','pvalue']
    t_df['term'] = term
    # t_df['freq'] = freq
    # t_df['season'] = season
    t_df['mean_a'] = np.mean(a)
    t_df['mean_b'] = np.mean(b)
    t_df['perc_diff_in_means'] = 100*(np.mean(a)-np.mean(b))/np.abs((np.mean(a)+np.mean(b))/2)

    # rounding to clean up output
    t_df.statistic = t_df.statistic.round(3)
    t_df.pvalue = t_df.pvalue.round(4)
    t_df.perc_diff_in_means = t_df.perc_diff_in_means.round(2)

    # if pvalue is insignificant then don't include
    t_df.loc[t_df.pvalue>=0.05,'perc_diff_in_means'] = '-'
    return(t_df)

# wet_months = [11,12,1,2,3,4]
# dry_months = [5,6,7,8,9,10]
# fall_months=[9,10,11]

def run_stats(wb, wb0, term, season=None, freq='monthly', plot=False):
    if season == 'Winter':
        months = [12,1,2]
    elif season== 'Spring':
        months = [3,4,5]
    elif season =='Summer':
        months = [6,7,8]
    elif season=='Fall':
        months=[9,10,11]
    if season is not None:
        wb = wb[wb.index.month.isin(months)]
        wb0 = wb0[wb0.index.month.isin(months)]
    if freq=='annual':
        a = wb0.resample('AS-Oct').sum()[term].values
        b = wb.resample('AS-Oct').sum()[term].values
    elif freq=='monthly':
        a = wb0.resample('MS').sum()[term].values
        b = wb.resample('MS').sum()[term].values
    elif freq=='daily':
        a = wb0.resample('D').sum()[term].values
        b = wb.resample('D').sum()[term].values

    t_df = run_ttest(a, b, term, freq)
    t_df['freq'] = freq
    t_df['season'] = season
    if plot:
        print('T-test statistic: %.2f' %t_df.statistic.iloc[0], 'and pvalue: %.4f' %t_df.pvalue.iloc[0])
        lin_reg_plt(a, b)
        plt.title(term)
        plt.show()

    return(t_df)

# need to fix as this doesn't account for sfr_rows, sfr_cols or plt_segs
def plot_vka(vka, ax, sfr_nodata, k_max):
    sfr_hk = vka[:k_max][:, sfr_rows, sfr_cols]
    sfr_hk = np.ma.masked_where(sfr_nodata[:k_max], sfr_hk)
    im = ax.imshow(sfr_hk, norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax), 
                   aspect='auto', cmap='viridis')
    # plt.xticks([]);
    ax.set_yticks(ticks = np.arange(1,k_max,5), labels=m.dis.botm.array[:,0,0][:k_max:5]);
    ax.set_xticks(ticks = np.arange(0, len(plt_segs),10), labels=np.arange(0, len(plt_segs),10), rotation=90)
    return sfr_hk, im


def clean_sfr_df(model_ws,  drop_iseg, dt_ref, m_name ='MF'):
    """ specific to OD parallle realizations because of 
    grid_sfr load and num_coarse specification
    """
    ## load sfr reach data ##
    grid_sfr = pd.read_csv(model_ws+'/grid_sfr.csv')
    # remove stream segments for routing purposes only
    grid_sfr = grid_sfr[~grid_sfr.iseg.isin(drop_iseg)]
    pd_sfr = grid_sfr.set_index(['iseg','ireach'])[['rchlen','strtop', 'facies', 'color', 'strthick']]
    pd_sfr['Total distance (m)'] = pd_sfr['rchlen'].cumsum()
    num_coarse = int(grid_sfr.facies.isin(['Gravel','Sand']).sum())
    
    ## load sfr out file ##
    sfrout = flopy.utils.SfrFile(join(model_ws, m_name+'.sfr.out'))
    sfrdf = sfrout.get_dataframe()
    sfrdf = sfrdf.join(dt_ref.set_index('kstpkper'), on='kstpkper').set_index('dt')
    # convert from sub-daily to daily using mean, lose kstpkper
    sfrdf = sfrdf.groupby(['segment','reach']).resample('D').mean(numeric_only=True)
    sfrdf = sfrdf.reset_index(['segment', 'reach'], drop=True)
    sfrdf[['row','column']] = sfrdf[['row','column']].astype(int) - 1 # convert to python
    
    ## join sfr out to reach data ##
    sfrdf = sfrdf.join(pd_sfr ,on=['segment','reach'],how='inner',lsuffix='_all')
    sfrdf['num_coarse'] = num_coarse
    
    ## data transformation for easier manipulation ##
    sfrdf['month'] = sfrdf.index.month
    sfrdf['WY'] = sfrdf.index.year
    sfrdf.loc[sfrdf.month>=10, 'WY'] +=1
    # create column to calculate days flowing
    sfrdf['flowing'] = 1
    sfrdf.loc[sfrdf.Qout <= 0, 'flowing'] = 0
    
    # create different column for stream losing vs gaining seeapge
    sfrdf['Qrech'] = np.where(sfrdf.Qaquifer>0, sfrdf.Qaquifer,0)
    sfrdf['Qbase'] = np.where(sfrdf.Qaquifer<0, sfrdf.Qaquifer*-1,0 )
    # booleans for plotting
    sfrdf['gaining'] = (sfrdf.gradient <= 0)
    sfrdf['losing'] = (sfrdf.gradient >= 0)
    sfrdf['connected'] = (sfrdf.gradient < 1)
    return(sfrdf)

