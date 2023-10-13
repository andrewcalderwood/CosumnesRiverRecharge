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
#%% Model development

from scipy.stats import ttest_rel, linregress


def run_stats(wb, wb0, term, season=None, freq='monthly', plot=False):
    if season == 'Wet':
        wet_months = [11,12,1,2,3,4]
        wb = wb[wb.index.month.isin(wet_months)]
        wb0 = wb0[wb0.index.month.isin(wet_months)]
    elif season =='Dry':
        dry_months = [5,6,7,8,9,10]
        wb = wb[wb.index.month.isin(dry_months)]
        wb0 = wb0[wb0.index.month.isin(dry_months)]
    elif season=='Fall':
        fall_months=[9,10,11]
        wb = wb[wb.index.month.isin(fall_months)]
        wb0 = wb0[wb0.index.month.isin(fall_months)]
    if freq=='annual':
        a = wb0.resample('AS-Oct').sum()[term].values
        b = wb.resample('AS-Oct').sum()[term].values
    elif freq=='monthly':
        a = wb0.resample('MS').sum()[term].values
        b = wb.resample('MS').sum()[term].values
    elif freq=='daily':
        a = wb0.resample('D').sum()[term].values
        b = wb.resample('D').sum()[term].values
        
    t_out = ttest_rel(a, b)

    t_df = pd.DataFrame([t_out.statistic, t_out.pvalue]).transpose()
    t_df.columns=['statistic','pvalue']
    t_df['term'] = term
    t_df['freq'] = freq
    t_df['season'] = season
    t_df['mean_a'] = np.mean(a)
    t_df['mean_b'] = np.mean(b)
    t_df['perc_diff_in_means'] = 100*(np.mean(a)-np.mean(b))/np.abs((np.mean(a)+np.mean(b))/2)

    # rounding to clean up output
    t_df.statistic = t_df.statistic.round(3)
    t_df.pvalue = t_df.pvalue.round(4)
    t_df.perc_diff_in_means = t_df.perc_diff_in_means.round(2)

    # if pvalue is insignificant then don't include
    t_df.loc[t_df.pvalue>=0.05,'perc_diff_in_means'] = '-'

    
    if plot:
        slope, intercept, r_value, p_value, std_err = linregress(a, b)
        print('T-test statistic: %.2f' %t_out.statistic, 'and pvalue: %.4f' %t_out.pvalue)

        plt.scatter(a, b)
        x_range = np.array([[np.min((a,b))], [np.max((a,b))]])
        plt.plot(x_range, slope*x_range + intercept, color='black', linewidth=1)
        plt.annotate('y = '+str(np.round(slope,3))+'x + '+ str(np.round(intercept,2)), (0.1,0.8), xycoords='axes fraction')
        plt.title(term)
        plt.ylabel('No Reconnection')
        plt.xlabel('Baseline')
    return(t_df)