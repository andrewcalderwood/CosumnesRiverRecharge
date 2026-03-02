"""
hyd_utility module. 
Different functions for working with numeric/tabular hydrologic data with general python functions
First iteration as a Module December 2023
Author: Andrew Calderwood
"""

import pandas as pd
import numpy as np  


def dt_2_wy(dt):
	"""Given a 1D array or series of dates, return the water year
	Input: 1D array or series of dates
	Output: array of water year"""
	dt = pd.Series(dt)
	yr = dt.dt.year
	wy = np.where(dt.dt.month>=10, yr+1, yr)
	return(wy)

def mark_outlier(df):
    """ Given a pd series flag those outside the whiskers of a box plot
    as outliers
    Input:
    df:Series of data
    Output:
    df: series of data with outliers flagged"""
    # plot quantiles on the line plot
    quart = df.quantile([.25,.75])
    # median = df.quantile([.5])
    # calculate whiskers
    iqr = quart.loc[0.75]-quart.loc[0.25]
    # 1.5 x the whole interquartile range
    whisker_75 = quart.loc[0.75] + iqr*1.5
    whisker_25 = quart.loc[0.25] - iqr*1.5
    # # where whisker is greater than max or min set as max or min
    whisker_75 =  np.min((df.max(), whisker_75))
    whisker_25 = np.max((df.min(), whisker_25))
    # add column to identify fliers
    df = pd.DataFrame(df).assign(flier= False)
    df.loc[(df.iloc[:,0]<whisker_25)|(df.iloc[:,0]>whisker_75), 'flier']=True
    return(df)


def hob_fit_scatter(hobout):
    fig, ax = plt.subplots(dpi=200)
    
    # get boundary values for plotting a 1:1
    hobmax = hobout.loc[:,['sim_val','obs_val']].max().min()
    hobmin = hobout.loc[:,['sim_val','obs_val']].min().max()
    
    hob_lin = np.array([hobmin, hobmax])
    
    # plot observed vs simulated values
    hobout.plot.scatter(x='obs_val', y='sim_val',marker='.',ax=ax)
    ax.plot(hob_lin,hob_lin,'black', linestyle='--', label='1:1 line')
    # plot buffer lines with 30 ft, 10 ft limit
    hob_lin_adj = 10*0.3048
    ax.plot([hobmin-hob_lin_adj/2, hobmax-hob_lin_adj/2], [hobmin+hob_lin_adj/2, hobmax+hob_lin_adj/2], 'gray')
    ax.plot([hobmin+hob_lin_adj/2, hobmax+hob_lin_adj/2], [hobmin-hob_lin_adj/2, hobmax-hob_lin_adj/2], 'gray', label='10 ft limit')

    hob_lin_adj = 30*0.3048
    ax.plot([hobmin-hob_lin_adj/2, hobmax-hob_lin_adj/2], [hobmin+hob_lin_adj/2, hobmax+hob_lin_adj/2], 'red')
    ax.plot([hobmin+hob_lin_adj/2, hobmax+hob_lin_adj/2], [hobmin-hob_lin_adj/2, hobmax-hob_lin_adj/2], 'red', label='30 ft limit')
    
    ax.set_xlabel('Observed Values (m)')
    ax.set_ylabel('Simulated Values (m)')
    plt.legend()
    # lim2 = hobout.loc[:,['obs_val']].max().min()
    # lim1 = hobout.loc[:,['obs_val']].min().max()
    # ax.set_ylim(lim1,lim2)
    
