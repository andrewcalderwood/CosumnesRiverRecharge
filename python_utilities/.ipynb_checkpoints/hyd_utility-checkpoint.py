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