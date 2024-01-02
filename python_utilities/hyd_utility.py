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
