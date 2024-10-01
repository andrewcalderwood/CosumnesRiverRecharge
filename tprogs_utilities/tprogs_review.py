"""
tprogs_cleaning module. 
Different functions for preparing data for use in MODFLOW and other groundwater modeling software.
First iteration as a Module Oct 2021
Author: Andrew Calderwood
"""

import sys
import numpy as np
import pandas as pd
from numpy import ma

def get_tprogs_quants(vka, params, tprogs_vals, tprogs_hist):
    """
    Find the upper and lower values of hydraulic parameters in 
    a 3D array format based on the histogram values 
    """
    tprogs_quants = 1-np.append([0], np.cumsum(tprogs_hist)/np.sum(tprogs_hist))
    vka_quants = pd.DataFrame(tprogs_quants[1:], columns=['quant'], index=tprogs_vals)
    # dataframe summarizing dominant facies based on quantiles
    vka_quants['vka_min'] = np.quantile(vka, tprogs_quants[1:])
    vka_quants['vka_max'] = np.quantile(vka, tprogs_quants[:-1])
    vka_quants['facies'] = params.loc[tprogs_vals].Lithology.values
    return(vka_quants)