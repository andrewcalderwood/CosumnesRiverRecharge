"""
report_cln module. 
Different functions for refining numeric/tabular data with general python functions
First iteration as a Module May 2023
Author: Andrew Calderwood
"""

import numpy as np
import matplotlib.pyplot as plt

def dt_2_wy(dt):
    dt = pd.Series(dt)
    yr = dt.dt.year
    wy = np.where(dt.dt.month>=10, yr+1, yr)
    return(wy)

def base_round(x, base=1):
    """ Round to different numeric bases (e.g.,5, 20)
    Ex: base_round(11, 5) = 10 or base_round(14,5) = 15
    """
    return base * round(x/base)

def magnitude(x):
    """ Return the order of magnitude of a number regardless of sign
    """
    return int(np.log10(np.abs(x)))

def fmt(x):
    """ format numeric contour labels to avoid excess decimals
    """
    # print with 1 decimal
    s = f"{x:.1f}"
    # if the decimal is a zero then print with no decimal
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} " if plt.rcParams["text.usetex"] else f"{s} "

def nse(targets,predictions):
    """ Calculate the Nash-Sutcliffe efficiency for target and observed data
    """
    return 1-(np.sum((targets-predictions)**2)/np.sum((targets-np.mean(predictions))**2))
