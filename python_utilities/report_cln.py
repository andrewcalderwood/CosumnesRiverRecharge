"""
report_cln module. 
Different functions for refining numeric/tabular data with general python functions
First iteration as a Module May 2023
Author: Andrew Calderwood
"""

import numpy as np
import matplotlib.pyplot as plt

def base_round(x, base=1):
    """ Round to different numeric bases (e.g.,5, 20)
    Ex: base_round(11, 5) = 10 or base_round(14,5) = 15
    """
    return base * round(x/base)


def fmt(x):
    """ format numeric contour labels to avoid excess decimals
    """
    # print with 1 decimal
    s = f"{x:.1f}"
    # if the decimal is a zero then print with no decimal
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} " if plt.rcParams["text.usetex"] else f"{s} "