"""
swb_utility module. 
Different functions for soil water budget set up and post-processing with general python functions
First iteration as a Module August 2023
Author: Andrew Calderwood
"""

import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
import pandas as pd

import h5py

# not sure how to properly apply global variables
# def prep_soil(soil_ag):
# 	""" Given a dataframe with SSURGO data summarized to the field or cell level
# 	return the prepared values for use in the soil water budget calculations"""
# 	global soil_Ks, soil_por, soil_eps, soil_CN
# 	global soildepth, soil_m, wc_f, wc_wp, Smax
# 	global irr_eff_mult
# 	# # when soil_K_low is missing using a substitute of Ksat/10
# 	soil_Ks = np.where(soil_ag.Ksat_Low==0, soil_ag.Ksat/10, soil_ag.Ksat_Low)
# 	soil_por = soil_ag.Porosity.values/100
# 	soil_eps = soil_ag.EPS.values
# 	soil_CN = soil_ag.CN.values
# 	#     irr_eff = soil_ag.Avg_eff.values/100
# 	irr_eff_mult = soil_ag.irr_eff_mult.values

# 	soildepth = soil_ag.SoilDepth.values
# 	psdi =  soil_ag.PSDI.values
# 	# parameter for Mualem, van Genuchten
# 	soil_m = psdi/(psdi+1)
# 	wc_f =  soil_ag.w3rdbar.values/100 #field capacity
# 	wc_wp =  soil_ag.w15bar.values/100 #wilting point 

# 	# Calculate total available water in the root zone
# 	taw = (wc_f - wc_wp)*soildepth 

# 	# for runoff, convert CN from fraction to CN
# 	Smax = (1000/soil_CN) - 10
# 	return()



def calc_S(wc, Smax, wc_f, soil_por):
    """ Given an array of water contents return potential soil retention
    wc: water content, Smax: max soil retention, wc_f: field water content, soil_por: soil porosity"""
    S = np.zeros(Smax.shape)
    # when water content is less than 1/2 field capacity, full retention
    S[wc < wc_f/2] = Smax[wc < wc_f/2]
    # wc > porosity then set as porosity for calculations (avoids S >1)
    wc_calc = np.where(wc<soil_por, wc, soil_por)
    # when water content is greater than 1/2 field capacity, partial retention 
    S1 = Smax * (1 - ((wc_calc - wc_f/2)/(soil_por - wc_f/2)))
    S[wc >= wc_f/2]= S1[wc >= wc_f/2]
    # convert S from inches to meters
    S *= (0.3048/12)
    return(S)

def calc_pc(wc, soil_por, soil_Ks, soil_m):
    """ explicit calculation of percolation assuming water content of prior to 
    percolation is representative of conditions controlling percolation
    wc: water content, Smax: max soil retention, wc_f: field water content, soil_por: soil porosity"""
    # calculate soil saturation as water content can be greater than porosity assuming it represents ponding
    sat = wc/soil_por
    sat  = np.where(sat>1, 1, sat)
    # explicit calculation of percolation
    pc = soil_Ks*(sat)*(1- (1-(sat)**(1/soil_m))**soil_m)**2
    return(pc)

def load_swb_data(strt_date, end_date, param):
    """Returns data in an array format with the first dimension for time
    and then either a 1D or 2D array representing fields or model grid cells"""
    nper_tr = (end_date-strt_date).days+1
#     nrow_p,ncol_p = (100,230)
#     param_arr = np.zeros((nper_tr, nrow_p,ncol_p))
    # years and array index, include extra year in index to fully index years
    years = pd.date_range(strt_date,end_date+pd.DateOffset(years=1),freq='AS-Oct')
    yr_ind = (years-strt_date).days
    years= years[:-1]
    # need separte hdf5 for each year because total is 300MB
    for n in np.arange(0,len(yr_ind)-1):
        fn = join(uzf_dir, 'basic_soil_budget',param+"_WY"+str(years[n].year+1)+".hdf5")
        with h5py.File(fn, "r") as f:
            # load first dataset then allocate array
            if n == 0:
                dim  = list(f['array']['WY'].shape)
                dim[0] = nper_tr
                param_arr = np.zeros(dim)
            arr = f['array']['WY'][:]
            param_arr[yr_ind[n]:yr_ind[n+1]] = arr
    return(param_arr)