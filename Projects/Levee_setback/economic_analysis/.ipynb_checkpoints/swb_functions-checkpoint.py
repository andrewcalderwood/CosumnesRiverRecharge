# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import numpy as np
import pandas as pd

# %%


# %%

# %%
# what if this code was moved to the pre-processing
# then the cleaned soil data could be avialable up front and for post-processing
def prep_soil_dict(soil_ag, etc_arr, var_crops):
    """Given a dataframe of soil properties create a dictionary for reference
    """
    # dictionary to be filled with soil data for Soil Water Budget
    soil = dict()
    soil['nfield'] = len(soil_ag)

    # # when soil_K_low is missing using a substitute of Ksat/10
    soil['Ks'] = np.where(soil_ag.Ksat_Low==0, soil_ag.Ksat/10, soil_ag.Ksat_Low)
    soil['por'] = soil_ag.Porosity.values/100
    soil['eps'] = soil_ag.EPS.values
    soil['CN'] = soil_ag.CN.values
    soil['psdi'] =  soil_ag.PSDI.values
    # parameter for Mualem, van Genuchten
    soil['m'] = soil['psdi']/(soil['psdi']+1)
    soil['wc_f'] =  soil_ag.w3rdbar.values/100 #field capacity
    soil['wc_wp'] =  soil_ag.w15bar.values/100 #wilting point 

    # calculate the minimum of soil/crop rooting depth
    zr = (var_crops.zr/12)*0.3048
    soil['depth'] = np.min([soil_ag.SoilDepth.values, np.repeat(zr,len(soil_ag))], axis=0)
    # soil['depth'] = soil_ag.SoilDepth.values
    
    # Calculate total available water in the root zone
    soil['taw'] = (soil['wc_f'] - soil['wc_wp'])*soil['depth']

    # for runoff, convert CN from fraction to CN
    soil['Smax'] = (1000/soil['CN']) - 10
    
    # p_table22: Soil water depletion fraction for no stress, requires ETc in mm/day
    soil['P'] = var_crops['p_table22'] + 0.04*((5-(etc_arr*1000))) # Calculate adjusted daily soil water depletion fraction for no stress
    soil['raw'] = soil['taw']*soil['P'] # Calculate readily available water in the root zone
    return(soil)


# %%
def calc_S(wc, Smax, wc_f, soil_por):
    """ Given an array of water contents return potential soil retention"""
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


# %%
def calc_pc(wc, soil_por, soil_Ks, soil_m):
    """ explicit calculation of percolation assuming water content of prior to 
    percolation is representative of conditions controlling percolation"""
    # calculate soil saturation as water content can be greater than porosity assuming it represents ponding
    sat = wc/soil_por
    sat  = np.where(sat>1, 1, sat)
    # explicit calculation of percolation
    pc = soil_Ks*(sat)*(1- (1-(sat)**(1/soil_m))**soil_m)**2
    return(pc)

# calc_pc(wc[100], soil_por, soil_Ks, soil_m)


# %%
def calc_yield(ETc, K_S, gen):
    """
    ETc : crop evapotranspiration
    K_S : crop water stress
    gen : dictionary with yield variables

    """
    # set up local variables
    K_Y = gen.K_Y
    y_max = gen.y_max
    yield_ind = gen.yield_ind
    
    ## Calculate daily crop outcomes 
    ETc_adj = np.transpose(K_S)*ETc # Calculate daily crop ET with soil water stress, pairwise functoin check?
    ## Calculate economic outcomes 
    arr1 = np.ones(K_S.shape)
    # average the yield scaling across the season
    # yield max changes during the season so update this to be fluid for alfalfa (mean for eaching cutting then sum means)
    Y_A_arr = np.zeros(len(y_max))
    for n in np.arange(0,len(yield_ind)-1):
        # subset the yield scaling by the growing period then multiply by the appropriate yield max for that period
        Y_A_arr[n] = y_max[n]*np.mean((arr1- (arr1 - (ETc_adj/ETc))*K_Y)[:, yield_ind[n]:yield_ind[n+1]])
    # the total yield is the sum of individual yields
    Y_A = Y_A_arr.sum()
    # can't sum yield if the crop price varies by season
    Y_A = np.copy(Y_A_arr)
    return(Y_A)


# %%

       
def calc_profit(Y_A, dtw_arr, irr_gw, irr_sw, gen):
    """
    Y_A : actual yield (tons)
    dtw_arr : depth to water (ft)
    irr_gw : groundwater irrigation (m)
    irr_sw : surface water irrigation (m)
    gen : dictionary with cost variables
    """
    # set up local variables
    p_c = gen.p_c
    p_e = gen.p_e # energy price, $/kWh
    p_o = gen.p_o # operating cost, $/acre
    p_sw = gen.p_sw # surface water cost, $/acre-in
    phi = gen.phi # energy req to raise unit of water per unit length, kWh/acre-in/ft

    in_2_m = (1/12)*0.3048 # convert from inches to meters
    c_gwtot = p_e*phi*(np.multiply(dtw_arr, irr_gw[:,0])/in_2_m) # Calculate total groundwater pumping costs for the season ($/acre)
    c_swtot = np.multiply(p_sw, irr_sw[:,0])/in_2_m # Calcualte total surface water costs for the season ($/acre)
    cost = c_gwtot+c_swtot
    # calculate profit (daily values must be summed for the seasonal value)
    pi = -((np.sum(p_c*Y_A - p_o) - np.sum(cost))) # Calculate profit ($/acre)
    # forced internal boundary to prevent negatives
    # if any(irr_lvl <0):
    #     # set a scalable penalty, assuming p_o would be a sizable penalty
    #     pi = irr_lvl[irr_lvl<0].sum()*-p_o*10
    return(pi)


# %%
def choose_water_source(dtw_arr, gen, mix_fraction = 1):
    """
    Determine if GW or SW is more efficient
    dtw_arr : depth to water (ft)
    gen : dictionary with cost variables
    """
    # set up local variables
    p_e = gen.p_e # energy price, $/kWh
    p_sw = gen.p_sw # surface water cost, $/acre-in
    phi = gen.phi # energy req to raise unit of water per unit length, kWh/acre-in/ft

    # use 1 inch of irrigation for the cost
    p_gw = p_e*phi*dtw_arr.mean() # Calculate total groundwater pumping costs for the season ($/acre)
    # return string that specifies water source
    # apply a mix fraction to require a minimum offset
    if p_gw <= p_sw*mix_fraction:
        water_source = 'gw'
    elif p_sw < p_gw*mix_fraction:
        water_source = 'sw'
    else:
        water_source = 'mixed'
    
    return(water_source)

# %%

    
def run_swb_model(soil, rain, ETc, irr_gw, irr_sw):
    """
    irr_lvl: depth of irrigation to apply
    soil: class object with soil parameter for the given field
    gen: class object with the general parameters for the irrigation cost model
    rain: 1D array of the rainfall depth for each step
    ETc: 1D array of the crop evapotranspiration depth for each step

    arrays: boolean to identify whether output arrays should be returned
        True: water budget arrays are the output
        False: the output is the net profit (used for optimization)
    """
    nfield = soil.nfield
    nper = ETc.shape[0]
    
    m2_ac = (1/0.3048**2)/43560 # convert from m2 to acres
    in_2_m = (1/12)*0.3048 # convert from inches to meters

    wc = np.zeros((nper+1, nfield)) # water content, add initial conditions with +1
    pc = np.zeros((nper, nfield)) # percolation
    rp = np.zeros((nper, nfield)) # runoff 
    ETa = np.zeros((nper, nfield)) # actual ET
    wb_sum= np.zeros((nper, nfield)) # water budget check
    # time units are days for everything

    D = np.zeros((nper+1, nfield)) # soil depletion, add initial conditions with +1
    K_S = np.zeros((nper, nfield)) # crop water stress

    # most frequently used variables can be locally defined
    soildepth = soil.depth
    taw = soil.taw
    
    # initial water content and root zone depletion are pulled from the last step of the previous run
    # -1 starts at IC for BC
    # WC/D starts at 0
    for ns, n in enumerate(np.arange(-1, nper-1)):
        ## Runoff ##
        S = calc_S(wc[ns+1], soil.Smax, soil.wc_f, soil.por)
        water_in = rain[n+1] 
        # calculate runoff only when there is rain, and rain is uniform
        if (water_in>0).any():
            rp[n+1] = ((water_in - 0.2*S)**2)/(water_in + 0.8*S)
        # where rainfall is less than initial abstraction (0.2S) set runoff as 0
        rp[n+1] = np.where(water_in<0.2*S, 0, rp[n+1])
        # add in irrigation after runoff (assume farm is set up to avoid runoff for irrigation season)
        water_in = water_in + irr_sw[n+1] + irr_gw[n+1]
        ## explicit percolation ##
        pc[n+1] = calc_pc(wc[ns], soil.por, soil.Ks, soil.m)
        # stepwise water budget, explicit to avoid iteration
        # add rain and take away runoff first
        wc[ns+1] = (wc[ns]*soildepth + (water_in - rp[n+1]))/soildepth
        # take away ET, add term to prevent going to zero
        ETa[n+1] = np.where(ETc[n+1] <= wc[ns+1]*soildepth, ETc[n+1], wc[ns+1]*soildepth - 1E-9)
        wc[ns+1] = wc[ns+1] + (-ETa[n+1])/soildepth
        # take away percolation
        pc[n+1] = np.where(pc[n+1] <= wc[ns+1]*soildepth, pc[n+1], wc[ns+1]*soildepth - 1E-9)
        wc[ns+1] = wc[ns+1] + (-pc[n+1])/soildepth
        # check water budget error
        wb_sum[n+1] = (wc[ns]-wc[ns+1])*soildepth + water_in - rp[n+1] - ETa[n+1] - pc[n+1] 
        if (wb_sum[n+1]>1E-3).any()|(wb_sum[n+1]<-1E-3).any():
            print('WB error exceeds 1E-3',n )
            ## additional code for optimizing irrigation
        # calculate soil depletion for irrigation decision (must use ETc to see how much should be depleted)
        D[ns+1] = D[ns] - water_in + ETc[n+1] + rp[n+1] + pc[n+1] 
        # root zone depletion can't be greater than TAW 
        D[ns+1] = np.min([D[ns+1], taw], axis=0)
        # root zone depletion should be greater than 0
        D[ns+1] = np.where(D[ns+1]<0,0, D[ns+1])
        # default value of water stress is 1 (none): # potentially unnecessary just fill in 1
        K_S[n+1] = 1
        # where rootzone depletion is greater than RAW there is water stress
        K_S_ws = (taw - D[ns+1])/((1 - soil.P[n+1])*taw);
        K_S[n+1] = np.where(D[ns+1]>soil.raw[n+1], K_S_ws, K_S[n+1])

        return(pi, pc, K_S, Y_A)




# %%
def run_swb(irr_lvl, soil, gen, rain, ETc, dtw_arr, irr_src='both', arrays = False):
    """
    dtw_arr: 1D array of the depth to groundwater for each step
    irr_src: string variable to identify the expected shape of irr_lvl
        'both': array of 2*n_irr with depths for GW and SW
        'sw': array of n_irr with depths for SW 
        'gw: array of n_irr with depths for GW
    """
    nper = gen.nper
    n_irr = gen.n_irr
    irr_days = gen.irr_days
    
    irr_sw = np.zeros((nper,nfield))
    irr_gw = np.zeros((nper,nfield))
    # updated code to correct if irr_lvl is only for one irrigation type
    if irr_src=='both':
        for i in np.arange(0,n_irr):
            irr_sw[irr_days[i]] = irr_lvl[i]
            irr_gw[irr_days[i]] = irr_lvl[i+n_irr]
    elif irr_src=='sw':
        for i in np.arange(0,n_irr):
            irr_sw[irr_days[i]] = irr_lvl[i]
    elif irr_src=='gw':
        for i in np.arange(0,n_irr):
            irr_gw[irr_days[i]] = irr_lvl[i]

    pi, pc, K_S, Y_A = run_swb_model(soil, rain, ETc, irr_gw, irr_sw)
    
    ## Calculate yield outcomes 
    Y_A = calc_yield(ETc, K_S, gen)
    # Y_A = calc_yield(ETc, K_S, K_Y, y_max, yield_ind,  nfield, nper)
    
    ## profit simplified to a function
    # pi = calc_profit(Y_A, p_c, p_e, phi, dtw_arr, irr_gw, p_sw, irr_sw, p_o)  
    pi = calc_profit(Y_A, dtw_arr, irr_gw,irr_sw, gen)  
    if wb_sum.sum(axis=1).mean() > 1E-6:
        print('Avg WB error was %.2E m' % wb_sum.sum(axis=(1)).mean())
    if arrays:
        # for secondary output need to also save deep percolation
        return(pi, pc, K_S, Y_A)
    return(pi)

# %%

    
# def run_swb(irr_lvl, soil, gen, rain, ETc, dtw_arr, irr_src='both', arrays = False):
#     """
#     irr_lvl: depth of irrigation to apply
#     soil: class object with soil parameter for the given field
#     gen: class object with the general parameters for the irrigation cost model
#     rain: 1D array of the rainfall depth for each step
#     ETc: 1D array of the crop evapotranspiration depth for each step
#     dtw_arr: 1D array of the depth to groundwater for each step
#     irr_src: string variable to identify the expected shape of irr_lvl
#         'both': array of 2*n_irr with depths for GW and SW
#         'sw': array of n_irr with depths for SW 
#         'gw: array of n_irr with depths for GW
#     arrays: boolean to identify whether output arrays should be returned
#         True: water budget arrays are the output
#         False: the output is the net profit (used for optimization)
#     """
#     nper = gen.nper
#     nfield = soil.nfield
#     n_irr = gen.n_irr
#     irr_days = gen.irr_days
#     m2_ac = (1/0.3048**2)/43560 # convert from m2 to acres
#     in_2_m = (1/12)*0.3048 # convert from inches to meters

#     irr_sw = np.zeros((nper,nfield))
#     irr_gw = np.zeros((nper,nfield))
#     # updated code to correct if irr_lvl is only for one irrigation type
#     if irr_src=='both':
#         for i in np.arange(0,n_irr):
#             irr_sw[irr_days[i]] = irr_lvl[i]
#             irr_gw[irr_days[i]] = irr_lvl[i+n_irr]
#     elif irr_src=='sw':
#         for i in np.arange(0,n_irr):
#             irr_sw[irr_days[i]] = irr_lvl[i]
#     elif irr_src=='gw':
#         for i in np.arange(0,n_irr):
#             irr_gw[irr_days[i]] = irr_lvl[i]
        
#     wc = np.zeros((nper+1, nfield)) # water content, add initial conditions with +1
#     pc = np.zeros((nper, nfield)) # percolation
#     rp = np.zeros((nper, nfield)) # runoff 
#     ETa = np.zeros((nper, nfield)) # actual ET
#     wb_sum= np.zeros((nper, nfield)) # water budget check
#     # time units are days for everything

#     D = np.zeros((nper+1, nfield)) # soil depletion, add initial conditions with +1
#     K_S = np.zeros((nper, nfield)) # crop water stress

#     # most frequently used variables can be locally defined
#     soildepth = soil.depth
#     taw = soil.taw
    
#     # initial water content and root zone depletion are pulled from the last step of the previous run
#     # -1 starts at IC for BC
#     # WC/D starts at 0
#     for ns, n in enumerate(np.arange(-1, nper-1)):
#         ## Runoff ##
#         S = calc_S(wc[ns+1], soil.Smax, soil.wc_f, soil.por)
#         water_in = rain[n+1] 
#         # calculate runoff only when there is rain, and rain is uniform
#         if (water_in>0).any():
#             rp[n+1] = ((water_in - 0.2*S)**2)/(water_in + 0.8*S)
#         # where rainfall is less than initial abstraction (0.2S) set runoff as 0
#         rp[n+1] = np.where(water_in<0.2*S, 0, rp[n+1])
#         # add in irrigation after runoff (assume farm is set up to avoid runoff for irrigation season)
#         water_in = water_in + irr_sw[n+1] + irr_gw[n+1]
#         ## explicit percolation ##
#         pc[n+1] = calc_pc(wc[ns], soil.por, soil.Ks, soil.m)
#         # stepwise water budget, explicit to avoid iteration
#         # add rain and take away runoff first
#         wc[ns+1] = (wc[ns]*soildepth + (water_in - rp[n+1]))/soildepth
#         # take away ET, add term to prevent going to zero
#         ETa[n+1] = np.where(ETc[n+1] <= wc[ns+1]*soildepth, ETc[n+1], wc[ns+1]*soildepth - 1E-9)
#         wc[ns+1] = wc[ns+1] + (-ETa[n+1])/soildepth
#         # take away percolation
#         pc[n+1] = np.where(pc[n+1] <= wc[ns+1]*soildepth, pc[n+1], wc[ns+1]*soildepth - 1E-9)
#         wc[ns+1] = wc[ns+1] + (-pc[n+1])/soildepth
#         # check water budget error
#         wb_sum[n+1] = (wc[ns]-wc[ns+1])*soildepth + water_in - rp[n+1] - ETa[n+1] - pc[n+1] 
#         if (wb_sum[n+1]>1E-3).any()|(wb_sum[n+1]<-1E-3).any():
#             print('WB error exceeds 1E-3',n )
#             ## additional code for optimizing irrigation
#         # calculate soil depletion for irrigation decision (must use ETc to see how much should be depleted)
#         D[ns+1] = D[ns] - water_in + ETc[n+1] + rp[n+1] + pc[n+1] 
#         # root zone depletion can't be greater than TAW 
#         D[ns+1] = np.min([D[ns+1], taw], axis=0)
#         # root zone depletion should be greater than 0
#         D[ns+1] = np.where(D[ns+1]<0,0, D[ns+1])
#         # default value of water stress is 1 (none): # potentially unnecessary just fill in 1
#         K_S[n+1] = 1
#         # where rootzone depletion is greater than RAW there is water stress
#         K_S_ws = (taw - D[ns+1])/((1 - soil.P[n+1])*taw);
#         K_S[n+1] = np.where(D[ns+1]>soil.raw[n+1], K_S_ws, K_S[n+1])

#     ## Calculate yield outcomes 
#     Y_A = calc_yield(ETc, K_S, gen)
#     # Y_A = calc_yield(ETc, K_S, K_Y, y_max, yield_ind,  nfield, nper)
    
#     ## profit simplified to a function
#     # pi = calc_profit(Y_A, p_c, p_e, phi, dtw_arr, irr_gw, p_sw, irr_sw, p_o)  
#     pi = calc_profit(Y_A, dtw_arr, irr_gw,irr_sw, gen)  
#     if wb_sum.sum(axis=1).mean() > 1E-6:
#         print('Avg WB error was %.2E m' % wb_sum.sum(axis=(1)).mean())

#     if arrays:
#         # for secondary output need to also save deep percolation
#         return(pi, pc, K_S, Y_A)
#     return(pi)

# %%
from scipy.optimize import Bounds, LinearConstraint


# %%

def mak_irr_con(n_irr, sw_con = 100, gw_con = 100):
    """ 
    Make simple constraints on SW and GW with seasonal limits (inches)
    The unconstrained version has very high limits (unreachable)
    """

    # Total surface water and groundwater available during the season (in)
    irr_tot = np.array([sw_con, gw_con]) 
    irr_tot = (irr_tot/12)*0.3048 # convert to meters
    # Coefficients for inequality constraints (first n_irr columns are for surface water; second n_irr columns are for groundwater)
    ACON = np.zeros((2,2*n_irr))
    ACON[0,:n_irr] = np.ones(n_irr)
    ACON[1,(n_irr):(2*n_irr)] = np.ones(n_irr)
    con_min = np.zeros(len(ACON)) 

    # make constraint
    linear_constraint = LinearConstraint(ACON, list(con_min), list(irr_tot))
    return linear_constraint


# %%
# want to make a version of constraints that adapt the number of parameters based on hard constraints

def mak_irr_con_adj(n_irr, sw_con = 100, gw_con = 100):
    """ 
    Make simple constraints on SW and GW with seasonal limits (inches)
    The unconstrained version has very high limits (unreachable)
    """

    # Total surface water and groundwater available during the season (in)
    irr_tot = np.array([sw_con, gw_con]) 
    # if one constraint is zero then remove from calculation
    if sw_con==0:
        irr_tot = np.array([gw_con])
    elif gw_con==0:
        irr_tot = np.array([sw_con])
    irr_tot = (irr_tot/12)*0.3048 # convert to meters
    # Coefficients for inequality constraints (first n_irr columns are for surface water; second n_irr columns are for groundwater)
    n_wt = len(irr_tot)
    ACON = np.zeros((n_wt,n_wt*n_irr))
    ACON[0,:n_irr] = np.ones(n_irr)
    if n_wt==2:
        ACON[1,(n_irr):(n_wt*n_irr)] = np.ones(n_irr)

    con_min = np.zeros(len(ACON)) 

    # make constraint
    linear_constraint = LinearConstraint(ACON, list(con_min), list(irr_tot))
    return linear_constraint


# %%

def mak_irr_con_old(soil_ag, n_irr, sw_con = 100, gw_con = 100):
    """ 
    Make simple constraints on SW and GW with seasonal limits (inches)
    The unconstrained version has very high limits (unreachable)
    """
    ## for no POD case the SW limit would be 0
    sw_scale = 1
    gw_scale = 1
    if soil_ag.pod.iloc[0]=='No Point of Diversion on Parcel'|sw_con==0:
        sw_scale = 0
        gw_scale = 2 # give groundwater twice as much availability

    # Total surface water and groundwater available during the season (in)
    irr_tot = np.array([sw_con*sw_scale, gw_con*gw_scale]) 
    irr_tot = (irr_tot/12)*0.3048 # convert to meters
    # Coefficients for inequality constraints (first n_irr columns are for surface water; second n_irr columns are for groundwater)
    ACON = np.zeros((2,2*n_irr))
    ACON[0,:n_irr] = np.ones(n_irr)
    ACON[1,(n_irr):(2*n_irr)] = np.ones(n_irr)
    con_min = np.zeros(len(ACON)) 

    # make constraint
    linear_constraint = LinearConstraint(ACON, list(con_min), list(irr_tot))
    return linear_constraint
