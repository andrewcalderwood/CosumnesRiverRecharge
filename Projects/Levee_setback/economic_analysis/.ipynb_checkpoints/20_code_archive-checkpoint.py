

# # Hydrology data
#

# +
# # load in pre-processed array of ETc for all time
# ETc_long = pd.read_hdf(join(uzf_dir, "dwr_ETc",'long_ETc_all_lu.hdf5'), key='variable')
# # identify crop to use
# etc_var = pd.Series(ETc_long.variable.unique())
# # etc_var[etc_var.str.contains(crop, case=False)]
# # crop_dict = {'Corn':'Corn (field & sweet)'}

# +
# etc_var[etc_var.str.contains('Grain', case=False)]
# # ETc_long[['Grain and hay crops']]
# ETc_long[ETc_long.variable=='Grain and hay crops']['value'][strt_date:end_date].plot()

# +
# ETc_long.variable[ETc_long.variable.str.contains('vineyards',case=False)]
# -

# While the database of Kc will work fine for the rest of the crops it may make sense to use Yusuke's specific Kc values.
# - with the updated version Yusuke has categorized the majority of crops appropriately now so the final mix may be grouped into an others class which would require a uniform Kc to be applied to align with cost estimates that are averagd.

# +
# # subset for crop
# ETc = ETc_long[ETc_long.variable==crop_dict[crop]]['value']
# # subset for model period
# ETc = ETc[strt_date:end_date]
# # fill in empty dates with interpolation (temporary since corn dates are slightly different)
# ETc = ETc.reindex(dates)
# ETc = ETc.bfill().ffill()

# # convert to an array for calculations
# ETc = ETc.values
# -

# # SWB testing

# +

    
def swb(irr_lvl):
    global wc, pc, rp, ETa, D, K_S
#     global c_gwtot, c_swtot
    m2_ac = (1/0.3048**2)/43560 # convert from m2 to acres
    in_2_m = (1/12)*0.3048 # convert from inches to meters
    nper = (end_date-strt_date).days +1

    tic = time.time()

    irr_sw = np.zeros((nper,nfield))
    irr_gw = np.zeros((nper,nfield))
    for i in np.arange(0,n_irr):
        irr_sw[irr_days[i]] = irr_lvl[i]
        irr_gw[irr_days[i]] = irr_lvl[i+n_irr]
        
    wc = np.zeros((nper+1, nfield)) # water content, add initial conditions with +1
    pc = np.zeros((nper, nfield)) # percolation
    rp = np.zeros((nper, nfield)) # runoff 
    ETa = np.zeros((nper, nfield)) # actual ET
    wb_sum= np.zeros((nper, nfield)) # water budget check
    # time units are days for everything

    D = np.zeros((nper+1, nfield)) # soil depletion, add initial conditions with +1
    K_S = np.zeros((nper, nfield)) # crop water stress
    
    # initial water content and root zone depletion are pulled from the last step of the previous run
    
    # -1 starts at IC for BC
    # WC/D starts at 0
    for ns, n in enumerate(np.arange(-1, nper-1)):
        ## Runoff ##
        S = calc_S(wc[ns+1], Smax, wc_f, soil_por)
        water_in = rain[n+1] 
        # calculate runoff only when there is rain, and rain is uniform
        if (water_in>0).any():
            rp[n+1] = ((water_in - 0.2*S)**2)/(water_in + 0.8*S)
        # where rainfall is less than initial abstraction (0.2S) set runoff as 0
        rp[n+1] = np.where(water_in<0.2*S, 0, rp[n+1])
        # add in irrigation after runoff (assume farm is set up to avoid runoff for irrigation season)
        water_in = water_in + irr_sw[n+1] + irr_gw[n+1]
        ## explicit percolation ##
        pc[n+1] = calc_pc(wc[ns], soil_por, soil_Ks, soil_m)
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
        K_S_ws = (taw - D[ns+1])/((1 - P[n+1])*taw);
        K_S[n+1] = np.where(D[ns+1]>raw[n+1], K_S_ws, K_S[n+1])

    ## Calculate daily crop outcomes 
    ETc_adj = np.transpose(K_S)*ETc # Calculate daily crop ET with soil water stress, pairwise functoin check?
    ## Calculate economic outcomes 
    arr1 = np.ones((nfield,nper))
    # average the yield scaling across the season
    # yield max changes during the season so update this to be fluid for alfalfa (mean for eaching cutting then sum means)
    Y_A_arr = np.zeros(len(y_max))
    for n in np.arange(0,len(yield_ind)-1):
        # subset the yield scaling by the growing period then multiply by the appropriate yield max for that period
        Y_A_arr[n] = y_max[n]*np.mean((arr1- (arr1 - (ETc_adj/ETc))*K_Y)[:, yield_ind[n]:yield_ind[n+1]])
    # the total yield is the sum of individual yields
    Y_A = Y_A_arr.sum()
    # # irrigation must be converted from meters to inches for calculating cost
    c_gwtot = p_e*phi*(np.multiply(dtw_arr, irr_gw[:,0])/in_2_m) # Calcualte total groundwater pumping costs for the season ($/acre)
    c_swtot = np.multiply(p_sw, irr_sw[:,0])/in_2_m # Calcualte total surface water costs for the season ($/acre)
    cost = c_gwtot+c_swtot
    # calculate profit (daily values must be summed for the seasonal value)
    # profit can be left as $/acre with area assigned after if each field is done separate
    # for alfalfa it would be helpful to see the profit by cutting
    pi = -((p_c*Y_A - np.sum(cost))- p_o) # Calculate profit ($/acre)
    # forced internal boundary to prevent negatives
    if any(irr_lvl <0):
        # set a scalable penalty
        pi = irr_lvl[irr_lvl<0].sum()*-1E4
    # # profit simplified to a function
    toc = time.time()
#     print('Run time was %.2f minutes' %((toc-tic)/60))
    if wb_sum.sum(axis=1).mean() > 1E-6:
        print('Avg WB error was %.2E m' % wb_sum.sum(axis=(1)).mean())
    return(pi)



# -

# If Ks is never less than 1 then there is never water stress (D < RAW). Water stress does appear when lower irrigation levels are used.
#
# I should double check transpose is doing what I expect with the first two cells

# We had discussed that we could adjust the irrigation volumes to assign only SW or GW irrigation if the price of GW or SW is cheaper. The issue is that the price of GW will fluctuate throughout the season so at the beginning it might be cheaper but get more expensive by the end so it could reasonably put GW at the beginning and SW at the end. To do this we need to specify which irrigation events are more or less cost.

# +
## error checking
c_gwtot = p_e*phi*(np.multiply(dtw_arr, irr_gw[:,0])/in_2_m) # Calcualte total groundwater pumping costs for the season ($/acre)
c_swtot = np.multiply(p_sw, irr_sw[:,0])/in_2_m # Calcualte total surface water costs for the season ($/acre)
cost = c_gwtot+c_swtot
print('GW Cost %.2f' %c_gwtot.sum(),'SW Cost %.2f' %c_swtot.sum())

# # the yield scale is negative sometimes
# ETc_adj = np.transpose(K_S)*ETc
# arr1 = np.ones((nfield,nper))
# # average the yield scaling across the season
# Y_A_arr = np.zeros(len(y_max))
# for n in np.arange(0,len(yield_ind)-1):
#     # subset the yield scaling by the growing period then multiply by the appropriate yield max for that period
#     Y_A_arr[n] = y_max[n]*np.mean((arr1- (arr1 - (ETc_adj/ETc))*K_Y)[:, yield_ind[n]:yield_ind[n+1]])
# # the total yield is the sum of individual yields
# Y_A = Y_A_arr.sum()
Y_A = calc_yield(ETc, K_S, y_max, yield_ind,  nfield, nper)
p_c*Y_A - p_o - np.sum(cost)
p_c*Y_A, np.sum(cost), p_o

# -



# profit for alfalfa (1,000) seems a little excessive but operating costs are low and if 7 cuttings are done then it's a lot
# vineyard profit is even greater (2,000)
# corn only nets around 200 
#
# Misc. grain and hay didn't find a positive solution (started -397 went -307) -> issue is the Kc values used were for summer crop not winter so there was no ETc to make yield

# # Solver testing

# +

def mak_irr_con(soil_ag):
    ## for no POD case the SW limit would be 0
    sw_scale = 1
    gw_scale = 1
    if soil_ag.pod.iloc[0]=='No Point of Diversion on Parcel':
        sw_scale = 0
        gw_scale = 2 # give groundwater twice as much availability

    # Total surface water and groundwater available during the season (in)
    # irr_tot = np.array([10*sw_scale, 10*gw_scale]) 
    # to represent unconstrained conditions set boundaries at 100, 100
    irr_tot = np.array([100*sw_scale, 100*gw_scale]) 
    irr_tot = (irr_tot/12)*0.3048 # convert to meters
    # Coefficients for inequality constraints (first n_irr columns are for surface water; second n_irr columns are for groundwater)
    ACON = np.zeros((2,2*n_irr))
    ACON[0,:n_irr] = np.ones(n_irr)
    ACON[1,(n_irr):(2*n_irr)] = np.ones(n_irr)
    con_min = np.zeros(len(ACON)) 

    # make constraint
    linear_constraint = LinearConstraint(ACON, list(con_min), list(irr_tot))
    return linear_constraint


# AEQCON = [] #  No equality contraints
# BEQCON = [] # No equality contraints
# I_WMAXLBCON = np.zeros((1,2*n_irr)) # Irrigation cannot be negative

# COBYLA can't use bounds so need extra inequality constraints
# pos_con = np.zeros((len(irr_lvl), len(irr_lvl)))
# np.fill_diagonal(pos_con,1)
# ACON = np.append(ACON, pos_con, axis=0)
# irr_tot = np.append(irr_tot, np.full(len(pos_con), 10))

bounds = Bounds(lb = 0)

linear_constraint = mak_irr_con(soil_ag)

# +
# # for SLSQP
# ineq_cons = {'type':'ineq',
#             'fun': lambda x: np.array([irr_tot[0] - np.sum(x[:n_irr]),
#                                 irr_tot[1] - np.sum(x[n_irr:])]),
#              'jac': lambda x: np.array([[10],
#                                       [10]])
#             }

# +
# the minimization with 'trust-constr' and no constraints doesn't solve and has increasing WB error
# out = minimize(swb, irr_lvl, method='SLSQP',
# #         constraints = [ineq_cons],
#         bounds=bounds,
# #          options={'verbose':1}
#         )
# -



# +
Nelder-Mead isn't faster, fails to solve with 300 maxiter (can be bounded)  
CG - no bounds available, tried to set internal boundary but model accepted the penalty or after adjusting internal constraint it failed to solve  
-> forcing an internal constraint only serves to create more non-linearity rather than setting the solve space  

SLSQP was able to solve the problem with an internal 0-bound but it found a solution with primarily GW and a little bit of SW which doesn't make sense since GW is cheaper. This was not improved with finer tolerance (1E-2 times more).  
COBYLA was able to evaluate successfully after increasing maxiter to 600. Same issue as SLSQP where a mix of GW and SW are used despite no constraints.  

# +
# out = minimize(swb, irr_lvl, method='COBYLA',
# #         constraints = [linear_constraint],
# #         bounds=bounds,
#          options={'disp':1, 'maxiter':600},
#                tol=0.0001
#         )
# -

Notes on solvers from:
- fmincon from MATLAB was suggested anecodotally as faster than scipy.minimize
- someone suggested CVXPY as an alternate which is what Jon Herman had us use for non-linear convex problems

# 'trust-constr', 
#'SLSQP' and 'COBYLA' -> don't seem to work easily



