# standard python utilities
import pandas as pd
import numpy as np
import os 

# run installed version of flopy or add local path
try:
    import flopy
    from flopy.discretization.structuredgrid import StructuredGrid
    from flopy.utils.reference import SpatialReference
except:
    import flopy
    fpth = os.path.abspath(os.path.join('..', '..'))
    sys.path.append(fpth)
    from flopy.discretization.structuredgrid import StructuredGrid
    from flopy.utils.reference import SpatialReference


# save modflow workspace file to WRDAPP sub folder to improve runtime calculations
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes_simple'
# load model with only DIS to reduce load time
# the model will run off of the .nam file connection so flopy doesn't need them
all_model_ws = loadpth + '/WEL_RIV_RCH_layercake'
model_ws = os.getcwd()
m = flopy.modflow.Modflow.load('MF.nam', model_ws=model_ws, 
                                exe_name='mf2005', version='mf2005')

###############################################################################
## LPF Package ##
deep_geology = np.loadtxt(all_model_ws+'/input_data/deep_geology.tsv', delimiter ='\t')
deep_geology = np.reshape(deep_geology, (m.dis.nlay,m.dis.nrow,m.dis.ncol))
# initial guess for hydraulic parameters
params = pd.read_csv('ZonePropertiesInitial.csv')
params = params.set_index('Zone')
# convert from m/s to m/d
params['K_m_d'] = params.K_m_s * 86400    

botm = m.dis.botm.array

# allocate arrays for Kx, Ky
hk = np.zeros(botm.shape)
vka = np.zeros(botm.shape)
sy = np.zeros(botm.shape)
ss = np.zeros(botm.shape)

# set value for unsaturated zone/layer 1
hk[0,:,:] = np.loadtxt(all_model_ws+'/input_data/unsat_hk.tsv',delimiter = '\t')
vka[0,:,:] = np.loadtxt(all_model_ws+'/input_data/unsat_vka.tsv',delimiter = '\t')
sy[0,:,:] = np.loadtxt(all_model_ws+'/input_data/unsat_sy.tsv',delimiter = '\t')
ss[0,:,:] = np.loadtxt(all_model_ws+'/input_data/unsat_s.tsv', delimiter = '\t')

# set values for second to bottom layer, Laguna formation
hk[-2,:,:] = params.loc[5,'K_m_d']
vka[-2,:,:] = params.loc[5,'K_m_d']/100 # assume 1/100 for Kx to Kz
sy[-2,:,:] = params.loc[5,'Sy']
ss[-2,:,:] = params.loc[5,'Ss']

# the deep_geology array shows where the mehrten formation comes out of the surface
hk[deep_geology[:,:,:].astype(bool)] = params.loc[6,'K_m_d']
vka[deep_geology[:,:,:].astype(bool)] = params.loc[6,'K_m_d']/100 # assume 1/100 for Kx to Kz
sy[deep_geology[:,:,:].astype(bool)] = params.loc[6,'Sy']
ss[deep_geology[:,:,:].astype(bool)] = params.loc[6,'Ss']

# set values for bottom layer, Mehrten formation
hk[-1,:,:] = params.loc[6,'K_m_d']
vka[-1,:,:] = params.loc[6,'K_m_d']/100 # assume 1/100 for Kx to Kz
sy[-1,:,:] = params.loc[6,'Sy']
ss[-1,:,:] = params.loc[6,'Ss']

# 0—indicates VKA is vertical hydraulic conductivity
layvka = 0
# no defined anisotropy between kx and ky, could set a value based on stream deposition
hani = 1 

# LAYTYP MUST BE GREATER THAN ZERO WHEN IUZFOPT IS 2
# 0 is confined, >0 convertible, <0 convertible unless the THICKSTRT option is in effect
# laytyp = [1,1,1,0,0,0,0,0]
laytyp = np.zeros(m.dis.nlay,dtype=int).tolist()

laywet=laytyp[:]
# Laywet must be 0 if laytyp is confined laywet = [1,1,1,1,1]
# laywet = 1 means layers can be rewetted.
#ipakcb = 53 means cell-by-cell budget is saved because it is non zero (default is 53)
# indicates that variable Ss and SS parameters are read as storage coefficient rather than specific storage. (default is False).
# removed sy while ss is active as storage coefficient rather than specific storage
# ss= sy for storage coefficient run
# using storage coefficient drastically increased run time from 10 to 25 minutes, storagecoefficient =True
lpf = flopy.modflow.ModflowLpf(model = m, hk =hk, chani = 0, hani = hani,
                               layvka = layvka, vka = vka, ss=ss, sy=sy,
                               laytyp=laytyp, laywet = laywet, ipakcb=53)
# overwrite the previous lpf file with updated version
lpf.write_file()
#################################################################
## GHB, RCH, WEL scaling factors ##
scaling_factors = pd.read_csv('GHB_UZF_WEL_scaling.csv', delimiter = ',')

###############################################################################
## GHB Package ##
# join top and botm for easier array referencing for elevations
top_botm = np.zeros((m.dis.nlay+1,m.dis.nrow,m.dis.ncol))
top_botm[0,:,:] = m.dis.top.array
top_botm[1:,:,:] = m.dis.botm.array
# load ghb dataframes
ghbse_spd = pd.read_csv(all_model_ws+'/input_data/ghbse_spd.csv')
ghbnw_spd = pd.read_csv(all_model_ws+'/input_data/ghbnw_spd.csv')
ghbdelta_spd = pd.read_csv(all_model_ws+'/input_data/ghbdelta_spd.csv')
# edit GHB conductance by scaling with the new K vs original
ghb_hk_nw = scaling_factors.loc[0,'K_nw']*86400
ghb_hk_se = scaling_factors.loc[0,'K_se']*86400
K_delta = scaling_factors.loc[0,'K_delta']*86400

# calculate new conductance
# only need to recalculate conductance in ucode
def recalc_cond(i,j,k,hk):
    distance = 5000
    delr = m.dis.delr.array.mean()
    cond = hk*(top_botm[k,i,j]-top_botm[k+1,i,j])*delr/distance
    return(cond)
ghbse_spd.cond = recalc_cond(ghbse_spd.i.values,ghbse_spd.j.values,ghbse_spd.k.values, ghb_hk_se)
ghbnw_spd.cond = recalc_cond(ghbnw_spd.i.values,ghbnw_spd.j.values,ghbnw_spd.k.values, ghb_hk_nw)
ghbdelta_spd.cond = recalc_cond(ghbdelta_spd.i.values,ghbdelta_spd.j.values,ghbdelta_spd.k.values, K_delta)
# lay, row, col for delta ghb
zxy = ghbdelta_spd.values[:,:3].astype(int)
# drop any delta ghb cells where cell bottom is below sea level
ghbdn_spd =  ghbdelta_spd.values[botm[zxy[:,0],zxy[:,1],zxy[:,2]]<0]
# join dataframes of 3 ghb boundaries together
ghb_spd = np.vstack((ghbdn_spd, ghbse_spd.values, ghbnw_spd.values))

# allocate empty dictionary
ghb_dict = {}
ghb_dict[0] = ghb_spd

# create GHB for flopy
ghb = flopy.modflow.ModflowGhb(model = m,stress_period_data =  ghb_dict)
# overwrite the previous ghb file with updated version
ghb.write_file()

###############################################################################
## RIV Package ##

###############################################################################

riv_df = pd.read_csv(all_model_ws+'/input_data/riv_input.csv')
# apply scaling factor to river conductance
riv_s = scaling_factors.loc[0,'RIV']
riv_df.cond *= riv_s
# set dictionary and set riv
riv_dict = {0:riv_df.values}
riv = flopy.modflow.ModflowRiv(model = m, stress_period_data=riv_dict)
# overwrite previous riv file
riv.write_file()

# run the modflow model
success, buff = m.run_model()