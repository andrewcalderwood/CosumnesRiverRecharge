# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
# standard python utilities
import os
from os.path import join, exists, dirname, basename, expanduser
import glob
import sys
import time
from importlib import reload
import h5py

import pandas as pd
import numpy as np

# standard geospatial python utilities
import shapely
import geopandas as gpd

import matplotlib.pyplot as plt


# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')

# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

# # flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# # other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)
from mf_utility import get_layer_from_elev, param_load

import Basic_soil_budget_monthly as swb
# from importlib import reload
# reload(swb)
# from Basic_soil_budget_monthly import calc_yield, calc_profit, prep_soil
# , calc_S, calc_pc


# %%
# resampled ground surface elevation
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_mean.tsv')


# %%
year = int(2015)
# crop='Grape'
crop='Corn'
# crop='Alfalfa'
crop='Pasture' # will require extra work due to AUM vs hay
# crop = 'Misc Grain and Hay'

# %%
var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)

# %%

#  get the dates for each yield cycle
yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
# get the total extent of the irrigation season (calculation period)
strt_date = yield_start.min()
end_date = yield_end.max()
dates = pd.date_range(strt_date, end_date, freq='D')
nper = (end_date-strt_date).days +1


# %%
# the model will run the irrigation optimizer on specified dates (multiple crops can be done at once or in sequence)
# the modflow model will be run for the periods between the specified irrigation optimizer dates
loadpth = 'F://WRDAPP/GWFlowModel/Cosumnes/Regional/'

model_ws = loadpth + 'crop_soilbudget/'+crop+'_'+str(strt_date.date())


# %%
# fn = join(data_dir,'static_model_inputs.xlsx')

Kc, Kc_dates = swb.load_Kc(year)

Kc_c = Kc[Kc.Crop==crop]
Kc_dates_c = Kc_dates[Kc_dates.Crop==crop]

Kc_df = swb.get_Kc_dates(Kc_dates_c, Kc_c)



# %%
rain, ETo_df = swb.load_hyd(year, dates)

# get the crop ET
ETc = ETo_df.values*Kc_df.Kc.values

# %%
## Load crop type predictions

# double check with Yusuke whether crop2016 or crop_cat2016 is the original or predicted
crop_wide = pd.read_excel(join(proj_dir, 'crop_prediction', 'wide_dataset.xls'))
# if a dataset has a diversion then enable decision between SW and GW
# very few parcels are considered to have a POD (0.5%)
frac_no_pod = (crop_wide.pod=='No Point of Diversion on Parcel').sum()/len(crop_wide)
print('%.2f parcels have no POD' %(frac_no_pod*100))
# crop_wide.waterbody_type.unique()
# filter to relevant columns for SWB
keep_cols = crop_wide.columns.str.extract(r'(crop\d{4})').dropna().iloc[:,0].tolist()
crop_wide = crop_wide.set_index(['parcel_id','pod'])[keep_cols]
# crop_wide = crop_wide.set_index(['parcel_id','pod'])[['crop2016']]

# %%
# long format then filter for year
crop_long = crop_wide.melt(ignore_index=False,value_name='name', var_name='crop_yr')
crop_long['year'] = crop_long.crop_yr.str.extract(r'(\d+)').astype(int)


# %%
d_ini = 0; # Initial root zone depletion (in)
y_max = var_crops[['y_max']].values # Maximum expected crop yield (tons/acre)

phi = var_gen['phi'] # Energy requirement to raise a unit of water by a unit of vertical distance (kWh/acre-in/ft)
p_c = var_crops['p_c'] # Crop price ($/ton)
p_sw = var_gen['p_sw'] # Surface water charges and fees ($/acre-in)
p_e = var_gen['p_e'] # Cost of energy for groundwater pumping ($/kWh)
p_o = var_crops['p_o'] # Variable operating costs per acre, excluding irrigation costs ($/acre)



# %%
# create dictionary of crop variables and join multiples into list
var_crops_dict = {}
for v in var_crops.index.unique():
    var_crops_dict[v] = np.asarray(var_crops[v].tolist())

# create dictionary of general parameters
gen_dict = {**var_gen.to_dict(), **var_crops_dict}

# %%
gen_dict['y_max'] = var_crops[['y_max']].values # y_max may be multiple so keep array structure
gen_dict['nper'] = (end_date-strt_date).days +1

# %%
# days to index for calculating yield impacts, add 1 to include end day
yield_ind = np.append([0], (yield_end-strt_date).dt.days.values +1)
gen_dict['yield_ind'] = np.append([0], (yield_end-strt_date).dt.days.values +1)

# %%
# create time series of daily yield response factors
var_yield['dt'] = swb.ymd2dt(year, var_yield.month, var_yield.day, var_yield.year_adj)
# adjust for dates in the next year
K_Y = var_yield.set_index('dt').reindex(dates).ffill()
K_Y = K_Y['K_Y'].values


# %%

gap_irr = var_crops['gap_irr'] # Number of days between irrigations
n_irr = np.floor(len(dates)/gap_irr).astype(int) + 1 # Calculate number of irrigations
irr_days = np.arange(0, (n_irr*gap_irr-1), gap_irr).astype(int) # Calculate days on which irrigation takes place
irr_days


# %%
gen_dict['K_Y'] = var_yield.set_index('dt')['K_Y'].reindex(dates).ffill().values
gen_dict['n_irr'] = n_irr
gen_dict['irr_days'] = np.arange(0, (n_irr*gap_irr-1), gap_irr).astype(int) # Calculate days on which irrigation takes place

# %% [markdown]
# ## Test SWB on one field

# %%
# decide later best way to identify but in future the predicted data will be available for each year
crop_in = crop_long[crop_long.year==2016].reset_index()

# %%
avg_irr_eff = pd.read_csv(join(proj_dir, 'model_inputs', 'avg_irr_eff_by_crop.csv'),index_col=0)
irr_eff_mult = 100/avg_irr_eff.loc[crop_dict[crop]].Avg_eff
gen_dict['irr_eff_mult'] = 100/avg_irr_eff.loc[crop_dict[crop]].Avg_eff

soil_crop = swb.load_soil(pred_dict[crop], crop_in)
nfield_crop = len(soil_crop)


# %%
soil_path = join(uzf_dir,'clean_soil_data')
# connection of ag fields to all grid cells
grid_soil = pd.read_csv(join(soil_path, 'parcel_field_to_cell.csv'),index_col=0)
# subset to fields for crop
grid_crop = grid_soil.merge(soil_crop[['UniqueID']])


# %%
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Regional/'
model_ws = loadpth+'historical_simple_geology_reconnection'
model_ws = loadpth+'strhc1_scale'
load_only = ['DIS','BAS6']
m = flopy.modflow.Modflow.load('MF.nam', model_ws=model_ws, 
                                exe_name='mf-owhm.exe', version='mfnwt', load_only=load_only)
botm = np.copy(m.dis.botm.array)

# %%
# also need shapefile of pumping well locations for each parcel
parcel_wells = gpd.read_file(join(gwfm_dir, 'WEL_data', 'parcels_to_wells', 'parcels_to_wells.shp'))
frow = parcel_wells.row-1
fcol = parcel_wells.column-1
# # parcel_wells layers (make 1-based
parcel_wells['layer'] = get_layer_from_elev(dem_data[frow,fcol] - parcel_wells.depth_m*0.9, botm[:, frow,fcol], m.dis.nlay) + 1
# subset to crops
parcel_wells.UniqueID = parcel_wells.UniqueID.astype(int)
crop_wells = soil_crop[['UniqueID']].merge(parcel_wells)
crop_dem = dem_data[crop_wells.row-1, crop_wells.column-1]

# %%
# sample groundwater elevations from model output
# it might be more time efficient to use the hob.out package for field centers (well locations)
# rather than sampling the hds which might take more time

# sample heads
# model_ws_last = loadpth + 'crop_modflow/'+str(all_run_dates.loc[m_per-1].date.date())
# hdobj = flopy.utils.HeadFile(model_ws_last + '/MF.hds')
hdobj = flopy.utils.HeadFile(model_ws + '/MF.hds')
# sp_last = hdobj.get_kstpkper()[-1]


# %%
m_strt = pd.to_datetime(str(year-1)+'-10-1')
m_strt
# m_end = all_run_dates.iloc[m_per+1].date

# %%
# it's going to be very slow to sample the time series for all wells for all stress periods and time steps
# it does seem more efficient to sample on a monthly scale like a farmer would (irrigation dates) and interpolate
well_kij = list(zip(parcel_wells.layer-1, parcel_wells.row-1, parcel_wells.column-1))
# hdobj.get_ts(well_kij)
head = hdobj.get_data(kstpkper = (0,0))
head[crop_wells.layer-1, crop_wells.row-1, crop_wells.column-1].shape

# %%
# for spd in np.arange(1,
crop_dtw = np.zeros((nper,nfield_crop))
spd_strt = (strt_date-m_strt).days
dtw_days = np.hstack((irr_days,[nper-1])) # irr_days alternate
for n, kper in enumerate(dtw_days):
# for kper in np.arange(0,nper): # daily sampling
    head = hdobj.get_data(kstpkper = (0,kper+spd_strt))
    crop_dtw[kper,:] = crop_dem - head[crop_wells.layer-1, crop_wells.row-1, crop_wells.column-1]
    # if n >=1: # not necessary except for plotting
    #     # fill in between as farmer might expect linear change
    #     crop_dtw[irr_days[n-1]:kper+1,:] = np.linspace(crop_dtw[irr_days[n-1]], crop_dtw[kper], kper-irr_days[n-1]+1)
# convert from meters to feet
crop_dtw /= 0.3048

# %%
plt.plot(crop_dtw[:,0])


# %%
# temporary input for dtw variables to calculate pumping cost
dtw_time = np.linspace(90,95, nper)
dtw_all = np.zeros((nper,nfield_crop))
for n in np.arange(nper):
    dtw_all[n,:] = dtw_time[n]

# %%
# plt.plot(dtw_time)

# %% [markdown]
# ## Iterate over each unique soil condition

# %%
import swb_functions
reload(swb_functions)
from swb_functions import prep_soil_dict, calc_S, calc_pc, calc_yield, calc_profit,  run_swb, mak_irr_con


# %%
# so if you have a dictionary d and want to access (read) its values with the syntax x.foo instead of the clumsier d['foo'], just do
# convert a dictionary to an object with object style referencing
class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)

# # example
# d = {'a':1}
# x = Bunch(d)
# x.a


# %%
# convert dictionary of variables to class for easier referencing, constant over different soil
gen = Bunch(gen_dict)

from scipy.optimize import Bounds
bounds = Bounds(lb = 0)
from scipy.optimize import minimize

# %%
ns=0
soil_ag = soil_crop.iloc[[ns]] #keep as dataframe for consistency 
nfield = soil_ag.shape[0]

dtw_arr = dtw_all[:,ns]

etc_arr = np.zeros((nper,nfield))
for n in np.arange(nper):
    etc_arr[n,:] = ETc[n]

# prep_soil(soil_ag, etc_arr, var_crops)
# conver to class for referencing easier
# benefit of the class is no indexes like pd.DataFrame
soil_dict = prep_soil_dict(soil_ag, etc_arr, var_crops)
soil = Bunch(soil_dict)


# %%
irr_lvl = np.zeros(2*n_irr); # Initial irrigation values for optimization
irr_lvl[:] = (2/12)*0.3048 # irrigate with 2 inches (convert to meters)

# %%
# p_gw = p_e*phi*dtw_arr[irr_days]
# # where the gw costs more then set gw as 0 and vice versa
# irr_lvl[:n_irr][p_gw > p_sw] = 0
# irr_lvl[n_irr:][p_gw < p_sw] = 0


# %%
# irr_sw = np.zeros((nper,nfield_crop))
# irr_gw = np.zeros((nper,nfield_crop))
# for i in np.arange(0,n_irr):
#     irr_sw[irr_days[i]] = irr_lvl[i]
#     irr_gw[irr_days[i]] = irr_lvl[i+n_irr]


# %%

# %%
etc_arr = np.zeros((nper))
for n in np.arange(nper):
    etc_arr[n] = ETc[n]

for ns in [0]:#np.arange(0,len(soil_ag_unique)):
    soil_ag = soil_crop.iloc[[ns]] #keep as dataframe for consistency 
    nfield = soil_ag.shape[0]
    
    # prep_soil(soil_ag, etc_arr, var_crops)
    soil_dict = prep_soil_dict(soil_ag, etc_arr, var_crops)
    soil = Bunch(soil_dict)
    
    # dtw_arr = dtw_all[:,ns]
    dtw_arr = crop_dtw[:,ns]
    irr_lvl = np.zeros(2*n_irr); # Initial irrigation values for optimization
    irr_lvl[:] = (2/12)*0.3048 # irrigate with 2 inches (convert to meters)
    # if no POD then no SW irrig
    if soil_ag.pod.iloc[0]=='No Point of Diversion on Parcel':
        irr_lvl[:n_irr] = 0
        irr_lvl[n_irr:] *= 2 # put double the irrigation to the GW

    irr_lvl[:] *= irr_eff_mult # scale by irrigation efficiency of the crop after optimizing

    pi = run_swb(irr_lvl, soil, gen, rain, ETc, dtw_arr)
    print('Profit %.2f ($/acre)' %(-pi))

# %%
# pi, pc, K_S = run_swb(irr_lvl, soil, gen, rain, ETc, dtw_arr, arrays=True)

# Y_A = calc_yield(ETc, K_S, gen)
# Y_A


# %%
# irr_lvl = np.copy(out.x)
# irr_sw = np.zeros((nper,nfield))
# irr_gw = np.zeros((nper,nfield))
# for i in np.arange(0,n_irr):
#     irr_sw[irr_days[i]] = irr_lvl[i]
#     irr_gw[irr_days[i]] = irr_lvl[i+n_irr]
# calc_profit(Y_A, dtw_arr, irr_gw, irr_sw, gen)  

# %%
from scipy.optimize import Bounds, LinearConstraint
from scipy.optimize import minimize

# %% [markdown]
# - Total surface water and groundwater available during the season (in)
# - to represent unconstrained conditions set boundaries at 100, 100
#

# %%


bounds = Bounds(lb = 0)

linear_constraint = mak_irr_con(soil_ag, n_irr, sw_con=100, gw_con=100)

# %% [markdown]
# Test minimizer tolerance for estimating optimal irrigation.

# %%
for tol in [0.00001, 0.0001, 0.001, 0.01, 0.1]:
    out = minimize(run_swb, irr_lvl, args = (soil, gen, rain, ETc, dtw_arr), method='trust-constr',
            constraints = [linear_constraint],
            bounds=bounds,
            tol = tol,
    #          options={'verbose':1}
            )
    print('For tol:'+str(tol)+', Irr was %.2f (in)' %(out.x.sum()*12/0.3048), 'for %.3f $' %(-out.fun),'and total time was %.2f' %out.execution_time )
    


# %% [markdown]
# From 1E-5 to 0.01 tolerance the irrigation and profit are within 0.1 which is pretty good. Time goes from 110 to 35 sec which is alright, but not great if there are 4,100 fields (all crops) because then it is 38 min

# %%
tol = 0.01
out = minimize(run_swb, irr_lvl, args = (soil, gen, rain, ETc, dtw_arr),method='trust-constr',
        constraints = [linear_constraint],
        bounds=bounds,
        tol = tol,
#          options={'verbose':1}
        )    
print('Optimization time %.2f sec' %out.execution_time)

# %%
tol = 0.01
t0 = time.time()
out_new = minimize(run_swb, irr_lvl, args=(soil, gen, rain, ETc, dtw_arr), 
         #jac=func_deriv, 
         constraints=[linear_constraint], 
         bounds = bounds,
         tol = 1E-4,
         method='SLSQP', options={'disp': True})

t1 = time.time()
print('Optimization time %.2f sec' %(t1-t0))


# %% [markdown]
# SLSQP runs slightly faster (77 to 58 sec), but requires a tighter tolerance because with 1E-2 it didn't finish solving and at 1E-4 it did (1E-3 was slightly off too).
# - not a clear way to estiamte jacobian (gradient) since it is a pretty nonlinear function

# %%
# pi, pc, K_S = run_swb(out.x, soil, gen, rain, ETc, dtw_arr, arrays=True)
# Y_A = calc_yield(ETc, K_S, gen)

# %% [markdown]
# Test irrigation efficiency impact. The solver can be run for all fields then irrigation efficiency can be applied before re-running all SWBs.

# %%
def check_true_swb(irr_out):
    irr_true = irr_out * irr_eff_mult # scale by irrigation efficiency of the crop after optimizing
    p_true = run_swb(irr_true, soil, gen, rain, ETc, dtw_arr)
    p = run_swb(irr_out, soil, gen, rain, ETc, dtw_arr)
    print('True profit %.2f' %(-p_true), 'profit no irr. inefficiency %.2f' %(-p))

check_true_swb(out.x)
check_true_swb(out_new.x)


# %%
import matplotlib.pyplot as plt


# %%
def summarize_out(irr_out):
    sw_irr_in= (irr_out[:n_irr]/0.3048).sum()*12
    gw_irr_in = (irr_out[n_irr:]/0.3048).sum()*12
    print('Irr (in) GW %.2f'%gw_irr_in, 'SW %.2f'%sw_irr_in)
    # print('Cost ($/acre) GW %.2f' %(c_gwtot.sum()),'SW %.2f'%(c_swtot.sum()))
    
    plt.bar(np.arange(0,n_irr)-.25/2,(irr_out[:n_irr]/0.3048)*12, label='SW', width=0.25)
    plt.bar(np.arange(0,n_irr)+.25/2,(irr_out[n_irr:]/0.3048)*12, label='GW', width=0.25)
    plt.legend()
    plt.ylabel('Irrigation (in)')
    plt.xlabel('Irrigation Event')
    plt.show()


# %%
summarize_out(out.x)
summarize_out(out_new.x)


# %% [markdown]
# # Solving over fields
# - tol = 0.01 is sufficient to have most solver runs work and tightening tolerance by orders of 10 is enough to avoid negatives, need to insert break for extreme tol which would lead to never ending. (Run time for corn was 30 min, another time was 20 min). After switching to new functions which require a dictionary then object it took 60 min and 48 min another time.
# - If we go to parallel solving then the irrigation optimizer can be grouped by crop then major soil group so that the starting irrigation rates are most similar.

# %%

# %%
t0 = time.time()
# save irrigation, function, time
irr_all =np.zeros((nfield_crop,len(irr_lvl)))
p_all = np.ones(nfield_crop) # make initial profit negative
t_all = np.zeros(nfield_crop)

# tol = 0.01
irr_lvl = np.zeros(2*n_irr); # Initial irrigation values for optimization
irr_lvl[:] = (2/12)*0.3048 # irrigate with 2 inches (convert to meters)
irr_lvl_base = np.copy(irr_lvl)

etc_arr = np.zeros((nper))
for n in np.arange(nper):
    etc_arr[n] = ETc[n]

# for ns in np.arange(0,nfield_crop):
for ns in np.arange(0,10):
    soil_ag = soil_crop.iloc[[ns]] #keep as dataframe for consistency 
    nfield = soil_ag.shape[0]

    # dtw_arr = dtw_all[:,ns]
    dtw_arr = crop_dtw[:,ns]

    # prep_soil(soil_ag, etc_arr, var_crops)
    soil_dict = prep_soil_dict(soil_ag, etc_arr, var_crops)
    soil = Bunch(soil_dict)
    
    if ns > 1:
        irr_lvl[:] = irr_all[ns-1]
    # if no POD then no SW irrig
    if soil_ag.pod.iloc[0]=='No Point of Diversion on Parcel':
        irr_lvl[:n_irr] = 0
        irr_lvl[n_irr:] *= 2 # put double the irrigation to the GW

    linear_constraint = mak_irr_con(soil_ag, n_irr) # will change for POD
    tol = 0.01
    # continue optimizing until profit is positive
    while p_all[ns] >0 :
        # the minimization with 'trust-constr' and no constraints doesn't solve and has increasing WB error
        out = minimize(run_swb, irr_lvl, args = (soil, gen, rain, ETc, dtw_arr), method='trust-constr',
                constraints = [linear_constraint],
                bounds=bounds,
        #          options={'verbose':1}
                       tol = tol
                )
        # decrease tolerance and reset starting irrigation to help solving
        if out.fun >0:
            tol /=10
            irr_lvl[:] = np.copy(irr_lvl_base)
        if tol < 1E-5:
            break # if tolerance gets too small then skip
        irr_all[ns] = out.x
        p_all[ns] = out.fun
        t_all[ns] = out.execution_time
    # print final results
    print('Soil ', str(ns),'%.2f' %(-out.fun),'$ ,in %.2f' %(out.execution_time/60),'min')
t1 = time.time()
print('Total time was %.2f min' %((t1-t0)/60))

# %% [markdown]
# Once the solution from the previous field is applied as a starting point it seems that the solutions run about 10 times faster so it should only take 3.8 min instead of 38 min. There was one solution that ended negative, early on and might have been because it went the wrong way with the previous solution so may need to have an if statement to run again with a tighter tolerance if result is negative.
#
# Tightening the tol to 0.001 still had quite a few negative values and slowed up the execution time from 0.04 to 0.07. It might be more important to have that if statement then with a reset of initial values. It could be worth comparing which field combinations have more consistent values

# %%
t_all.sum()/60/60

# %%
lines =plt.plot(np.transpose(irr_all))


# %%
fig,ax = plt.subplots(5,1,sharex=True)
ax[0].plot(rain)
# ax[0].plot(irr_sw[:,0])
ax[0].bar(x=np.arange(0,len(irr_sw))-0.25, height=irr_sw[:,0], width=1,label='SW')
ax[0].bar(x=np.arange(0,len(irr_gw))+0.25, height=irr_gw[:,0], width=.5, label='GW')
ax[1].plot(rp.mean(axis=(1)))
ax[2].plot(ETa.mean(axis=(1)))
ax[3].plot(pc.mean(axis=(1)))
ax[4].plot(wc.mean(axis=(1)))

for n, l in enumerate(['Rain/\nIrrigation','Runoff','ET','Perc','WC']):
    ax[n].set_ylabel(l)
ax[0].legend()

# %% [markdown]
# ## Save Output

# %%

import h5py

def crop_arr_to_h5(arr, crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "w") as f:
        grp = f.require_group('array') # makes sure group exists
        grp.attrs['units'] = 'meters/day'
        grp.attrs['description'] = 'Each layer of the array is a day in the water year'
        dset = grp.require_dataset(crop, arr.shape, dtype='f', compression="gzip", compression_opts=4)
        dset[:] = arr
    


# %%
# years and array index 
years = pd.date_range(strt_date,end_date,freq='AS-Oct')
yr_ind = (years-strt_date).days


# %%
crop

# %%

# # need separte hdf5 for each year because total is 300MB
# for n in np.arange(0,len(yr_ind)-1):
#     arr = pc[yr_ind[n]:yr_ind[n+1]]
#     fn = join(model_ws,"percolation_WY"+str(years[n].year+1)+".hdf5")
#     arr_to_h5(arr, fn)

# need separte hdf5 for each year because total is 300MB, group by crop in array
for n in np.arange(0,len(yr_ind)-1):
#     arr = pc[yr_ind[n]:yr_ind[n+1]]
    arr = pc_arr[yr_ind[n]:yr_ind[n+1]]
    fn = join(model_ws, "field_percolation_WY"+str(years[n].year+1)+".hdf5")
    arr_to_h5(arr, fn)

    # applied water 
    arr = AW[yr_ind[n]:yr_ind[n+1]]
    fn = join(model_ws, "field_applied_water_WY"+str(years[n].year+1)+".hdf5")
    arr_to_h5(arr, fn)
