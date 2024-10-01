# +

import sys
from os.path import basename, dirname, join, exists, expanduser
import glob

import pandas as pd
import numpy as np

# when dealing with larger data sets, it may be worthwhile using parallel subprocess
import shutil, os

# -

doc_dir = join(expanduser('~'), 'Documents')


def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# +
ucode_fxn_dir = doc_dir+'/GitHub/CosumnesRiverRecharge/ucode_utilities'
add_path(ucode_fxn_dir)
import ucode_input

# from importlib import reload
# reload(ucode_input)
# -

loadpth =  'C:/WRDAPP/GWFlowModel/Cosumnes/Regional'
model_ws = join(loadpth, 'historical_simple_geology_reconnection')
model_ws = join(loadpth, 'historical_geology_cal')


# # Load parameter data

# load geologic parameters
params = pd.read_csv(join(model_ws, 'ZonePropertiesInitial.csv'))
params = params.set_index('Zone')
# # convert from m/s to m/d
params['K_m_d'] = params.K_m_s * 86400   

# load other model parameters
bc_params = pd.read_csv(join(model_ws, 'BC_scaling.csv'))
bc_params = bc_params.set_index('ParamName')
if 'Scenario' in bc_params.columns:
    bc_params = bc_params[bc_params.Scenario != 'No reconnection'].drop(columns=['Scenario'])

# # Pgroup data

name = 'MF_ucode'
# name = 'MF_ucode_zone'
pgroup = pd.read_csv(model_ws+'/'+name+'.pgroup', delimiter=r'\s+', index_col='GroupName',
                     skipfooter=1, skiprows=2, engine='python')
adj_pgroup = pgroup[pgroup.Adjustable=='yes'].index.values


adj_pgroup

params_long = ucode_input.make_gel_p_long(params)
pdata = ucode_input.prep_gel_pdata(params_long)
pdata_bc = ucode_input.prep_gel_pdata(bc_params.reset_index())
pdata_std = pd.concat((pdata, pdata_bc))

pdata_zone = ucode_input.pdata_by_facies(pdata, params)
pdata_zone = pd.concat((pdata_zone, pdata_bc))

ucode_input.write_pdata(pdata_zone, model_ws, 'MF_ucode_zone.pdata')
ucode_input.write_pdata(pdata_std, model_ws, 'MF_ucode.pdata')

# # JTF

# +
# Write out jtf file
p_out = params.drop(columns=['K_m_d'])
p_out.K_m_s = '@'+('Kx_'+p_out.index.astype(str)).str.ljust(20)+'@'
# p_out.vani = '@'+('vani_'+p_out.index.astype(str)).str.ljust(20)+'@'
p_out.Sy = '@'+('Sy_'+p_out.index.astype(str)).str.ljust(20)+'@'
p_out.Ss = '@'+('Ss_'+p_out.index.astype(str)).str.ljust(20)+'@'
p_out.vani = '@'+('vani_'+p_out.index.astype(str)).str.ljust(20)+'@'

with open(model_ws+'/ZonePropertiesInitial.csv.jtf', 'w',newline='') as f:
    f.write('jtf @\n')
    p_out.to_csv(f,index=True, mode="a")
    
scaling_jtf = bc_params.reset_index().copy()
# scaling_jtf = scaling_factors_all.copy()
# Write out jtf file
scaling_jtf.StartValue = '@'+scaling_jtf.ParamName.str.ljust(20)+'@'

# with open(model_ws+'/GHB_UZF_WEL_scaling.csv.jtf', 'w',newline='') as f:
with open(model_ws+'/BC_scaling.csv.jtf', 'w',newline='') as f:
    f.write('jtf @\n')
    scaling_jtf.to_csv(f,index=False, mode="a")
# -

# # Observation data

# ## HOB jif

all_obs = pd.read_csv(join(model_ws, 'input_data','all_obs_grid_prepared.csv'), index_col=0)
# all_obs

# +
hobout = pd.read_csv(model_ws+'/MF.hob.out',delimiter = r'\s+')

# here std deviation represents the actual value one expects
# for a well the accuracy is 0.01 ft at best based on measuring tape scale
all_obs['Statistic'] = 0.01
all_obs['StatFlag'] = 'SD'
# locations with significant difference between RPE GSE and the DEM should have additional uncertainty included
all_obs['Statistic'] += np.round(np.abs(all_obs.dem_wlm_gse),4)
# correct extreme statistics with outlier
upper_limit = all_obs.Statistic.quantile(.75) +1.5*(all_obs.Statistic.quantile(0.75)-all_obs.Statistic.quantile(0.25))
all_obs.loc[all_obs.Statistic>upper_limit, 'Statistic'] = upper_limit
all_obs['Weight'] = 1/(all_obs.Statistic**2)

hobout_in = hobout.join(all_obs.set_index('obs_nam')[['Statistic','StatFlag']],on=['OBSERVATION NAME'])
# temporary fix for misjoin for single observation HOB nodes
hobout_in.loc[hobout_in.Statistic.isna(),'Statistic'] = 0.01 
# hobout_in['StatFlag'] = 'SD'

ucode_input.write_hob_jif_dat(model_ws, hobout_in, statflag=True)
# -

# ## Gage jif

flow_obs = pd.read_csv(join(model_ws, 'input_data','flow_obs_gage.csv'), index_col=0)


# +

# mcc_d = mcc_d[(mcc_d.index>=strt_date)&(mcc_d.index<=end_date)]
# ObsName ObsValue Statistic StatFlag GroupName
flow_obs['ObsName'] = (flow_obs.GroupName+'_'+pd.Series(np.arange(0,len(flow_obs)).astype(str)).str.zfill(5)).values
# make sure units are flow in m^3/day
flow_obs = flow_obs.rename(columns={'flow_cmd':'ObsValue'})

cols_out = ['ObsName','ObsValue','Statistic','StatFlag','GroupName']

header = 'BEGIN Observation_Data Table\n'+\
    'NROW= '+str(len(flow_obs))+' NCOL= 5 COLUMNLABELS\n'+\
    ' '.join(cols_out)

footer = 'End Observation_Data'
# get array of just strings
flow_arr = flow_obs[cols_out].values
# pull out observed value and name of obs
# np.savetxt(model_ws+'/flow_obs_table.dat', flow_arr,
#            fmt='%s', header = header, footer = footer, comments = '' )

# +
# need to pull in start date to set up gage file jif
# strt_date = 
# or need to add a column to flow_obs_gage for the time since start

# +
# gagenam = model_ws+'/MF.gage1.go'
# gage = pd.read_csv(gagenam,skiprows=1, delimiter = r'\s+', engine='python')
# # clean issue with column name read in
# cols = gage.columns[1:]
# gage = gage.dropna(axis=1)
# gage.columns = cols
# # set datetime for joining with flow obs data
# gage['dt'] = strt_date + (gage.Time-1).astype('timedelta64[D]')
# gage = gage.set_index('dt').resample('D').mean()
# gage_jif = gage[['Time','Flow']].join(mcc_d)
# # if I leave Nan values then ucode gets upset, Cab used the filler dum which I Ucode identifies to skip
# gage_jif.loc[gage_jif.ObsName.isna(),'ObsName'] = 'dum'
# for gage file JIFs need to specify which flows are used by specify the observation name
# for the correct row (time) and filling the rest with a dummy variable (Cab used dum)
# -

# def write_flw_jif(model_ws, gagout):
#     # skip 2 rows, use 3rd column values for 1345 values for std MF gage out file
header = 'jif @\nStandardFile 2 3 '+str(len(gage_jif))
# header = 'jif @\n'+'StandardFile  1  1  '+str(len(obsoutnames))
# obsoutnames.to_file(m.model_ws+'/MF.hob.out.jif', delimiter = '\s+', index = )
np.savetxt(model_ws+'/MF.gage1.go.jif', gage_jif.ObsName.values,
           fmt='%s', delimiter = r'\s+',header = header, comments = '')

# # Parallel ucode
# The 6 year run is only 1 hr 20 min by itself but in parallel on the laptop it slowed way down. taking likely 6 hrs. Better to just run SA with only 3 or 4 parallel runs.  

# find pdata that will be adjusted
adj_pdata = pdata_std[pdata_std.GroupName.isin(adj_pgroup)].copy()


# n_nodes = ucode_input.get_n_nodes(4)
n_nodes = adj_pdata.shape[0]


# 2400 seconds is about 40 minutes which is avg run time
# 2 hr 15 min is 8100
run_time = 3*3600 # 
# may need to extend upward to 3 hours (10800) for slow runs
ucode_input.write_parallel(model_ws, n_nodes, run_time) 

# +
# # copy mf files except cbc and hds
mf_files = pd.Series(glob.glob(model_ws+'/MF.*'))
mf_files = mf_files[~mf_files.str.contains('cbc|hds|chk').values].tolist()
jtfs = glob.glob(model_ws+'/*.jtf')
run = glob.glob(model_ws+'/*py*')

files = mf_files+jtfs+run
# files = glob.glob(model_ws+'/*.jtf')
# files = glob.glob(model_ws+'/*.sfr')
# files = glob.glob(model_ws+'/*.csv')

# files

# +


for n in np.arange(0, n_nodes).astype(str):
    folder = '/r'+ n.zfill(3)+'/'
    os.makedirs(model_ws+folder,exist_ok=True)
    for f in files:
        shutil.copy(f, model_ws+folder)


# +
# replace oc file with simplified version that only prints the budget monthly
f = glob.glob(model_ws+'/MF_parallel.oc')[0]

for n in np.arange(0, n_nodes).astype(str):
    folder = '/r'+ n.zfill(3)+'/'
    os.makedirs(model_ws+folder,exist_ok=True)
    shutil.copy(f, model_ws+folder+'/MF.oc')
# -

# write out just the updated python write file
write_file = glob.glob(model_ws+'/*.py')
for n in np.arange(0,n_nodes).astype(str):
    folder = '/r'+ n.zfill(3)
    shutil.copy(write_file[0], model_ws+folder)


