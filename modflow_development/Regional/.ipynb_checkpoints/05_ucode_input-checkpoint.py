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

all_obs = pd.read_csv(join(model_ws, 'input_data','all_obs_grid_prepared.csv'), index_col=0)
# all_obs

# +
hobout = pd.read_csv(model_ws+'/MF.hob.out',delimiter = r'\s+')

# here std deviation represents the actual value one expects
# for a well the accuracy is 0.01 ft at best based on measuring tape scale
# all_obs['Statistic'] = 0.01
# all_obs['StatFlag'] = 'SD'
# locations with significant difference between RPE GSE and the DEM should have additional uncertainty included
# all_obs['Statistic'] += np.round(np.abs(all_obs.dem_wlm_gse),4)
# hob_gpd['Weight'] = 1/(hob_gpd.Statistic**2)

hobout_in = hobout.join(all_obs.set_index('obs_nam')[['Statistic','StatFlag']],on=['OBSERVATION NAME'])
# temporary fix for misjoin for single observation HOB nodes
hobout_in.loc[hobout_in.Statistic.isna(),'Statistic'] = 0.01 
hobout_in['StatFlag'] = 'SD'

ucode_input.write_hob_jif_dat(model_ws, hobout_in, statflag=True)
# -

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
# -


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
