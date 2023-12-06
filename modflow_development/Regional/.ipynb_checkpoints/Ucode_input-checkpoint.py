# %%
from os.path import join, expanduser
import sys
import glob

import pandas as pd
import numpy as np


# %%

# %%
loadpth = 'F:/WRDAPP/GWFlowModel/Cosumnes/Regional/'
# model_ws = loadpth+'historical_simple_geology_reconnection'
model_ws = loadpth+'steadystate'


# %% [markdown]
# # UCODE Input

# %%
name = 'MF_ucode'
# name = 'MF_ucode_zone'
pgroup = pd.read_csv(model_ws+'/'+name+'.pgroup', delimiter=r'\s+', index_col='GroupName',
                     skipfooter=1, skiprows=2, engine='python')
pgroup.loc['GHB','Adjustable']
pgroup.index
# if pgroup.loc['UZF','Adjustable'] =='yes':
#     print('True')

# %% [markdown]
# ## Format parameter data (pdata) file

# %%
def get_magnitude(x):
    return(10.0**(np.log10(x).astype(int)))

def prep_gel_pdata(params):
    # melt parameter data and rename columns to fit UCODE format for .pdata
    pdata = params.copy().rename(columns={'K_m_s':'Kx'})[['Kx','vani','Ss','Sy']]
    pdata = pdata.melt(ignore_index=False)
    pdata['ParamName'] = pdata.variable + '_' + pdata.index.astype(str)
    pdata = pdata.rename(columns={'variable':'GroupName','value':'StartValue'}).reset_index(drop=True)
    # join scaling factors to hydraulic parameters
    pdata = pd.concat((pdata, bc_params.reset_index()))
    
    # default values for pdata input
    pdata['LowerValue'] = 1E-38
    pdata['UpperValue'] = 1E38
    
    # local adjustment based on typical parameter scaling (start value scaled by a range)
    # need to find a better rounding function
    grps = pdata.GroupName.isin(['Kx','Ss','vani','GHB'])
    pdata.loc[grps,'LowerValue'] = get_magnitude(pdata.loc[grps,'StartValue']) *1E-3
    pdata.loc[grps,'UpperValue'] = get_magnitude(pdata.loc[grps,'StartValue']) *1E3
    grps = pdata.GroupName.isin(['Sy'])
    pdata.loc[grps,'LowerValue'] = get_magnitude(pdata.loc[grps,'StartValue']) *1E-2
    pdata.loc[grps,'UpperValue'] = 1
    grps = pdata.ParamName.str.contains('rch_')
    pdata.loc[grps,'LowerValue'] = get_magnitude(pdata.loc[grps,'StartValue']) *1E-3
    pdata.loc[grps,'UpperValue'] = 2
    
    # assume constraints align with expected range
    pdata['Constrain'] = 'No'
    pdata['LowerConstraint'] = pdata.LowerValue
    pdata['UpperConstraint'] = pdata.UpperValue
    return(pdata)
    
pdata = prep_gel_pdata(params)


# %% [markdown]
# Version for grouping by aquifer unit.

# %%
def pdata_by_facies(pdata):
    pdata_zone = pdata[pdata.GroupName.isin(['Kx','vani','Ss','Sy'])].copy()
    # alternate pdata where group is the lithology
    pdata_zone['Zone'] = pdata_zone.ParamName.str.extract(r'(\d)')
    pdata_zone.Zone = pd.to_numeric(pdata_zone.Zone)
    pdata_zone = pdata_zone.join(params[['Lithology']], on='Zone')
    pdata_zone['GroupName'] = pdata_zone.Lithology
    pdata_zone.loc[pdata_zone.Lithology.isin(['Gravel','Sand','Sandy Mud','Mud']),'GroupName'] = 'tprogs'
    pdata_zone = pdata_zone.drop(columns=['Zone','Lithology'])
    pdata_zone = pdata_zone.dropna(subset='GroupName')
    return(pdata_zone)
pdata_zone = pdata_by_facies(pdata)


# %%
def write_pdata(pdata, name):
    """ Write out pdata file """
    with open(join(model_ws, name), 'w',newline='') as f:

        # 27 before rch_1 to rch_12, 6 more for vani
        f.write('BEGIN Parameter_Data Table\n')
        f.write('NROW='+str(pdata.shape[0])+' NCOL='+str(pdata.shape[1])+' COLUMNLABELS\n')
        f.write(pdata.columns.str.ljust(12).str.cat(sep = ' '))
        f.write('\n')
        for n in np.arange(0, len(pdata)):
    #         f.write(pdata_zone.iloc[n].str.cat())
            f.write(pdata.iloc[n].astype(str).str.ljust(12).str.cat(sep=' '))
            f.write('\n')
        f.write('END Parameter_Data Table')
    print('Wrote pdata file')
write_pdata(pdata_zone, 'MF_ucode_zone.pdata')
write_pdata(pdata, 'MF_ucode.pdata')

# %%
# pdata.ParamName.str.contains('rch_')
# pdata.ParamName.isin([r'rch_.'])


# %% [markdown]
# ## JTF files

# %%
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
# Write out jtf file
scaling_jtf.StartValue = '@'+scaling_jtf.ParamName.str.ljust(20)+'@'

# with open(model_ws+'/GHB_UZF_WEL_scaling.csv.jtf', 'w',newline='') as f:
with open(model_ws+'/BC_scaling.csv.jtf', 'w',newline='') as f:
    f.write('jtf @\n')
    scaling_jtf.to_csv(f,index=False, mode="a")

# %% [markdown]
# ## Observation data

# %%

# %%
doc_dir = join(expanduser('~'), 'Documents')
ucode_fxn_dir = doc_dir+'/GitHub/CosumnesRiverRecharge/ucode_utilities'
if ucode_fxn_dir not in sys.path:
    sys.path.append(ucode_fxn_dir)
# sys.path
import ucode_input

# reload(ucode_input)

# %%
# save cleaned data to input data for post-processing
all_obs = pd.read_csv(model_ws+'/input_data/all_obs_grid_prepared.csv')

# %%
hobout = pd.read_csv(model_ws+'/MF.hob.out',delimiter = r'\s+')

# here std deviation represents the actual value one expects
# for a well the accuracy is 0.01 ft at best based on measuring tape scale
all_obs['Statistic'] = 0.01
all_obs['StatFlag'] = 'SD'
# locations with significant difference between RPE GSE and the DEM should have additional uncertainty included
all_obs['Statistic'] += np.round(np.abs(all_obs.dem_wlm_gse),4)

hobout_in = hobout.join(all_obs.set_index('obs_nam')[['Statistic','StatFlag']],on=['OBSERVATION NAME'])
# temporary fix for misjoin for single observation HOB nodes
hobout_in.loc[hobout_in.Statistic.isna(),'Statistic'] = 0.01 
hobout_in['StatFlag'] = 'SD'

ucode_input.write_hob_jif_dat(model_ws, hobout_in, statflag=True)

# %%
sfr_dir = gwfm_dir+'/SFR_data/'

# %%
# data for obs table
mcc_d = pd.read_csv(sfr_dir+'MCC_flow_obs_all.csv', index_col='DATE TIME', parse_dates=True)

mcc_d = mcc_d[(mcc_d.index>=strt_date)&(mcc_d.index<=end_date)]
# ObsName ObsValue Statistic StatFlag GroupName
mcc_d['ObsName'] = ('mcc_'+pd.Series(np.arange(0,len(mcc_d)).astype(str)).str.zfill(5)).values
# make sure units are flow in m^3/day
mcc_d = mcc_d.rename(columns={'flow_cmd':'ObsValue'})

cols_out = ['ObsName','ObsValue','Statistic','StatFlag','GroupName']

header = 'BEGIN Observation_Data Table\n'+\
    'NROW= '+str(len(mcc_d))+' NCOL= 5 COLUMNLABELS\n'+\
    ' '.join(cols_out)

footer = 'End Observation_Data'
# get array of just strings
flow_arr = mcc_d[cols_out].values
# pull out observed value and name of obs
np.savetxt(model_ws+'/flow_obs_table.dat', flow_arr,
           fmt='%s', header = header, footer = footer, comments = '' )

# for gage file JIFs need to specify which flows are used by specify the observation name
# for the correct row (time) and filling the rest with a dummy variable (Cab used dum)

# %%
## not set up quite??

# gagenam = model_ws+'/MF_mcc.go'
# gage = pd.read_csv(gagenam,skiprows=1, delimiter = r'\s+', engine='python')
# # clean issue with column name read in
# cols = gage.columns[1:]
# gage = gage.dropna(axis=1)
# gage.columns = cols
# # set datetime for joining with flow obs data
# gage['dt'] = strt_date + (gage.Time-1).astype('timedelta64[D]')
# gage = gage.set_index('dt').resample('D').mean()
# gage_jif = gage[['Time','Flow']].join(mcc_d)
# # if I leave Nan values then ucode gets upset, Cab used the filler dum which I think Ucode identifies
# gage_jif.loc[gage_jif.ObsName.isna(),'ObsName'] = 'dum'

# %%
# def write_flw_jif(model_ws, gagout):
#     # skip 2 rows, use 3rd column values for 1345 values for std MF gage out file
header = 'jif @\nStandardFile 2 3 '+str(len(gage_jif))
# header = 'jif @\n'+'StandardFile  1  1  '+str(len(obsoutnames))
# obsoutnames.to_file(m.model_ws+'/MF.hob.out.jif', delimiter = '\s+', index = )
np.savetxt(model_ws+'/MF_mcc.go.jif', gage_jif.ObsName.values,
           fmt='%s', delimiter = r'\s+',header = header, comments = '')

# %% [markdown]
# ## Parallel ucode

# %%
n_nodes = ucode_input.get_n_nodes(13)

# %%
# 2400 seconds is about 40 minutes which is avg run time
# 2 hr 15 min is 8100
# may need to extend upward to 3 hours (10800) for slow runs
# steady state is 4 min so max about 20 min
ucode_input.write_parallel(model_ws, n_nodes,1200) 

# %%
# # copy mf files except cbc and hds
mf_files = pd.Series(glob.glob(model_ws+'/MF.*'))
mf_files = mf_files[~mf_files.str.contains('cbc|hds|.hob.out|.sfr.out').values].tolist()
jtfs = glob.glob(model_ws+'/*.jtf')
run = glob.glob(model_ws+'/*py*')

files = mf_files+jtfs+run
# files = glob.glob(m.model_ws+'/*.jtf')

# %%
# when dealing with larger data sets, it may be worthwhile using parallel subprocess
import shutil, os


for n in np.arange(0, n_nodes).astype(str):
    folder = '/r'+ n.zfill(3)+'/'
    os.makedirs(model_ws+folder,exist_ok=True)
    for f in files:
        shutil.copy(f, model_ws+folder)

# %%
# replace oc file with simplified version that only prints the budget monthly
f = glob.glob(model_ws+'/MF_parallel.oc')[0]

for n in np.arange(0, n_nodes).astype(str):
    folder = '/r'+ n.zfill(3)+'/'
    os.makedirs(model_ws+folder,exist_ok=True)
    shutil.copy(f, model_ws+folder+'/MF.oc')

# %%
import shutil, os

# write out just the updated python write file
write_file = glob.glob(model_ws+'/*.py')
# write_file = glob.glob(model_ws+'/*.hob')
for n in np.arange(0,n_nodes).astype(str):
    folder = '/r'+ n.zfill(3)
    shutil.copy(write_file[0], model_ws+folder)
