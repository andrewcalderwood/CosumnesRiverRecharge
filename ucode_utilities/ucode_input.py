"""
ucode_input module. 
Different functions for preparing data from MODFLOW (e.g. HOB output) into input files for UCODE
First iteration as a Module June 2022
Author: Andrew Calderwood
"""

import sys
import numpy as np
from numpy import ma

def get_n_nodes(n_params):
    ''' Returns number of cpu nodes that should be used for parallel processing'''
    import multiprocessing
    max_n_nodes = multiprocessing.cpu_count() - 1

    if n_params < max_n_nodes:
        n_nodes = n_params
    elif n_params >= max_n_nodes:
        n_nodes = max_n_nodes
    return(n_nodes)

def write_hob_jif_dat(model_ws, hobout, statflag=False):
    ''' Create the JIF and DAT files for the HOB package based on standard MODFLOW output
    and create a StandardFile for UCODE'''
    obsoutnames = hobout['OBSERVATION NAME']
    obs_vals = hobout['OBSERVED VALUE']

    header = 'jif @\n'+'StandardFile  1  1  '+str(len(obsoutnames))
    # obsoutnames.to_file(m.model_ws+'/MF.hob.out.jif', delimiter = '\s+', index = )
    np.savetxt(model_ws+'/MF.hob.out.jif', obsoutnames.values,
               fmt='%s', delimiter = r'\s+',header = header, comments = '')
    # ucode wants the observed value already written out to the obs_table
    # the simulated equivalent and obs name must be referenced 
    # in the hob.out.jif file
    # file1.writelines('  NROW='+str(len(obs_vals))+' NCOL=3 COLUMNLABELS \n')
    # file1.writelines('  ObsName GroupName ObsValue'+'\n')
    cols = ['OBSERVATION NAME','group','OBSERVED VALUE']
    stat=''
    if statflag==True:
        stat = ' Statistic StatFlag'
        cols +=['Statistic','StatFlag']

    header = 'BEGIN Observation_Data Table\n'+\
        'NROW='+str(len(obs_vals))+' NCOL='+str(len(cols))+' COLUMNLABELS \n'+\
        'ObsName GroupName ObsValue'+stat

    footer = 'End Observation_Data'
    # add column for group
    hobout['group'] = 'Heads'
    # get array of just strings
    hob_arr = hobout.loc[:,cols].values
    # pull out observed value and name of obs
    np.savetxt(model_ws+'/hob_obs_table.dat', hob_arr,
               fmt='%s', header = header, footer = footer, comments = '' )

def write_parallel(model_ws, n_nodes,exp_runtime):
    ''' Function to create the batch file that starts the runners needed
    for parallel processing with UCODE2014 and writes a text file that has the needed 
    code block to specify parallel computing in the UCODE main.in file
    PARAMETERS:
    model_ws = directory where parallel runs will exist
    n_nodes = number of parallel runs (based on processor count)
    exp_runtime = expected/maximum length of a model run in seconds for UCODE to expect
    '''
    exp_runtime = ' '+str(exp_runtime)
    # Open a file in write mode to open all runners
    f = open(model_ws+'/00_runner_all_ucode.bat', 'w')
    ft = open(model_ws+'/00_runner_all_ucode_table.txt', 'w')
    # write parallel runs
    #100
    ft.write('BEGIN Parallel_Runners Table\n')
    ft.write('  NROW='+str(n_nodes)+' NCOL=3 COLUMNLABELS\n')
    ft.write('  RunnerName RunnerDir RunTime\n')
    n= str(0)
    folder = 'r'+ n.zfill(3)+'\\'
    f.write('cd '+ folder+'\n')
    f.write('Start "Runner '+n.zfill(3)+'" /min runner'+'\n')
    ft.write('  Runner'+ n.zfill(3)+' '+folder+exp_runtime+'\n')
    for n in np.arange(1,n_nodes).astype(str):
        folder = 'r'+ n.zfill(3)+'\\'
        f.write('cd ..\\'+ folder+'\n')
        f.write('Start "Runner '+n.zfill(3)+'" /min runner'+'\n')
    # RunnerName, RunnerDir Run Time
        ft.write('  Runner'+ n.zfill(3)+' '+folder+exp_runtime+'\n')

    # Start "Runner 1"/min runner#
        # close batch file now that all is written
    f.close()
    ft.write('END Parallel_Runners\n')
    ft.close()