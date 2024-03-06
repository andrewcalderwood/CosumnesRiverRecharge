from multiprocessing import Pool
import multiprocessing
import os
from os.path import join,basename,dirname
import glob
import subprocess
import time
import numpy as np

ext_dir = 'F:/WRDAPP'
c_dir = 'C:/WRDAPP'

if os.path.exists(ext_dir):
    loadpth = ext_dir 
elif os.path.exists(c_dir):
    loadpth = c_dir 

tprogs_id = '_no_conditioning'

loadpth = loadpth +'/GWFlowModel/Cosumnes/levee_setback/setback_distance_analysis/'
model_ws = loadpth+'Connec3d' + tprogs_id+'/'

gwfm_dir = 'C:/Users/ajcalder/Box/research_cosumnes/GWFlowModel'
mf_tprogs_dir = gwfm_dir+'/UPW_data/tprogs_final' + tprogs_id+'/'
tprogs_files = glob.glob(mf_tprogs_dir+'*')

def f(t):
    tprogs_line = np.loadtxt(tprogs_files[t])
    # convert any negatives representing input data to same value
    tprogs_arr = np.abs(np.reshape(tprogs_line, (320, 100,230)))
    conn_arr = np.zeros((320, 100,230))
    # new array where sand, gravel (1,2) are 1-coarse and sandy mud, mud are 0-fine
    conn_arr[(tprogs_arr == 1)|(tprogs_arr == 2)] = 1
    # convert to z,y, x order
    conn_line = np.reshape(conn_arr, (320*100*230) ) 

    dat_nam = 'r'+str(t).zfill(3)+'.DAT'
    np.savetxt(join(model_ws,  dat_nam), conn_line.astype(int))


def main():
    pool = Pool(processes=multiprocessing.cpu_count()-2)  # set the processes max number to the number of cpus
    result = pool.map(f, range(1))
    pool.close()
    pool.join()
    print(result)
    print('end')



if __name__ == "__main__":
    tic = time.time()
    main()
    toc = time.time()
    print('Total time: %.2f minutes' %((toc-tic)/60))