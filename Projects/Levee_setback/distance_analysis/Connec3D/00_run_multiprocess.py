from multiprocessing import Pool
import multiprocessing
import os
from os.path import basename, dirname, join
import glob
import subprocess
import time

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

def f(n):
    file = 'r'+ str(n).zfill(3)+'.PAR'
    print(file)
    runline = 'echo '+ file +' | '+'Connec3D_O3.exe'
    rv = subprocess.run(runline,  shell=True, check=True, capture_output=True, cwd = model_ws)
    return rv


def main():
    pool = Pool(processes=multiprocessing.cpu_count()-2)  # set the processes max number to the number of cpus
    result = pool.map(f, range(100))
    pool.close()
    pool.join()
    print(result)
    print('end')


if __name__ == "__main__":
    tic = time.time()
    main()
    toc = time.time()
    print('Total time: %.2f minutes' %((toc-tic)/60))