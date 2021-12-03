from multiprocessing import Pool
import multiprocessing
import os
import subprocess
import time

ext_dir = 'F:/WRDAPP'
c_dir = 'C:/WRDAPP'

if os.path.exists(ext_dir):
    loadpth = ext_dir 
elif os.path.exists(c_dir):
    loadpth = c_dir 

loadpth = loadpth +'/GWFlowModel/Cosumnes/levee_setback/setback_distance_analysis/'
model_ws = loadpth+'Permeameter_for_velocity'


def f(n):
    folder = '/realization'+ str(n).zfill(3)+'/'
    rv = subprocess.run('mf2005.exe MF.nam', shell=True, check=True, capture_output=True, cwd = model_ws + folder)
    return rv


def main():
    pool = Pool(processes=multiprocessing.cpu_count())  # set the processes max number to the number of cpus
    result = pool.map(f, range(100))
    pool.close()
    pool.join()
    print(result)
    print('end')



if __name__ == "__main__":
    tic = time.time()
    main()
    toc = time.time()
    print('Total time: %.2f minutes' (toc-tic)/60)