from os.path import join, basename, dirname, exists
import pandas as pd

from multiprocessing import Pool
import multiprocessing
import os
import subprocess
import time


run_dir = 'C:/WRDAPP/GWFlowModel'
run_dir = 'F://WRDAPP/GWFlowModel'
# run_dir = 'D://WRDAPP/GWFlowModel'

loadpth = run_dir +'/Cosumnes/Economic/'

folders = list()

model_nam = 'input_write_2014_2022'
for scen in ['R200']:
# scen = 'R200'
    # it is probably better to create a slightly different file name then to copy these over for a set scenario
    econ_model_ws = join(loadpth, model_nam+'_'+scen, 'crop_modflow')
    
    all_run_dates = pd.read_csv(join(econ_model_ws, 'all_run_dates.csv'))
    
    econ_model_ws = join(loadpth, model_nam+'_'+scen, 'crop_modflow')
    # get subfolder basenames
    folders_run = all_run_dates.date[:-1]
    # set complete filepaths to have more flexibility on the command
    folders_run = (econ_model_ws+'/'+folders_run).tolist()

    folders = folders + folders_run
# print(folders, sep='\n')

# %%
def f(folder):
    print(f"Starting in {folder}")
    try:
        rv = subprocess.run('mf-owhm.exe MF.nam > MF_log.txt', shell=True, check=True, capture_output=True, cwd = folder)
        print(f"Finished {folder}.")
        # print(f"Output: {rv.stdout}.")
        return(rv)
    except subprocess.CalledProcessError as e:
        print(f"Error in {folder}: {e.stderr}")
    except Exception as e:
        print(f"An unexpected error occurred in {folder}: {e}")



def main(folders):
    # set the processes max number to the number of cpus/2 to avoid overtaxing
    pool = Pool(processes=int(multiprocessing.cpu_count()/2))  
    result = pool.map(f, folders)
    pool.close()
    pool.join()
    print(result)
    print('end')

# %%
if __name__ == "__main__":
    tic = time.time()
    main(folders)
    toc = time.time()
    print('Total time: %.2f minutes' %((toc-tic)/60))


# %%
import subprocess
import concurrent.futures
import os

def run_script_in_dir(directory, script_command):
    """
    Function to run a shell command in a specific directory.

    Args:
        directory (str): The path to the working directory.
        script_command (list or str): The command to run (list is safer).
    """
    print(f"Starting in {directory}")
    try:
        # Use subprocess.run with the cwd argument
        # 'shell=True' might be needed for some batch file specifics or chained commands (e.g., 'cd X && mycommand')
        # However, passing the command as a list and using 'cwd' is generally preferred for safety.
        # For a simple batch file, you might use ['my_script.bat']
        result = subprocess.run(script_command, cwd=directory, check=True, capture_output=True, text=True, shell=True)
        print(f"Finished {directory}. Output: {result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"Error in {directory}: {e.stderr}")
    except Exception as e:
        print(f"An unexpected error occurred in {directory}: {e}")

# %%

# command = ['mf-owhm.exe MF.nam']*len(folders)
# # specify run folder and corresponding command
# tasks = list(zip(folders, command))

# # Max number of processes to run concurrently (e.g., number of CPU cores)
# MAX_WORKERS = 4 

# # Use ProcessPoolExecutor to run tasks in parallel
# with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
#     # Submit tasks to the executor
#     # map(function, *iterables)
#     # The 'map' method applies the function to each item in the iterables
#     executor.map(run_script_in_dir, [task[0] for task in tasks], [task[1] for task in tasks])

# print("All processes have finished.")