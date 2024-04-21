import os
from os.path import join, expanduser
import glob
from datetime import datetime
import shutil

# set up directories
usr_dir = expanduser('~')
gwfm_dir = join(usr_dir,'Box','research_cosumnes','GWFlowModel')
out_dir = join(gwfm_dir, 'Oneto_Denier', 'output','notebook_html')
# identify html files to move to Box for saving
files = glob.glob("*.html")
print(files)
date = datetime.today().strftime('%Y-%m-%d')
for f in files:
	shutil.move(f, join(out_dir, date+'_'+f))

# delete the ipynb associated with notebooks?