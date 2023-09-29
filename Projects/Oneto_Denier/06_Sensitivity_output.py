# +
# standard python utilities
import os
from os.path import join, basename,dirname, exists, expanduser
import sys
import glob
import pandas as pd
import numpy as np
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns

# standard geospatial python utilities
# import pyproj # for converting proj4string
# import shapely
import geopandas as gpd
# import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator


# +
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
    
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study
proj_dir = gwfm_dir + '/Oneto_Denier/'
dat_dir = proj_dir+'Stream_level_data/'

fig_dir = proj_dir+'/Streambed_seepage/figures/'
hob_dir = join(gwfm_dir, 'HOB_data')
sfr_dir = gwfm_dir+'/SFR_data/'



# -

sa_dir = join(fig_dir, 'sensitivity_analysis')
os.makedirs(sa_dir, exist_ok=True)


# +
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)

add_path(doc_dir+'/GitHub/flopy')
import flopy 
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)
from mf_utility import get_dates, get_layer_from_elev
from map_cln import gdf_bnds, plt_cln
# -

# ## Specify model

# +
loadpth =  'C:/WRDAPP/GWFlowModel/Cosumnes/Stream_seepage'
scenario=''
scenario='_no_reconnection'

upscale = 'upscale4x_'
model_nam = 'oneto_denier_'+upscale+'2014_2020'
model_ws = join(loadpth,model_nam+scenario) +'_ucode'
m_nam = 'MF_ucode'
# -

u_name = 'r11_full'
# u_name = 'r11_missing_flow'
ucode_file = join(model_ws, u_name)
sa_run_dir = join(sa_dir, u_name+scenario)
os.makedirs(sa_run_dir, exist_ok=True)


# +

def clean_ucode(filename):
    cleaned_columns = pd.read_csv(filename,header=None,nrows=1, 
                                  sep= '" "', engine='python').replace('"','', regex=True).values[0]
    df = pd.read_csv(filename,sep=r'\s+',engine='python',header=None,skiprows=1)
    df.columns = cleaned_columns
    return(df)


# -

# ss

# +
{'scgrp':'Grouped CSS', 'sc': 'CSS', 'sc_svd':'CSS by SVD'
}
ext = ['sc', 'scgrp']
n=0

filename = join(ucode_file, m_nam+'._'+ext[n])
print(filename)
ss = clean_ucode(filename)
ss['NORMALIZED COMPOSITE SCALED SENSITIVITY'] = ss.loc[:,'COMPOSITE SCALED SENSITIVITY'] / ss.loc[:,'COMPOSITE SCALED SENSITIVITY'].max()


# +
# ss

# +
# fig, ax = plt.subplots(2,1, sharex=True)
# ax_n = ax[0]
# ss.plot(x='PARAMETER NAME',y='COMPOSITE SCALED SENSITIVITY',kind='bar',legend=False, ax=ax_n)
# ax_n.set_yscale('log')
# ax_n = ax[1]
fig, ax_n = plt.subplots()
ss.plot(x='PARAMETER NAME',y='COMPOSITE SCALED SENSITIVITY',kind='bar',legend=False, ax=ax_n)

fig.supylabel('NORMALIZED COMPOSITE \n SCALED SENSITIVITY')
plt.savefig(join(sa_run_dir, 'noramlized_CSS.png'),bbox_inches='tight',dpi=600)

# -








