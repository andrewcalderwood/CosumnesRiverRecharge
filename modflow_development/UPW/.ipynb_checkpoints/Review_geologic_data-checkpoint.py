# +
# standard python utilities
import os
import sys
from os.path import basename, dirname, join, exists, expanduser
import glob
# import time

import pandas as pd
import numpy as np
import numpy.ma as ma
from scipy.stats import gmean

# standard python plotting utilities
# import matplotlib as mpl
# import matplotlib.pyplot as plt

# # standard geospatial python utilities
# # import pyproj # for converting proj4string
# import geopandas as gpd
# from osgeo import gdal
# import rasterio

# # mapping utilities
# import contextily as ctx
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as fm

# +
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

gel_dir = gwfm_dir+'/UPW_data/'
uzf_dir = gwfm_dir+'/UZF_data/'


# +
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
tprogs_fxn_dir = doc_dir+'/GitHub/CosumnesRiverRecharge/tprogs_utilities'
add_path(tprogs_fxn_dir)
import tprogs_cleaning as tc

from importlib import reload
reload(tc)
tprogs_info = [80, -80, 320]
# -

mf_tprogs_dir = gwfm_dir+'/UPW_data/tprogs_final/'
tprogs_files = glob.glob(mf_tprogs_dir+'*')


(nrow_p, ncol_p) = (100, 230)

# +

dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_52_9_200m_mean.tsv')
# check typical extent of model top
np.quantile(dem_data, [0,0.25,0.5,0.75,0.9,0.95, 1])
# -

t=0
tprogs_line = np.loadtxt(tprogs_files[t])
# filter elevation by regional model
# masked_tprogs= tc.tprogs_cut_elev(tprogs_line, dem_data_p, tprogs_info)
# if want to keep full hk, vka then don't crop elevation
masked_tprogs= tc.tprogs_cut_elev(tprogs_line, np.full((nrow_p,ncol_p),80), tprogs_info)

# identify conditioning well locations
cond_locs = pd.DataFrame(np.transpose(np.where(masked_tprogs<0)), columns=['k','i','j'])

cond_locs.k.hist()
