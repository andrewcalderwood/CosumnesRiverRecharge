# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Prepare any necessary reference files for the MF model updates

# %%
# standard python utilities
import os
from os.path import basename, dirname, join, exists, expanduser
import sys
from importlib import reload
import glob

import pandas as pd
import numpy as np
import time

# standard python plotting utilities
# import seaborn as sns
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates

# standard geospatial python utilities
import geopandas as gpd
import rasterio

# mapping utilities
# import contextily as ctx
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as fm
# from matplotlib.ticker import MaxNLocator


# %%
doc_dir = join(expanduser('~'), 'Documents')

git_dir = join(doc_dir, 'GitHub')
## Set up directory referencing
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

bas_dir = join(gwfm_dir, 'BAS6')
# proj_dir = join(gwfm_dir,'Projects','EcoFIP')
# gis_dir = join(proj_dir,'GIS')
plt_dir = join(proj_dir,'figures/')

# %%
lwa_dir = join(usr_dir, 'LWA Dropbox','01_Project-Teams')
proj_dir = join(lwa_dir, '669.03 - DWR Cosumnes Floodplain Recharge')
main_concept_dir = join(proj_dir, 'Concepts')
gis_dir = join(main_concept_dir,'GIS')
concept_dir = join(concept_dir,'Blodgett_Dam')


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

# from mf_utility import get_layer_from_elev
from map_cln import gdf_bnds, plt_cln
from mf_utility import get_dates

# %% [markdown]
# # Preprocessing pumping related files

# %%
# load dataframe to reference spatial extent of fields
fields = gpd.read_file(join(gwfm_dir, 'UZF_data', 'county_landuse','domain_ag_lu_2018.shp'))
fields = fields[['geom_id', 'name','geometry']]
fields['full_area'] = fields.area
# load dataframe to reference fields to wells
well2field = gpd.read_file(join(gwfm_dir, 'WEL_data','ag_fields_to_wells','ag_fields_to_wells.shp'))
well2field = well2field[['geom_id','row','column']]

# merge the field with the corresponding well location
fields = fields.merge(well2field)

fields.to_file(join(gis_dir, 'modflow_related','ag_fields_with_well_loc.shp'))


# %% [markdown]
# # Stream/floodplain update
