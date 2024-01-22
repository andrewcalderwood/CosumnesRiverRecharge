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
# # Expand floodplain at Oneto-Denier
# Only inundating the inner floodplain at Oneto-Denier led to not a large amount of baseflow being developed. How much land would need to be added to create more continuous baseflow. Or is it a constraint of the model?

# %% editable=true slideshow={"slide_type": ""}
# standard python utilities
import os
import sys
from os.path import basename, dirname, join, exists, expanduser
import glob
import time

import pandas as pd
import numpy as np

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import geopandas as gpd
import rasterio


# %% editable=true slideshow={"slide_type": ""}
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study
proj_dir = gwfm_dir + '/Oneto_Denier/'
gis_dir = join(proj_dir, 'GIS')
dat_dir = proj_dir+'Stream_level_data/'

sfr_dir = gwfm_dir+'/SFR_data/'
uzf_dir = gwfm_dir+'/UZF_data/'


# %% editable=true slideshow={"slide_type": ""}
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
from mf_utility import get_layer_from_elev, param_load

from map_cln import gdf_bnds, plt_cln

# %%
ext_dir = 'F:/WRDAPP'
c_dir = 'C:/WRDAPP'

if os.path.exists(ext_dir):
    loadpth = ext_dir 
elif os.path.exists(c_dir):
    loadpth = c_dir 
upscale=8   
loadpth += '/GWFlowModel/Cosumnes/Stream_seepage/'
model_nam = 'oneto_denier_upscale4x_2014_2020'
# model_nam = 'oneto_denier'

model_ws = loadpth+ model_nam 
#+'_'+ str(strt_date.year)+'_'+str(end_date.year)
# if scenario != '':
#     model_ws += '_' + scenario

# %%
load_only = ['DIS','BAS6','SFR', 'UPW']
m = flopy.modflow.Modflow.load(model_ws+'/MF.nam', load_only=load_only)

# %%
model_grp = 'inset_oneto_denier'
grid_dir = join(gwfm_dir, 'DIS_data/streambed_seepage/grid')
grid_fn = join(grid_dir, model_grp,'rm_only_grid.shp')
grid_p = gpd.read_file(grid_fn)
grid_p.crs='epsg:32610'

delr, delc = 100,100

# %%
# load other model parameters
ucode_dir = join(gwfm_dir, 'UCODE')
bc_params = param_load(model_ws, ucode_dir, 'BC_scaling.csv')
# bc_params = pd.read_csv(join(model_ws,'BC_scaling.csv'))
bc_params = bc_params.set_index('ParamName')

# %%
vka = m.upw.vka.array
botm = m.dis.botm.array
nlay, nrow,ncol = m.dis.botm.shape
dem_data = np.loadtxt(join(gis_dir, 'local_subset_dem_52_9_200m_mean.tsv'))


# %% [markdown]
# ## Prepare Lake bathymetry

# %%
lak_extent = gpd.read_file(join(gis_dir,'floodplain_expansion.shp'))


# %%
order = 7
# subset for only the zones that are already floodplain-like is <=4
# # subset to include everything but the zones to the south is <=5
# don't subset to include all regions
lak_extent = lak_extent.loc[lak_extent.order<=order]


# %%
fn = join(gis_dir,"extended_floodplain_"+str(order)+"_crop.tif")
if not exists(fn):
    # create clipped raster of just lake area
    dem_dir = join(gwfm_dir,'DEM_data')
    raster_name = dem_dir+'/mwt_peri_2_3.tif/mwt_peri_2_3_clipped.tif'
    import rasterio.mask
    with rasterio.open(raster_name) as src:
        out_image, out_transform = rasterio.mask.mask(src, lak_extent.geometry.values, crop=True)
        out_meta = src.meta
    # write output
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    with rasterio.open(fn, "w", **out_meta) as dest:
        dest.write(out_image)

# prepare bathymetry file
lakeRst = rasterio.open(fn)
lakeBottom = lakeRst.read(1)
noDataValue = np.copy(lakeBottom[0,0])
#replace value for np.nan
lakeBottom[lakeBottom==noDataValue]= np.nan

# get raster minimum and maximum 
minElev = np.nanmin(lakeBottom)
maxElev = np.nanmax(lakeBottom)
print('Min bottom elevation %.2f m., max bottom elevation %.2f m.'%(minElev,maxElev))

# steps for calculation
nSteps = 151
# lake bottom elevation intervals
elevSteps = np.round(np.linspace(minElev,maxElev,nSteps),2)

# definition of volume function
def calculateVol_A(elevStep,elevDem,lakeRst, conv=1):
    tempDem = elevStep - elevDem[elevDem<elevStep]
    tempArea = len(tempDem)*lakeRst.res[0]*conv*lakeRst.res[1]*conv
    tempVol = tempDem.sum()*lakeRst.res[0]*conv*lakeRst.res[1]*conv
    return(tempVol, tempArea)
# calculate volumes, areas for each elevation
volArray = [0]
saArray = [0]
for elev in elevSteps[1:]:
    tempVol,tempArea = calculateVol_A(elev,lakeBottom,lakeRst)
    volArray.append(tempVol)
    saArray.append(tempArea)

# print("Lake bottom elevations %s"%elevSteps)
volArrayMCM = round(volArray[-1]/1000000,2) 
print("Lake volume in million of cubic meters %s"%volArrayMCM)

# %%
# lak_buf = lak_extent[['OID_','geometry']].copy()
# lak_buf.geometry = lak_buf.buffer(10)
lak_grid = gpd.overlay(grid_p, lak_extent[['order','name','geometry']], how='intersection')
# check if more than 50% of cell is covered by the lake, avoid conflicts with sfr
lak_grid = lak_grid.loc[lak_grid.geometry.area > (delr*delc*0.5)]


# %% [markdown]
# # LAK Package

# %% [markdown]
# When using a rectangular grid for TPROGs the lake package must be more carefully defined in terms of layering.
#
# - expanded lake is to allow for more recharge
# - the simplest version is to add it to the existing singular defined lake as quick test
# - the advanced version would be to make the expanded lakes connected to first lake to allow more varied depths

# %%
sfr_rch = pd.DataFrame(m.sfr.reach_data)
# sfr_rch

# %%
# remove lake grid cells that overlap with sfr cells
lak_grid_cln = lak_grid.join(sfr_rch.set_index(['i','j'])[['iseg']], on=['row','column'])
lak_grid_cln = lak_grid_cln[lak_grid_cln.iseg.isna()]
lak_row, lak_col = lak_grid_cln.row.values-1, lak_grid_cln.column.values-1
# find the layers above the dem
lak_k = get_layer_from_elev(dem_data[lak_row, lak_col], botm[:,lak_row, lak_col], m.dis.nlay)
# the lake should include the layer below th dem as well, fix issue with min lake elev in MF
lak_k += 1

# Set empty array of zeros for nonlake cells
lakarr = np.zeros((nlay, nrow,ncol))
# Each lake is given a different integer, and needs to be specified depending on the layer
# may need to decide if lake should be in more than 1 layer
# lakarr[lak_k, lak_row, lak_col] = 1
# for lakarr I think I may need to assign all cells from and above the lake because the documentation
# example shows the layers above as well
for n in np.arange(0,len(lak_row)):
    lakarr[:lak_k[n], lak_row[n], lak_col[n]] = 1
    
# set Ksat same as vertical conductivity, 
lkbd_thick = 2
lkbd_K = np.copy(vka)
lkbd_K[lakarr==0] = 0 # where lake cells don't exist set K as 0
# leakance is K/lakebed thickness, reduce by 1/10 for cloggin
# bdlknc = (lkbd_K/lkbd_thick)/10 #, accounted for in seep_vka
bdlknc = (lkbd_K/lkbd_thick)/bc_params.loc['bdlknc_scale', 'StartValue']



# %%
lak_active = (np.sum(lakarr,axis=0)>0) # cells where lake is active

# %%
# Exactly 151 lines must be included within each lake bathymetry input file and each line must contain 1 value 
#  of lake stage (elevation), volume, and area (3 numbers per line) if the keyword “TABLEINPUT” is specified in item 1a.
# A separate file is required for each lake. 
# initial lake stage should be dry (below lake bottom)
# stages = minElev - lkbd_thick - 0.1 # causes lake to remain dry for entire simulation
stages = minElev +0.01

# (ssmn, ssmx) max and min stage of each lake for steady state solution, there is a stage range for each lake
# so double array is necessary
stage_range = [[minElev, maxElev]]

# lake stage (elevation), volume, and area (3 numbers per line)
bathtxt = np.column_stack((elevSteps, volArray, saArray))

## Need to specify flux data
# Dict of lists keyed by stress period. The list for each stress period is a list of lists,
# with each list containing the variables PRCPLK EVAPLK RNF WTHDRW [SSMN] [SSMX] from the documentation.
flux_data = {}
flux_data[0] = {0:[0,0,0,0]} # default to no additional fluxes
# if ss_bool == True:
#     flux_data[0] = {0:[precip_lake.mean(), et_lake.mean(),0,0]}
# for j in np.arange(time_tr0, nper):
#     flux_data[j] = {0: [precip_lake[j-1], et_lake[j-1], 0, 0] } 


# if scenario != 'no_reconnection':
    # 1 1000 1E-5 0.02 - taken from mt shasta
    # filler value for bdlknc until soil map data is loaded by uzf
lak = flopy.modflow.ModflowLak(model = m, lakarr = lakarr, bdlknc = bdlknc,  stages=stages, 
                               stage_range=stage_range, flux_data = flux_data,
                               theta = 1, nssitr = 1000, sscncr = 1E-5, surfdepth = 0.02, # take from Shasta model
                               tabdata= True, tab_files='MF.bath', tab_units=[57],ipakcb=55)

lak.options = ['TABLEINPUT']
# # need to reset tabdata as True before writing output for LAK
lak.tabdata = True


# %%
#lak.write_file()


# %%
# if scenario != 'no_reconnection':
flopy.modflow.mfaddoutsidefile(model = m, name = 'DATA',extension = 'bath',unitnumber = 57)


# %%
# # numgage is total number of gages
# # gage_data (list, or array), includes 2 to 3 entries (LAKE UNIT (OUTTYPE)) for each LAK entry
# #  4 entries (GAGESEG< GAGERCH, UNIT, OUTTYPE) for each SFR package entry
# # gage_data = rm_xs[['iseg','reach','unit', 'outtype']].values.tolist()
# gage_file = ['MF.gage']
# gag_out_files = ('MF_gage_' + XSg.index.astype(str) +'.go').values.tolist()

# lak_gage_data = [[-1, -37, 1]]
# gage_file = ['MF.gage']
# gag_out_files = ['MF_lak.go']
# gag = flopy.modflow.ModflowGage(model=m,numgage= 1,gage_data=lak_gage_data, 
#                                 filenames =gage_file+gag_out_files)

# %%
ibound = m.bas6.ibound.array.copy()
# #lake cells must be set as inactive
ibound[lakarr>0] = 0

# %%
# ibound < 0 is constant head
# ibound = 0 is inactive cell
# ibound > 0 is active cell
# strt is array of starting heads
# add option: STOPERROR 0.01 to reduce percent error when OWHM stops model
# if solver criteria are not met, the model will continue if model percent error is less than stoperror
bas = flopy.modflow.ModflowBas(model = m, ibound=ibound, strt = m.bas6.strt.array, stoper = None) #

# %%
# update model work space and write input
new_ws = join(loadpth, 'oneto_denier_upscale4x_2014_2020_floodplain'+str(order))
m.model_ws = new_ws
m.write_input()
np.savetxt(new_ws+'/MF.bath', bathtxt, delimiter = '\t')


# %%
# # copy over other relevant files
# # copy mf files except cbc and hds
mf_files = pd.Series(glob.glob(model_ws+'/MF.*'))
pks_rem = 'cbc|hds|list|.hob.out|upw|sfr|lak|bath|chk|bas'
mf_files = mf_files[~mf_files.str.contains(pks_rem).values].tolist()

files = mf_files
files

# %%
import shutil
for f in files:
    shutil.copy(f, new_ws)

# %%
