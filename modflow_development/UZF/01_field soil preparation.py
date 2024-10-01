# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
# standard python utilities
import os
from os.path import exists, join, basename, dirname
import sys
import glob
import pandas as pd
import numpy as np
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import geopandas as gpd
from rasterstats import zonal_stats
import richdem
import rasterio

# mapping utilities
import contextily as ctx
# import osmnx as ox # open street map
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm



# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
gwfm_dir
dis_dir = join(gwfm_dir, 'DIS_data')
uzf_dir = join(gwfm_dir,'UZF_data')

# %%
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')

m_domain = gpd.read_file(gwfm_dir+'/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp')


# %%
ag = gpd.read_file(join(uzf_dir, 'county_landuse','domain_ag_lu_2018.shp'))
# subset for relevant columns
ag = ag[['geom_id','name','irr_name','geometry']]

native = gpd.read_file(join(uzf_dir, 'county_landuse','domain_native_lu_2018.shp'))
# subset for relevant columns
native = native[['geom_id','name','irr_name','geometry']]
native.loc[native.irr_name.isna(), 'irr_name'] = 'no irrig'

# %%
fields = pd.concat((ag, native))

# %% [markdown]
# # Soil data
# Join gridded data to parcel data to calculate the soil parameters at the parcel scale for irrigated lands only.
#

# %%
uzf_dir = gwfm_dir+'/UZF_data'
soil_path = join(uzf_dir,'clean_soil_data')


# %%
gpd_soil = gpd.read_file(join(soil_path, 'cleaned_spatial_soil.shp'))
gpd_soil = gpd_soil.drop(columns=['AREASYMBOL', 'MUSYM']) # drop extra string columns

# %%
# only keep numeric columns that can be averaged (won't need texture class or HydGroup)
# soil_ag = gpd.overlay(gpd_soil, ag)
soil_ag = gpd.overlay(gpd_soil, fields)

# %%
# soil_ag[soil_ag.group=='native']
# soil_ag_cln[soil_ag_cln.group=='native']

# %%
# create a numeric equivalent hydgroup
HG_ref = pd.DataFrame({'num_hg':[1,2,3,4], 'hg':['A','B','C','D']})
soil_ag['Num_hydgroup'] = HG_ref.set_index('hg').loc[soil_ag.HydGroup, 'num_hg'].values

# %%
from scipy.stats import mode
# calculate average soil parameters
soil_ag_cln = soil_ag.drop(columns=['Texture','HydGroup'])
soil_ag_cln = soil_ag_cln.dissolve(['geom_id','name','irr_name'], np.mean)
# for objects use the mode
soil_ag_mode = soil_ag.dissolve(['geom_id','name','irr_name'], pd.Series.mode)

# add string cols back (texture, hydgroup)
soil_ag_cln = soil_ag_cln.join( soil_ag_mode[['Texture','HydGroup']]).reset_index()

# given the mean of the calculate HydGroup return to the nearest letter
soil_ag_cln.Num_hydgroup = soil_ag_cln.Num_hydgroup.round(0).astype(int) # clean up to nearest integer
soil_ag_cln.HydGroup = HG_ref.set_index('num_hg').loc[soil_ag_cln.Num_hydgroup, 'hg'].values


# %%
soil_ag_cln['field_area_m2'] = soil_ag_cln.geometry.area

# %%
soil_ag_cln.drop(columns='geometry').to_csv(join(soil_path, 'soil_for_ag_fields.csv'), index=False)

# %%
# identify the grid cells that will be connected to each ag field
grid_soil = gpd.overlay(grid_p, soil_ag_cln)

# if there are multiple ag field in a cell then the recharge/pumping needs to be scaled by the ag field area within the cell
grid_soil['cell_field_area_m2'] = grid_soil.area

# %%
# save ag field - grid cell connection and area
out_cols = ['row','column', 'geom_id','name','irr_name','cell_field_area_m2']
grid_soil[out_cols].to_csv(join(soil_path, 'ag_field_to_cell.csv'),index=False)

# %%
# there are some non-unique rows in the ag field dataset
soil_ag_cln.shape, ag.shape, ag.geom_id.unique().shape

# %%

# %% [markdown]
# ## Get field slope to assist curve numbers

# %%
raster_name = gwfm_dir+"/DEM_data/USGS_ten_meter_dem/modeldomain_10m_transformed.tif"

# load DEM for slope/aspect analysis
rd_dem = richdem.LoadGDAL(raster_name)

dem = rasterio.open(raster_name)
affine = dem.meta['transform']


# %%
def attr_stats(rd_dem, gdf, attr, folder):
    # calculate slope/aspect
    rd_attr = richdem.TerrainAttribute(rd_dem, attr)
    # takes several minutes
    zs_parcels = zonal_stats(gdf, rd_attr, affine=affine, nodata=rd_attr.no_data,
                stats=['min', 'max', 'mean', 'median', 'majority', 'std'])
    # zs_parcels = zonal_stats(parcels, raster=raster_name, stats=['min', 'max', 'mean', 'median', 'majority','std'])
    # convert to dataframe
    zs_df = pd.DataFrame(zs_parcels)
    # join zone stats of DEM to parcel data
    zs_df = gdf.join(zs_df)
    # save to shapefile
    zs_df.to_file(join(uzf_dir, folder, attr+'_statistics.shp'))
    return(zs_df)



# %%
# attr = 'slope_percentage'
for attr in ['slope_percentage']:
    slp = attr_stats(rd_dem, soil_ag_cln[['geom_id','name','irr_name','geometry']], 
                     attr, 'clean_soil_data/field_zonalstats')

# %% [markdown]
# ## Curve numbers

# %%
# load data for slope
zs_df = gpd.read_file(join(uzf_dir, 'clean_soil_data/field_zonalstats','slope_percentage_statistics.shp'))


# %%
soil_ag_CN = pd.read_csv(join(soil_path, 'soil_for_ag_fields.csv'))
# soil_ag_CN = soil_ag.copy()
soil_ag_CN['lu']='Cultivated land'
# assume that native lands are approximately like pasture for runoff, limited runoff
soil_ag_CN.loc[soil_ag_CN.irr_name=='no_irr', 'lu']='Pasture'

# %%
# overlay crops geometry so cell area can be taken into account
C_gpd = pd.merge(soil_ag_CN, zs_df[['geom_id','name','irr_name', 'mean']], how='left')
C_gpd.loc[C_gpd['mean'].isna(),'mean'] = 0 # if slope wasn't sampled then assume 0

# %%
# load curve number data and clean up
CN = pd.read_excel(join(gwfm_dir,'UZF_data','curve_numbers.xlsx'), comment='#')
CN = CN.rename(columns={'Cover type':'lu', 'Impervious':'impervious', 'Hydrologic Condition':'HydCond'})

CN_long = CN.melt(id_vars = ['lu','impervious','HydCond'], var_name='HydGroup', value_name='CN')
# columns for joining
CN_long = CN_long[['HydGroup','lu','impervious', 'HydCond', 'CN']]

# %%
# join spatial data to ID data for Curve Numbers
CN_out = C_gpd.join(CN_long.set_index(['HydGroup','lu']), on=['HydGroup','lu'], how='inner')

# for pasture there is an option of fair or poor hydrologic condition
# poor condition should be associated with slopes greater than 3%, less than 50% cover in hills
CN_pasture = CN_out[CN_out.lu=='Pasture']
hills_pasture = CN_pasture[(CN_pasture['mean']>=3)&(CN_pasture.HydCond=='Poor')]
flat_pasture = CN_pasture[(CN_pasture['mean']<3)&(CN_pasture.HydCond=='Fair')]
CN_pasture = pd.concat((hills_pasture, flat_pasture))
# add pasture back to out file
CN_out = CN_out[CN_out.lu!='Pasture']
CN_out = pd.concat((CN_out, CN_pasture))

# %%
CN_out.shape, soil_ag_CN.shape, C_gpd.shape

# %%
CN_out[['geom_id','name','irr_name','CN']].to_csv(join(soil_path, 'ag_field_CN.csv'))
CN_out.to_csv(join(soil_path, 'ag_field_properties_all.csv'))
