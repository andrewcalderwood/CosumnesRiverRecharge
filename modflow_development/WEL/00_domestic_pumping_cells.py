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

# %% [markdown]
# Goal: Preprocess the WCR information further based on well age, depth and location to prepare a general dataset useful for the Cosumnes modeling work. This was originally in the main model code, but was set aside to reduce clutter.

# %%
import os
from os.path import basename, dirname, exists, join
import glob

import numpy as np
import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt

# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
gwfm_dir
wel_dir = join(gwfm_dir, 'WEL_Data')

# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')

# %%
m_domain = gpd.GeoDataFrame(pd.DataFrame([0], columns=['id']),geometry = [grid_p.unary_union], crs=grid_p.crs)


# %%
nrow = grid_p.row.max()
ncol = grid_p.column.max()

# %%
wells_grid = pd.read_csv(wel_dir+'/wells_grid.csv')
wells_grid = gpd.GeoDataFrame(wells_grid, 
                              geometry= gpd.points_from_xy(wells_grid.easting, wells_grid.northing), crs = m_domain.crs)
domestic_pts = wells_grid[wells_grid.Simple_type=='domestic']

# %%
domestic_well_depth_arr = np.loadtxt(wel_dir+'/domestic_well_depth_arr.tsv', delimiter='\t')

# %% [markdown]
# # Identify the total extent of domestic pumping
# Review of residential land use against domestic pumping well point locations found that the points underestimate the total number of ag-residential parcels. This will result in an underestimate of pumping which is most prominent south of the Cosumnes River where it leads to largely oversimulated groundwater elevations. Parcel data is neded to subdivide the land use classification polygons. 

# %%
uzf_dir = join(gwfm_dir,'UZF_Data')

# %% [markdown]
# Prepare parcel data for model domain.

# %%
f_out = join(uzf_dir, 'Parcels','domain_parcels.shp')
if exists(f_out):
    parcels = gpd.read_file(f_out)
else:
    # very slow to run because of large file size
    fp = glob.glob(join(uzf_dir,'Parcels', 'Parcels - *'))
    parcels_all = gpd.GeoDataFrame()
    for f in fp:
        parcels_add = gpd.read_file(join(f, 'Parcels.shp'))
        parcels_all = pd.concat((parcels_all, parcels_add.to_crs(parcels.crs)))
    # crop to domain and save file
    parcels = gpd.overlay(parcels_all, m_domain.to_crs(parcels.crs)).to_crs(m_domain.crs)
    parcels.to_file(join(uzf_dir, 'Parcels','domain_parcels.shp'))

# %%
# land use from county data (Sac:2015, SJ:2018)
lu = gpd.read_file(join(uzf_dir, 'county_landuse','domain_urban_lu_2018.shp'))

# %% [markdown]
# Upon review of parcels and land use in QGIS, the rural parcels are not always totally covered by the urban land use class as the classification sometimes only covers the footprint of the home so the area of a parcel is underestimated.  
# Also regions like Rancho Murieta are not classified as urban when they are, and there is no public well identified in that area from WCRs so it would be better practice to assume they are drawing water from each parcel even though they are likely on one or two wells.
#
# **By not cropping rural parcels by land use classification the area of ag-residential doubles!**

# %%
# break up the landuse data into individual parcels
# need to use a spatial join to keep the area of each parcel as the land use doesn't fully cover parcels sometimes
lu_parcels = gpd.sjoin(parcels, lu)

# calculate area of the parcel
lu_parcels['area_m2'] = lu_parcels.geometry.area.copy()
lu_parcels['area_acres'] = lu_parcels.area_m2*(1/0.3048**2)/43560
lu_parcels = lu_parcels[['LU_GENERAL','LU_SPECIF','name','LU_DETAIL', 'CITY','area_m2','area_acres', 'geometry']]


# %%
# divide parcels into cities for aggregating numbers later
# cities = parcels.dissolve('CITY')

# %%
# with parcels LU_GENERAL 90% are covered by Residential when cropped to county landuse
perc_res1 = lu_parcels[lu_parcels.LU_GENERAL.isin(['Residential'])].shape[0]/lu_parcels.shape[0]
# 93.4% coverage using general residential str in LU name
perc_res2 = lu_parcels[lu_parcels.name.str.contains('residential')].shape[0]/lu_parcels.shape[0]
perc_res1, perc_res2

# %%
# filter to those only classified as residential, not urban which would be covered by municipal
# lu_res = lu_parcels[lu_parcels.name.str.contains('residential')]
lu_res = lu_parcels[lu_parcels.LU_GENERAL.str.contains('Residential', na=False)]

# %%
lu_res.LU_SPECIF.unique(), lu_res.LU_DETAIL.unique()
# Rural tag under LU_DETAIL
lu_rural = lu_parcels[lu_parcels.LU_DETAIL.str.contains('Rural', na=False)]


# %%
# lu_muni = lu_parcels[lu_parcels.LU_DETAIL.str.contains('Family|Neighborhood|Story', na=False)]
# lu_parcels.LU_DETAIL.unique()
lu_muni= lu_res[~lu_res.LU_DETAIL.str.contains('Rural', na=False)]
muni_pts = wells_grid[wells_grid.Simple_type=='public']


# %% [markdown]
# The number of parcels is much greater than the number of domestic wells because parcels in the Elk Grove/Galt are connected to municipal water systems and are accounted for in the municipal well pumping. Alternatively, I could remove the pumping based on domestic wells and municipal wells as is and instead estimate all residential pumping based on parcel size and use. 
#
# - The column LU_DETAIL from the parcel database is useful to identify rural residential parcels (those dependent on domestic wells). Water use at these locations should be about 2 AF/year because of high outdoor water use and the daily pumping should be scaled by the percent the daily ET is of the year after a minimum of 150 gpd for household use.  
#     - The LU_DETAIL 'Rural' caegories all showed only classificiations as LU_GENERAL:'Residential', LU_SPECIF:'Single-Family'  
# - This leaves the choice to update municipal pumping based on parcels or continue guessing the number of connections, realistically it should be updated to water use at each parcel since we don't know the number of connections or pumping data. For municipal locations water use should be kept at a constant 150 gpd for household use.

# %%
lu_muni.shape[0]*0.5, muni_pts.shape[0]*5000*365/7.48/43560, 0.5*43560*7.48/365

# %%
fig,ax= plt.subplots(1,2)
# lu_res.plot(ax=ax)
lu_muni.plot(ax=ax[1], color='yellow')
muni_pts_buf = muni_pts.copy()
muni_pts_buf.geometry = muni_pts_buf.buffer(1E3)
muni_pts_buf.plot(ax=ax[1], color='purple', alpha=0.4)

lu_rural.plot(ax=ax[0],color='green')
domestic_pts.plot(ax=ax[0], markersize=0.3, color='red')

# %%
lu_rural.shape[0]/domestic_pts.shape[0]
# using parcels nearly triples the number of ag-res which is what we expected

# %% [markdown]
# Each ag-residential parcel should sum on annual scale to 2AF, with a minimum daily value of 150 gpd
#

# %%

# spatial join with WCR points to keep information on well age and depth
lu_rural_out = gpd.sjoin(lu_rural,
                         domestic_pts[['WCRNumber','depth_m','geometry']], 
                         how='left').drop(columns=['index_right'])
# join parcels to the model grid assuming the well is in the centroid of the parcel
lu_rural_out.geometry = lu_rural_out.geometry.centroid
# rural_grid = gpd.sjoin(lu_rural_out, grid_p).drop(columns='index_right')
rural_grid = gpd.sjoin(grid_p, lu_rural_out).drop(columns='index_right')
# drop extra large parcels that shouldn't be considered residential, 20+ acres
rural_grid = rural_grid[rural_grid.area_acres < rural_grid.area_acres.quantile(0.99)]

# %% [markdown]
# Groupby parcels into area categories based on calculated area. By irrigating, it more generally means that water usage is scaled to property size.
# - 0-1.5 acres, at minimum each parcel should use about 2AF/year, so any p
# - 1.5-5, assume landowners are irrigating 80% of their parcel
# - 5-10 assume landowners are irrigating 60% of their parcel
# - 10-20 assume landowners are irrigating 40% of their parcel

# %%
rural_grid['usage_scale'] = 1/rural_grid.area_acres # default is 1 for area less than 1.5 acres
rural_grid.loc[rural_grid.area_acres>1.5,'usage_scale'] = 0.8
rural_grid.loc[rural_grid.area_acres>5,'usage_scale'] = 0.6
rural_grid.loc[rural_grid.area_acres>10,'usage_scale'] = 0.4
# calculated scaled area
rural_grid['used_area_acres'] = rural_grid.usage_scale * rural_grid.area_acres


# %%
# rural_grid.dissolve('CITY', aggfunc='sum').plot('used_area_acres', legend=True)
rural_grid.dissolve('CITY', aggfunc='sum').plot('area_acres', legend=True)

# %%
# fill in depth of wells based on interpolated domestic well depths
rural_grid['fill_depth_m'] = domestic_well_depth_arr[rural_grid.row-1, rural_grid.column-1]
rural_grid.loc[~rural_grid.depth_m.isna(),'fill_depth_m'] = rural_grid.loc[~rural_grid.depth_m.isna(), 'fill_depth_m']

# %%
rural_grid.drop(columns=['geometry']).to_csv(join(wel_dir, 'ag_res_parcel_domestic_wells.csv'))
