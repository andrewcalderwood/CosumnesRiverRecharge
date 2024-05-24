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
from os.path import join, exists, dirname, basename
import sys
import glob
import pandas as pd
import numpy as np
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
# import pyproj # for converting proj4string
# import shapely
import geopandas as gpd
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
    
# dir of all gwfm data
gwfm_dir = join(dirname(doc_dir),'Box/research_cosumnes/GWFlowModel')
dis_dir = join(gwfm_dir,'DIS_data')
print(gwfm_dir)

flopy_dir = doc_dir+'/GitHub/flopy'
if flopy_dir not in sys.path:
    sys.path.insert(0, flopy_dir)
import flopy 


# %%
from flopy.utils import Raster # for ET resampling


# %%
# New model domain 52.9 deg
m_domain = gpd.read_file(join(dis_dir,'NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp'))


# %%
xll, yll = list(m_domain.geometry.values[0].exterior.coords)[0]
#Maribeth's model parameters, had to switch nrow and ncol due to her issue in xul, yul
nrow=100
ncol=230
delr=np.repeat(200,ncol)
delc=np.repeat(200,nrow)
rotation=52.9
modelgrid = flopy.discretization.StructuredGrid(xoff=xll, yoff=yll, proj4='EPSG:32610', angrot=rotation,
                                   delr=delr, delc=delc, nrow=nrow,ncol=ncol)

# %% [markdown]
# ### Crop coefficients and ETo

# %%

uzf_dir = join(gwfm_dir,'UZF_data')

crop_path = join(uzf_dir,'county_landuse')
crop_shp_names = glob.glob(crop_path+'/*.shp')
# crop_shp_names


# %% [markdown]
# ### Reference crop class letter to the name

# %%
# land use referencing for DWR data (class, subclass, irrigation)
lu_class = pd.read_excel(join(uzf_dir,'DWR_landuse_ref.xlsx'), sheet_name = 'class' )
# clean up file, needed before saved over
# lu_class.name = lu_class.name.str.split('(', expand=True)[0]
# lu_class.name = lu_class.name.str.strip()

# files already clean
lu_subclass = pd.read_excel(join(uzf_dir,'DWR_landuse_ref.xlsx'), sheet_name = 'subclass' )
lu_irrig = pd.read_excel(join(uzf_dir,'DWR_landuse_ref.xlsx'), sheet_name = 'irrigation' )
lu_irrig.code = lu_irrig.code.str.strip()
lu_eff = pd.read_excel(join(uzf_dir,'DWR_landuse_ref.xlsx'), sheet_name = 'irr_efficiency', comment='#')
# calculate average irrigation efficiency from the range
lu_eff['Avg_eff'] = (lu_eff.Low_eff + lu_eff.High_eff)/2
# add irrigation efficiencies to type list
lu_irrig = lu_irrig.join(lu_eff.set_index('Irrigation Method'), on='irr_eff_name')

# # join with subclass file so if subclass doesn't exist then the bulk class is used
lu_class = pd.concat((lu_class,lu_subclass))

# %% [markdown]
# ### Prepare spatial data
# Due to the location, there are multiple county datasets required for the model domain. Additionally there are different datasets for certain years which need to be joined, these will later be filtered depending on the year of climate data.  
# The spatial data will be joined to the model grid before joining to the time series data because the model grid will be used to dissolve the data to each node.  
# \* the majority of parcels are single cropped, but still need to resolve those that are not by removing or accounting for crop sequences. For single crops the value to the cell is 100%, for intercropped need to check the percentage.
# - Now with the parcel data it is easier to spatial join the shapefile to the location of WCRs in the domain and filter to one or two to apply the pumping.

# %%
# join all the crop shapefile together
lu_all = gpd.GeoDataFrame()
# lu_yr = gpd.read_file(glob.glob(crop_path+'/*Sacramento2000.shp')[0])
for n in np.arange(0, len(crop_shp_names)):
    # load land use file
    lu_yr = gpd.read_file(crop_shp_names[n])
    # crop to model domain
    lu_yr = gpd.overlay(lu_yr, m_domain.to_crs(lu_yr.crs))
    # add id to use when pivoting dataframe
    lu_yr['geom_id'] = np.arange(0, len(lu_yr))
    # append to joint dataframe
    lu_all = pd.concat((lu_all, lu_yr.to_crs(m_domain.crs)))


# %%
def crop_number_split(lu_all):
    """
    Take land use data from county/DWR surveys and change from wide to long format,
    e.g., CLASS1, CLASS2, CLASS3, PCNT1, PCNT2  becomes two columns
    variable column (var_nam) identifies CLASS, IRR_TYP, PCNT
    value column has the corresponding value
    the variable column is split up into the main name and the crop number (crop_num)
    
    """
    # find columns that use numbers as identifiers
    melt_cols = lu_all.columns[lu_all.columns.str.contains(fr'\D\d')]
    melt_cols = melt_cols[~melt_cols.str.contains('AREA')]
    id_cols = lu_all.columns[~lu_all.columns.isin(melt_cols)]
    # melt the data to create clean columns broken up by crop numbers
    lu_long = lu_all.melt(value_vars = melt_cols,
               id_vars = id_cols
                         ) 

    # split into pieces
    # main variable name (e.g. CLASS or PCNT)
    lu_long['var_nam'] = lu_long.variable.str.extract(r'(\D+)')
    # identify crop number 1, 2, or 3
    lu_long['crop_num'] = lu_long.variable.str.extract(r'(\d+)')
    # extra id, relevant to irrigation with PA (irrig bool) and PB (irrig typ)
    lu_long['irr_id'] = lu_long.variable.str.extract(r'\D+\d+(\D+)')
    # add irrigation id back, to avoid NAs when pivoting wide as each crop needs both variables
    lu_long.loc[~lu_long.irr_id.isna(), 'var_nam'] += '_'+lu_long.loc[~lu_long.irr_id.isna(), 'irr_id']
    return lu_long


# %%
lu_long = crop_number_split(lu_all)

# %%
# columns to keep as indexers for crops
index_cols = ['WATERSOURC', 'MULTIUSE', 'crop_num', 
              'SURVEYAREA','SURVEYYEAR',
              'geom_id','geometry']


# %%
def make_wide_lu_df(lu_long):
    """
    The long format doesn't have individual columns for CLASS, CROPTYP, IRR_TYP_PA/B, PCNT, SPECOND, SUBCLASS
    Pivot wide for these from var_nam to make more easily referenceable
    """
    # must use pivot which requires non-duplicated entries, pivot_table aggregates
    lu_wide = lu_long.pivot(columns='var_nam', index=index_cols, values='value').reset_index()
    
    # convert subclass to numeric with coercion (make ** to NaN)
    # it is still not clear what ** are for
    lu_wide.SUBCLASS= pd.to_numeric(lu_wide.SUBCLASS, errors='coerce')
    # make -1 subclass for when one doesn't exist (avoids dealing with NAs)
    lu_wide.loc[lu_wide.SUBCLASS.isna(), 'SUBCLASS'] = -1
    return lu_wide


# %%
lu_wide = make_wide_lu_df(lu_long)


# %%
def get_full_lu_names(lu_wide, lu_class, lu_irrig):
    """
    Convert value codes into descriptions or numeric values
    INPUT:
    lu_wide is the land use data in wide format with one column for class, subclass, irr_typ, pcnt
    lu_class: dataframe to match each class, subclass with a detailed name
    lu_irrig: dataframe to match each IRR_TYP_PB with a detailed name
    OUTPUT:
    lu_wide with columns identifying the detailed name for key columns
    and strings of numbers converted to numeric 
    
    """
    # join by class and subclass to get crop name
    lu_classed = lu_wide.join(lu_class.set_index(['class', 'subclass']), on=['CLASS', 'SUBCLASS'])
    # many are still unknown but key point is knowing irrigation or not
    lu_classed = lu_classed.join(lu_irrig.set_index('code'),on='IRR_TYP_PB')
    # convert back to geodataframe
    lu_classed = gpd.GeoDataFrame(lu_classed, crs = lu_long.crs)
    # for some reason, 2015 came in as a string
    lu_classed.SURVEYYEAR = pd.to_numeric(lu_classed.SURVEYYEAR)
    
    # drop cells that were outside the survey area
    lu_classed = lu_classed[lu_classed.CLASS !='Z']
    # clean up survey area name
    lu_classed['county'] = lu_classed.SURVEYAREA.str.replace(' COUNTY','').copy()
    # simple name for plot
    lu_classed['name_plot'] = lu_classed.name.str.split(' ',expand=True)[0]
    # drop crops that don't have a name as they won't have a crop coefficient
    # these are mostly urban which can be filtered separately to identify areas of runoff, etc.
    lu_classed = lu_classed.dropna(subset='name')
    # convert PCT column to numeric and if single/double cropped the percentage should be 100
    # I, M, S have PCNT >20
    # 00 represents 100%, while ** means not used
    lu_classed.PCNT = lu_classed.PCNT.str.replace('00','100')
    lu_classed['PCNT'] = pd.to_numeric(lu_classed.PCNT, errors='coerce')
    return lu_classed


# %%
lu_classed = get_full_lu_names(lu_wide, lu_class, lu_irrig)

# %%
# subset for urban and native classes
lu_urban = lu_classed[lu_classed.name.str.contains('urban', case=False)]
lu_native = lu_classed[lu_classed.name.str.contains('native', case=False)]
# filter for those crops that are considered irrigated
lu_crops = lu_classed.dropna(subset='irr_name')
non_irrig = ['Unknown or not mapped', 'Wild Flooding']
lu_crops = lu_crops[~lu_crops.irr_name.isin(non_irrig)]

# %%
lu_plt = lu_crops[lu_crops.SURVEYYEAR==2015].copy()
# lu_plt.plot('irr_name',legend=True, legend_kwds={'ncol':1, 'loc':(1,0.2)})
# plot key categories
lu_plt['area_m2'] = lu_plt.geometry.area.copy()
lu_perc = 100*lu_plt.groupby('name').sum(numeric_only=True)['area_m2']/lu_plt.area_m2.sum()
lu_plt[lu_plt.name.isin(lu_perc[lu_perc>1].index)].plot('name',legend=True, legend_kwds={'ncol':1, 'loc':(1,0.2)})

# %%
# a 0.1% land area cutoff covers 99% of land use
cutoff = 0.2
print('%.1f %% land covered' %(lu_perc[lu_perc>cutoff].sum()), 'with %.1f %% cuttoff' % cutoff)
print('And drops %i crops' %(lu_perc.shape[0] - lu_perc[lu_perc>cutoff].shape[0]))

# %%
# S = single,  D = double and T = triple (cropped following each other), 
# I = inter cropped (orchards with annual grasses), M = multi,
# sac_lu_2000.MULTIUSE.unique()

print('%.1f%% are single cropped ' %(lu_crops[lu_crops.MULTIUSE=='S'].shape[0]*100/lu_crops.shape[0])) 

# %%
# save typical land use types to shapefile for quick reference of cultivated, urban or pasture
# from years just after new surveys
for y in [2001, 2018]:
    lu_crops['yr_diff'] = (y - lu_crops.SURVEYYEAR.copy() )
    # find the best year for land use data for each year
    pick_yrs = lu_crops[lu_crops['yr_diff']>0].groupby('county').min(numeric_only=True)
    pick_yrs['year'] = (y-pick_yrs.yr_diff)
    pick_yrs = pick_yrs[['year']].set_index('year', append=True)
    # pull out data for year
    lu_yr = lu_crops.join(pick_yrs, on=['county','SURVEYYEAR'], how='inner')
    lu_native_yr = lu_native.join(pick_yrs, on=['county','SURVEYYEAR'], how='inner')
    lu_urban_yr = lu_urban.join(pick_yrs, on=['county','SURVEYYEAR'], how='inner')
    lu_yr.to_file(join(crop_path, 'domain_ag_lu_'+str(y)+'.shp'))
    lu_native_yr.to_file(join(crop_path, 'domain_native_lu_'+str(y)+'.shp'))
    lu_urban_yr.to_file(join(crop_path, 'domain_urban_lu_'+str(y)+'.shp'))

# %%
# the land use data needs to stay in spatial format for the grid join
# saves storage by not pre-processing grid-join
lu_crops.to_file(join(crop_path, 'cleaned','lu_crops.shp'))
lu_native.to_file(join(crop_path, 'cleaned','lu_native.shp'))
