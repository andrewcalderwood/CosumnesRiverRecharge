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

# standard geospatial python utilities
import geopandas as gpd
import rasterio


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
def make_wide_lu_df(lu_long, index_cols):
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
def get_full_lu_names(lu_wide, lu_class, lu_irrig, crs):
    """
    Convert value codes into descriptions or numeric values
    INPUT:
    lu_wide is the land use data in wide format with one column for class, subclass, irr_typ, pcnt
    lu_class: dataframe to match each class, subclass with a detailed name
    lu_irrig: dataframe to match each IRR_TYP_PB with a detailed name
    crs: crs of geometry in lu_wide
    OUTPUT:
    lu_wide with columns identifying the detailed name for key columns
    and strings of numbers converted to numeric 
    
    """
    # join by class and subclass to get crop name
    lu_classed = lu_wide.join(lu_class.set_index(['class', 'subclass']), on=['CLASS', 'SUBCLASS'])
    # many are still unknown but key point is knowing irrigation or not
    lu_classed = lu_classed.join(lu_irrig.set_index('code'),on='IRR_TYP_PB')
    # convert back to geodataframe
    lu_classed = gpd.GeoDataFrame(lu_classed, crs = crs)
    if 'SURVEYYEAR' in lu_classed.columns:
        # for some reason, 2015 came in as a string
        lu_classed.SURVEYYEAR = pd.to_numeric(lu_classed.SURVEYYEAR)
    
    # drop cells that were outside the survey area
    lu_classed = lu_classed[lu_classed.CLASS !='Z']
    if 'SURVEYAREA' in lu_classed.columns:
        # clean up survey area name
        lu_classed['county'] = lu_classed.SURVEYAREA.str.replace(' COUNTY','').copy()
    elif 'COUNTY' in lu_classed.columns:
        lu_classed = lu_classed.rename(columns={'COUNTY':'county'})
    # simple name for plot
    lu_classed['name_plot'] = lu_classed.name.str.split(' ',expand=True)[0]
    # drop crops that don't have a name as they won't have a crop coefficient
    # these are mostly urban which can be filtered separately to identify areas of runoff, etc.
    # lu_classed = lu_classed.dropna(subset='name')
    # convert PCT column to numeric and if single/double cropped the percentage should be 100
    # I, M, S have PCNT >20
    # 00 represents 100%, while ** means not used
    lu_classed.PCNT = lu_classed.PCNT.str.replace('00','100')
    lu_classed['PCNT'] = pd.to_numeric(lu_classed.PCNT, errors='coerce')
    return lu_classed