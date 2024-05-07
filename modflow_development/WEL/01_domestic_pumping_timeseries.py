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
# rural_grid = pd.read_csv(join(wel_dir, 'ag_res_parcel_domestic_wells.csv'))

# %% [markdown]
# # Calculate water use time-series
# Assume baseline water usage for an ag-residential parcel of 1 acre is 2AF per year with a minimum daily usage of 50 gpd per capita with 3 people per parcel and the remaining daily water usage scaled by the ETc.

# %%
import h5py
nrow_p, ncol_p = (100,230)
def dwr_etc(strt_date, end_date):
    nper_tr = (end_date-strt_date).days+1
    natETc = np.zeros((nper_tr,nrow_p,ncol_p))
    agETc = np.zeros((nper_tr,nrow_p,ncol_p))

    per_n = 0 
    for y in np.arange(strt_date.year, end_date.year+1):
        # set start and end date for range for the year to be iterated over
        yr_strt = pd.to_datetime(str(y)+'-01-01')
        yr_end = pd.to_datetime(str(y)+'-12-31')
        # get the length of the date range needed for that year
        yearlen = len(pd.date_range(yr_strt, yr_end))
        if yr_strt < strt_date:
            yr_strt = strt_date
        if yr_end > end_date:
            yr_end = end_date
        yr_len = len(pd.date_range(yr_strt, yr_end))
        # load hdf5 files
        f_irr = h5py.File(join(uzf_dir, "dwr_ETc/irrigated_"+str(y)+".hdf5"), "r")
        agETc[per_n:per_n+yr_len,:,:] = f_irr['array'][str(y)][:][yr_strt.dayofyear-1:yr_end.dayofyear,:,:]
        f_irr.close()
        f_nat = h5py.File(join(uzf_dir, "dwr_ETc/native_"+str(y)+".hdf5"), "r")
        natETc[per_n:per_n+yr_len,:,:] = f_nat['array'][str(y)][:][yr_strt.dayofyear-1:yr_end.dayofyear,:,:]
        f_nat.close()
        per_n += yr_len
    # make sure the return value is separate from the loop
    return(agETc, natETc)


# %%
# strt_date = pd.to_datetime('1998-10-01')
# end_date = pd.to_datetime('2020-9-30')

# %%
uzf_dir = join(gwfm_dir, 'UZF_data')

# %%
# available years of data
years = pd.Series(os.listdir(join(uzf_dir, 'dwr_ETc'))).str.extract(r'(\d{4})').dropna().astype(int).values
strt_date = pd.to_datetime(str(years.min())+'-1-1')
end_date = pd.to_datetime(str(years.max())+'-12-31')

# %%
agETc, natETc = dwr_etc(strt_date, end_date)
# net ETc should be ETc from ag and native plants joined
ETc = agETc + natETc

# %%
years = pd.date_range(strt_date, end_date+pd.DateOffset(days=1), freq='AS-Oct')
year_ind = (years-strt_date).days

# %%
# take average daily ET across the domain to use for scaling the domestic use
mean_ETc = agETc.mean(axis=(1,2))
scaled_ETc = mean_ETc.copy()
# group ET on an annual scale to determine the contribution each day makes to the year
for n in np.arange(0, len(years)-1):
    yr_scaled = mean_ETc[year_ind[n]:year_ind[n+1]]/mean_ETc[year_ind[n]:year_ind[n+1]].sum()
    scaled_ETc[year_ind[n]:year_ind[n+1]] = yr_scaled

# %%
# the average minimum daily pumping in a single family household is 50 gpd for ~3 people (0.16 AF/yr)
min_pump = 150*(1/7.48)*(0.3048**3)
# typical max usage is 2AF/year for domestic
annual_domestic_irr = 2*43560*(0.3048**3) - min_pump*365
# the total annual pumping is downscaled to daily based on ET and minimum pumping is added back
domestic_use = scaled_ETc*annual_domestic_irr+min_pump


# %%
plt.plot(domestic_use)

# %%
# save to output with dates as ID
dom_use_out = pd.DataFrame(domestic_use, columns=['flux_m3d'])
dom_use_out.index = pd.date_range(strt_date, end_date)
dom_use_out.to_csv(join(wel_dir, 'domestic_water_use.csv'))

# %%
# for n in np.arange(0, len(year_ind)-1):
#     print(domestic_use[year_ind[n]:year_ind[n+1]].sum()*(1/0.3048**3)/43560)

# %%
