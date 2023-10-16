# ---
# jupyter:
#   jupytext:
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

# %%
import os
from os.path import basename, dirname, exists, join
import sys

import geopandas as gpd
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %%
## Set up directory referencing
# Package data
usr_dir = os.getcwd()
while basename(usr_dir)!='Users':
    temp = basename(usr_dir)
    usr_dir = dirname(usr_dir)
usr_dir = join(usr_dir, temp)
gwfm_dir = join(usr_dir,'Box/research_cosumnes/GWFlowModel')
doc_dir = join(usr_dir, 'Documents')

map_dir = join(gwfm_dir,'Mapping')
sfr_dir = join(gwfm_dir,'SFR_data')

plt_dir = join(map_dir,'figures')

# %%
from importlib import reload


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

from mf_utility import get_layer_from_elev
from map_cln import gdf_bnds, plt_cln

# %%
# reload(map_cln)

# %%
grid_sfr = gpd.read_file(join(sfr_dir, 'final_grid_sfr/grid_sfr.shp'))


# %%
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')


# %%
m_domain = gpd.read_file(gwfm_dir+"/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp")

rivers = gpd.read_file(join(sfr_dir,"NHD_MajorRivers/NHD_MajorRivers.shp")).to_crs(m_domain.crs)
rivers_clip = gpd.clip(rivers, m_domain)

soam = gpd.read_file(map_dir+"/so_am_subbasin/so_am_subbasin.shp").to_crs(m_domain.crs)
cos = gpd.read_file(map_dir+"/cos_subbasin/cos_subbasin.shp").to_crs(m_domain.crs)
ca = gpd.read_file(map_dir+"/ca_state_boundary/CA_State_TIGER2016.shp").to_crs(m_domain.crs)


# %%
def regional_scale_arrow(ax):
    x, y, arrow_length = 0.7, 0.15, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20,
                xycoords=ax.transAxes)

    ax.ticklabel_format(style='plain')


    fontprops = fm.FontProperties(size=18)
    scalebar = AnchoredSizeBar(ax.transData,
                               10000, '10 km', 'lower right', 
                               pad=0.1,sep=2, color='black',
                               frameon=False,size_vertical=5E2, fontproperties=fontprops)
    ax.add_artist(scalebar)


# %%
fig,ax = plt.subplots(figsize=(6.5,6.5), dpi=300)

m_domain.plot(ax=ax,color="none",edgecolor='black')
gdf_bnds(m_domain,ax)

rivers.loc[rivers.GNIS_Name.isin(['Mokelumne River','South Mokelumne River'])].plot(ax=ax)
cr = rivers.loc[rivers.GNIS_Name=='Cosumnes River']
cr.plot( ax=ax,label='Cosumnes River')

ax.annotate(text='Cosumnes\n River', xy=list(cr.geometry.iloc[10].centroid.coords)[0], 
            xytext = (6,6), textcoords = 'offset pixels',
            bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2))


ctx.add_basemap(ax, source = ctx.providers.Esri.WorldImagery, crs='epsg:26910', alpha = 0.8, attribution=False)

regional_scale_arrow(ax=ax)
plt_cln(ax=ax)

# first one is CA map
axins = inset_axes(ax, width="80%", height="80%", bbox_to_anchor=(0, .5, .5, .5),
                  bbox_transform=ax.transAxes, loc=2)
# second one is smaller inset
axins2 = inset_axes(ax, width="40%", height="40%", bbox_to_anchor=(0.18, .5, .5, .5),
                  bbox_transform=ax.transAxes, loc=2)
axins.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)
axins2.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

ca.plot(ax = axins,alpha = 0.2)
soam.plot(ax = axins, edgecolor = 'black')
cos.plot(ax = axins, edgecolor = 'black')

cos.plot(ax = axins2, edgecolor = 'black')
soam.plot(ax = axins2, edgecolor = 'black')
m_domain.plot(ax = axins2, edgecolor = 'black', color = 'none')


plt.savefig(join(plt_dir, 'regional_domain_map.png'),  bbox_inches='tight')



# %% [markdown]
# # Plot Major Land Uses (GW Pumping)

# %%
uzf_dir = join(gwfm_dir, 'UZF_data')

lu_ag = gpd.read_file(join(uzf_dir, 'county_landuse', 'domain_ag_lu_2018.shp'))
lu_ag = gpd.overlay(lu_ag, m_domain)
lu_ag['area_m2'] = lu_ag.geometry.area
wel_dir = join(gwfm_dir, 'WEL_data')

# load data of locations of domestic wells
dom_loc = pd.read_csv(join(wel_dir, 'ag_res_parcel_domestic_wells.csv'), index_col=0)
# make row,column 0 based
dom_loc.row = (dom_loc.row-1).astype(int)
dom_loc.column = (dom_loc.column -1).astype(int)
# aggregate to the cell level, summing area will keep water usage scaling correct
dom_loc = dom_loc.groupby(['node','row','column', 'CITY']).sum(numeric_only=True).reset_index()
# join to grid for plotting
dom_grid = grid_p.merge(dom_loc,on=['row','column'])

# %%
# not super helpful since dometsic well parcel mapping will cover most
# lu_urban = gpd.read_file(join(uzf_dir, 'county_landuse', 'domain_urban_lu_2018.shp'))


# %%
lu_crop_sum = lu_ag.groupby('name').sum(numeric_only=True).reset_index()

min_area = 1000*43560*0.3048**2 # minimum area to include in plot is 1000 acres
min_area = 200*200*100*230*0.01 # minimum area is 1% of the domain
print('Min area %.2e m^2' %min_area,'vs number of cells %.i' %(min_area/(200**2)), 
      'vs percent of domain %.2f' %(min_area*100/(200**2)/(100*230) ))
main_crops = lu_crop_sum[lu_crop_sum.area_m2 > min_area].name.unique()

# %%
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import cm


# %%
fig,ax = plt.subplots(figsize=(6.5,6.5), dpi=300)

#normalize item number values to colormap
norm = mpl.colors.Normalize(vmin=0, vmax=len(main_crops))

m_domain.plot(ax=ax,color="none",edgecolor='black')
gdf_bnds(m_domain,ax)


custom_lgd = []
for nc, c in enumerate(main_crops):
    lu_ag[lu_ag.name.isin([c])].plot(ax=ax, color=cm.viridis(norm(nc)), label=c)
    custom_lgd += [Patch(facecolor=cm.viridis(norm(nc)), edgecolor='None', label=c)]

# rural residences
nc +=1
dom_grid.plot(ax=ax, label='Rural residences', color=cm.viridis(norm(nc)))
custom_lgd += [Patch(facecolor=cm.viridis(norm(nc)), edgecolor='None', label='Rural residences')]
# stream grid cells
grid_sfr.plot(ax=ax)
custom_lgd += [Line2D([0], [0], color='tab:blue', lw=4, label='Stream Cells')]

plt_cln(ax=ax)

ax.legend(handles=custom_lgd, facecolor='lightgray', framealpha=0.6, loc='upper left')

plt.savefig(join(plt_dir, 'regional_major_land_use_map.png'),  bbox_inches='tight')


# %% [markdown]
# # Plot of geology (XS or 3D view)
# I have a 2D cross-section with starting heads from the model in the Obs_plotting script.

# %%
