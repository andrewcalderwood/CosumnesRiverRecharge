# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
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
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

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


fig_dir = join(map_dir,'figures')

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
from map_cln import gdf_bnds, plt_cln, make_multi_scale

# %%
import map_cln
reload(map_cln)
from map_cln import gdf_bnds, pnt_2_tup, lab_pnt, plt_cln, make_multi_scale
from mf_utility import get_dates, get_layer_from_elev, clean_wb

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
prms_gis = join(gwfm_dir,'PRMS','GIS')
ws_gdf = gpd.read_file(join(prms_gis, "NHD_H_18040013_HU8_Shape\Shape\WBDHU8.shp"))
ws_gdf=ws_gdf.to_crs(m_domain.crs)

# %%
teichert = gpd.read_file(gwfm_dir+'/Mapping/Kautz_shapefiles/Kautz Property.shp')
rooney = gpd.read_file(gwfm_dir+'/Mapping/Kautz_shapefiles/Kautz Property2.shp')
mosher = gpd.read_file(gwfm_dir+'/Mapping/Kautz_shapefiles/MosherProperty.shp')
vineyard = pd.concat((teichert, rooney, mosher))

vineyard=vineyard.to_crs('epsg:32610')
teichert = teichert.to_crs('epsg:32610')

# %%
# scenario files
proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')

floodplain_clean = gpd.read_file(join(proj_dir, 'scenarios', 'R5_floodplain_MAR_clean_outline.shp'))
floodplain = gpd.read_file(join(proj_dir, 'scenarios', 'R4_floodplain_MAR_outline.shp'))


# %%
# def regional_scale_arrow(ax):
#     x, y, arrow_length = 0.7, 0.15, 0.1
#     ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
#                 arrowprops=dict(facecolor='black', width=5, headwidth=15),
#                 ha='center', va='center', fontsize=20,
#                 xycoords=ax.transAxes)

#     ax.ticklabel_format(style='plain')


#     fontprops = fm.FontProperties(size=18)
#     scalebar = AnchoredSizeBar(ax.transData,
#                                10000, '10 km', 'lower right', 
#                                pad=0.1,sep=2, color='black',
#                                frameon=False,size_vertical=5E2, fontproperties=fontprops)
#     ax.add_artist(scalebar)

# %%

def ref_map(axins, fontsize=10, lw=1):
    ca.plot(ax = axins,alpha = 0.2, edgecolor='black')
    axins.annotate(text='California', xy=lab_pnt(ca), 
                xytext = (6,-30), textcoords = 'offset points', fontsize=fontsize, 
                bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=lw), zorder=1)
    # second one is smaller inset
    axins.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

    ws_gdf.plot(color='none',ax=axins, linestyle='--')
    # ws_gdf.plot(color='none',ax=axins, linestyle='--')
    m_domain.plot(color="none",edgecolor='black',ax=axins)
    return axins



# %%
def regional_arrow(ax, xoff, yoff):
    x, y, arrow_length = xoff, yoff, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20,
                xycoords=ax.transAxes)
    

def main_map(ax):

    m_domain.plot(ax=ax,color="none",edgecolor='black')
    gdf_bnds(m_domain,ax)
    
    rivers.loc[rivers.GNIS_Name.isin(['Mokelumne River','South Mokelumne River'])].plot(ax=ax)
    cr = rivers.loc[rivers.GNIS_Name=='Cosumnes River']
    cr.plot( ax=ax,label='Cosumnes River')
    
    ax.annotate(text='Cosumnes\n River', xy=list(cr.geometry.iloc[10].centroid.coords)[0], 
                xytext = (6,6), textcoords = 'offset pixels',
                bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2))
    # plot floodplain first since larger shape
    floodplain_clean.plot(ax=ax, color='none', alpha = 0.7, edgecolor='blue', hatch='||')
    ax.annotate(text='Floodplain\nMAR', xy=list(floodplain_clean.iloc[0].geometry.centroid.coords)[0], 
                xytext = (-100,-200), textcoords = 'offset pixels',
                bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2))
    
    # plot vineyards on top
    vineyard.plot(ax=ax, color='green', alpha = 0.7,edgecolor='black')
    ax.annotate(text='Vineyard\nMAR', xy=list(vineyard.iloc[0].geometry.centroid.coords)[0], 
                xytext = (60,-20), textcoords = 'offset pixels',
                bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2))


    ctx.add_basemap(ax, source = ctx.providers.Esri.WorldImagery, crs='epsg:26910', alpha = 0.8, attribution=False)



# %%
def arr_lab(gdf, text, ax, offset = (0,0), arrow=False, exterior = False, fontsize=10):
    xy = gdf.geometry.unary_union.centroid.coords[0]
    lw = 1
    if exterior:
        xy = gdf.geometry.unary_union.exterior.representative_point().centroid.coords[0]
    if arrow:
        ax.annotate(text=text, xy=xy, ha='center', va = 'bottom', xytext = offset, textcoords='offset pixels', fontsize = fontsize, 
                    arrowprops = {'shrinkA':1,'arrowstyle':'simple', 'color':'black'},
                    bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=lw))
    else:
        ax.annotate(text=text, xy=xy, ha='center', va = 'bottom', xytext = offset, textcoords='offset pixels', fontsize = fontsize, 
            bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=lw))

def xy_lab(xy, text, offset = (0,0), lw=1, fontsize=10, bbox=True, fc='white', ec='black'):
    if bbox:
        ax.annotate(text=text, xy=xy, ha='center', va = 'bottom', xytext = offset, textcoords='offset pixels', fontsize = fontsize, 
                    bbox=dict(boxstyle="square,pad=0.3", fc=fc, ec=ec, lw=lw))
    else:
        ax.annotate(text=text, xy=xy, ha='center', va = 'bottom', xytext = offset, textcoords='offset pixels', fontsize = fontsize)


# %%
lgd_short = [
    Patch(facecolor='red', edgecolor='red',alpha=0.5, label='Levee removal and\nexcavation areas'),
    Line2D([0], [0], marker='.', linestyle='', color='black', label='Monitoring Well'),
]

# %%
fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=300)

main_map(ax=ax)
fontsize=10


regional_arrow(ax, 0.65, 0.15)
make_multi_scale(ax, 0.75,0.1, dist=2E3)
# regional_arrow(ax, 0.05, 0.15)
# make_multi_scale(ax, 0.1, 0.1, dist = 2E3, scales = [4,2,1])
plt_cln(ax=ax)

## inset map
axins = inset_axes(ax, width="35%", height="35%", 
                   loc='upper left',
                   bbox_to_anchor=(0, 0, 1, 1), 
                   # bbox_to_anchor=(-.01, .01, 1, 1),
                  bbox_transform=ax.transAxes, 
                   # loc=2
                  )

ref_map(axins, fontsize=8)
arr_lab(m_domain, 'Cosumnes\nWatershed', axins, offset = (125, 40), arrow=False, fontsize=8)

# plt.savefig(join(fig_dir, 'regional_domain_map.png'),  bbox_inches='tight')

