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
# standard python utilities
import os
from os.path import dirname, basename, exists, join
import sys
import glob
import pandas as pd
import numpy as np
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import shapely
import geopandas as gpd
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

from importlib import reload

# %%
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
    
git_dir = join(doc_dir,'GitHub','CosumnesRiverRecharge')
# dir of all gwfm data
gwfm_dir = join(dirname(doc_dir),'Box/research_cosumnes/GWFlowModel')




# %%
def add_path(fxn_dir):
    if fxn_dir not in sys.path:
        sys.path.append(fxn_dir)
        
add_path(doc_dir+'/GitHub/CosumnesRiverRecharge/python_utilities')

# %%
from map_cln import gdf_bnds, pnt_2_tup, lab_pnt, plt_cln


# %%
# set box directory for output figures and data
box_dir = gwfm_dir+'/Levee_setback/levee_setback_distance_analysis/'

# tprogs_id = '' # original tprogs with conditioning data in output tsim
# tprogs_id = '_no_conditioning'
tprogs_id = '_no_cond_c3d'

data_dir = box_dir+ tprogs_id+'/data_output/'
# fig_dir = box_dir+tprogs_id+'/figures/'
fig_dir = box_dir+'/figures/'

chan_dir = box_dir+'channel_data/'
gis_dir = chan_dir+'GIS/'

# %%
grid_sfr = gpd.read_file(gwfm_dir+'/SFR_data/final_grid_sfr/grid_sfr.shp')


# %%

import functions.setback_plotting
reload(functions.setback_plotting)
from functions.setback_plotting import sfr_setback, floodplain_annotate, arr_plt_label_stream

# %%
# vmax for depth/recharge are manually specified, #'viridis_r', 
#'YlOrRd': hard to see yellow, # plasma was nice, magma/inferno are intense with black


# %%

def ex_cbar(fig,ax, im_K, im_d, im_rch):
    cbar0 = ax[0].ravel()
    cbar1 = ax[1].ravel()
    cbar2 = ax[2].ravel()
    shrink=0.6
    plt.colorbar(im_K, ax=cbar0, orientation='vertical', label='Vertical Conductivity\n(m/d)', 
                 shrink=shrink,location='right')
    plt.colorbar(im_d, ax=cbar1, orientation='vertical', label='Depth (m)', shrink=shrink,location='right')
    plt.colorbar(im_rch, ax=cbar2, orientation='vertical', label='Total Recharge (MCM)', 
                         shrink=shrink,location='right')


# %%
s_list = [3,6]


# %%

# save files for Helen to use
soil_K = np.loadtxt(join(fig_dir, 'to_Helen','Figure05','vertical_conductivity_m_d.txt'))

s_list = [3,6]
d_out = np.empty(np.append(len(s_list), (soil_K).shape))
rch_out = np.empty(np.append(len(s_list), (soil_K).shape))
for n,s in enumerate(s_list):
    d_out[n] = np.loadtxt(join(fig_dir, 'to_Helen','Figure05', 'flood_depth_m_array_'+str(s*200)+'.txt'))
    rch_out[n] = np.loadtxt(join(fig_dir, 'to_Helen','Figure05', 'total_recharge_mcm_array_'+str(s*200)+'.txt'))

setback_outer = np.empty(np.append(len(s_list), (soil_K).shape))
for n, s in enumerate(s_list):
    setback_outer[n] = np.loadtxt(join(fig_dir,'to_Helen', 'setback_array_'+str(s*200)+'m.txt'))

setback_outer = np.ma.masked_where(setback_outer==0,setback_outer)
sfr_array = np.loadtxt(join(fig_dir,'to_Helen', 'cosumnes_streamline_array.txt'))
# %%
wet_frac = pd.read_csv(join(fig_dir,'to_Helen', 'Fig07_inundated_fraction_by_setback.csv'))
def plt_wet_frac(wet_frac, ax):
    setbacks = wet_frac.setback_dist_m.drop_duplicates().sort_values().values
    colors = ['darkviolet', 'darkorange', 'lime']
    
    for nf, ft in enumerate(ft_in):
        # plot flooded area divided by setback area with setback distance for x-axis
        wet_frac[wet_frac.ft==ft].plot(x='setback_dist_m', y='inundated_fraction', ax=ax,
                                       label='Type '+str(ft),color=colors[nf])
    # set x labels for boxplots 
    ax.set_xticks(setbacks[::2], setbacks[::2], rotation=0, ha="right", rotation_mode="anchor")
    ax.set_xticks(setbacks[1::2],minor=True)
    ymin, ymax = wet_frac.inundated_fraction.quantile([0,1]).round(1)
    ax.set_yticks(np.linspace(ymin, ymax,4))
    ax.set_xlabel('Setback Distance (m)')#,y=-0.04)
    ax.set_ylabel('Inundated Fraction')#, x=-0.01)
    ax.legend()
    # plt.savefig(join(fig_dir, 'inundated_fraction_by_setback.tif'), bbox_inches='tight')
# fig,ax = plt.subplots()
# line_ax = inset_axes(ax, width="100%", height="100%", loc='lower left',
#             bbox_to_anchor=(0.2,0.2, .7,.7), bbox_transform=ax.transAxes)
# plt_wet_frac(wet_frac, line_ax)



# %%
ft_in = [1,2,3]
d_ft_6 = np.empty(np.append(len(ft_in), (soil_K).shape))
d_ft_12 = np.empty(np.append(len(ft_in), (soil_K).shape))
for n,ft in enumerate(ft_in):
    s=3
    d_ft_6[n] = np.loadtxt(join(fig_dir, 'to_Helen','Figure06', 'flood_depth_m_array_Type'+str(ft)+'_'+str(s*200)+'m.txt'))
    s=6
    d_ft_12[n] = np.loadtxt(join(fig_dir, 'to_Helen','Figure06', 'flood_depth_m_array_Type'+str(ft)+'_'+str(s*200)+'m.txt'))
    

# %%
def cbaxes_label(ax_n, label, im_arr, frame=False):
    # cbaxes = inset_axes(ax_n, width="5%", height="20%", loc='upper center',
    #                 bbox_to_anchor=(0.3,-0.1, 1,1), bbox_transform=ax_n.transAxes) #, loc=('lower left')
    cbaxes = inset_axes(ax_n, width="5%", height="20%", loc='upper center',
                    bbox_to_anchor=(0.35,-0.2, 1,1), bbox_transform=ax_n.transAxes) #
    plt.colorbar(im_arr, cax=cbaxes,  #label='Recharge\n(MCM)', 
                         shrink=shrink, orientation='vertical',
                )
    # find bbox of the ytick labels for the colorbar
    bbox = cbaxes.get_yticklabels()[-1].get_window_extent()
    labx,_ = cbaxes.transAxes.inverted().transform([bbox.x1,0])
    ax_n.annotate(label, xy=(0,1.2), xycoords=cbaxes.transAxes,
                 ha = "center", va="bottom",
                  bbox={'facecolor': 'white', 'alpha': 0.9, 'pad': 1, 'edgecolor':'white'}
                  )
    if frame:
        bbox = cbaxes.get_position()
        rect = plt.Rectangle((bbox.x0, bbox.y0), bbox.width*1.5, bbox.height*1.5, facecolor='white', edgecolor='black',alpha=0.9, zorder=-1)
        ax_n.add_artist(rect)
    return(cbaxes)



# %%
import matplotlib.patches as patches

# %%
fig,ax = plt.subplots(ny, nx, figsize=(2.3*nx*ax_l, ny*ax_l), dpi=300, layout='constrained')
ax_n = ax[0,0]
im_K = ax_n.imshow(soil_K, norm= mpl.colors.LogNorm(),cmap='viridis', alpha=0.7)

cbaxes = cbaxes_label(ax_n, 'Vertical hydraulic\nconductivity (m/d)', im_K, frame=True)
# bbox = cbaxes.get_position()
# rect = patches.Rectangle((20, 180), 30, 50, facecolor='white', edgecolor='black',alpha=0.9, zorder=4)
# ax[0,0].add_patch(rect)

# %%
np.round(np.nanmax(d_ft_6),1),np.round(np.nanmax(d_ft_12),1), np.nanmax(d_out)

# %% [markdown]
# After reviewing, I realized that the same depth data was plotted twice essentially for one flow type. There might have been a slight difference in the averaging vs peak but the mapped results are near identical. These two subplots can be removed once leaving 10 instead of 12. Then the conductivity subplot and inundation line plot can be at the top with the rest below.

# %%
shrink=0.6
ax_l=2
nx,ny = (2,5)
# the plots should have a width 2.3 times width based on number of cells
fig,ax = plt.subplots(ny, nx, figsize=(2.3*nx*ax_l, ny*ax_l), dpi=300, layout='constrained')
print('figsize:',2.3*2*l, 6*l)
# ax.set_axis_off()
plt.setp(ax, xticks=[], yticks=[])
## map hydraulic conductivity
ax_n = ax[0,0]
im_K = ax_n.imshow(soil_K, norm= mpl.colors.LogNorm(),cmap='viridis', alpha=0.7)

cbaxes_label(ax_n, 'Vertical hydraulic\nconductivity (m/d)', im_K, frame=True)
ax_n.plot(grid_sfr.column-1, grid_sfr.row-1, color='blue', linewidth=0.5, linestyle='-')

## line plot of inundated fraction
# need to create an inset axis within the main axis
ax_n = ax[0,1]
line_ax = inset_axes(ax_n, width="100%", height="100%", loc='lower left',
            bbox_to_anchor=(0.15, 0.2, .8,.7), bbox_transform=ax_n.transAxes)
plt_wet_frac(wet_frac, line_ax)
ax_n.axis('off')

## map of depth then recharge for two setbacks
for n in [0,1]:
    ax_n = ax[3+n,0]
    # recharge 
    temp = rch_out[n]
    im_rch = ax_n.imshow(np.ma.masked_where(temp==0, temp), vmax=np.round(rch_out.max(),1), cmap = 'viridis' ) 
    cbaxes_label(ax_n, 'Recharge (MCM)', im_rch)
    # depth
    ax_n = ax[1+n,0]
    im_d = ax_n.imshow(np.ma.masked_invalid(d_out[n]), vmax=np.round(np.nanmax(d_out),1), cmap = 'plasma' ) 
    cbaxes_label(ax_n, 'Depth (m)', im_d)
## map of depth for two flood types and two setbacks
for n, ft in enumerate([0,2]):
    ax_n = ax[1+n*2,1]
    im_d = ax_n.imshow(np.ma.masked_invalid(d_ft_6[ft]), vmax=np.round(np.nanmax(d_ft_6),1), cmap = 'plasma' ) 
    cbaxes_label(ax_n, 'Depth (m)', im_d)
for n, ft in enumerate([0,2]):
    ax_n = ax[2+n*2,1]
    im_d = ax_n.imshow(np.ma.masked_invalid(d_ft_12[ft]), vmax=np.round(np.nanmax(d_ft_12),1), cmap = 'plasma' ) 
    cbaxes_label(ax_n, 'Depth (m)', im_d)  
    
## add annotation of floodplain to plots as neeeded
ann_600 = [(1,0), (3,0), (1,1),(3,1)]
for n, l in enumerate(ann_600):
    ax_n = ax[l]
    ax_n.imshow(setback_outer[0], cmap='gray') # works when dpi=600
    # ax_n.plot(grid_sfr.column-1, grid_sfr.row-1, color='black', linewidth=0.5, linestyle='-.')
ann_600 = [(0,0), (2,0), (4,0),(2,1),(4,1)]
for n, l in enumerate(ann_600):
    ax_n = ax[l]
    ax_n.imshow(setback_outer[1], cmap='gray') 
title = ['A','B','C','D','E']
title2 = ['F','G','H','I','J']

for n in np.arange(0, len(title)):
    ax_n = ax[n,0]
    ax_n.annotate(text = title[n], xy=(0.025,0.9), xycoords='axes fraction',
                  # bbox={'facecolor': 'lightgray', 'alpha': 0.9, 'pad': 2}
                 )
    ax_n = ax[n,1]
    ax_n.annotate(text = title2[n], xy=(0.025,0.9), xycoords='axes fraction',
                  # bbox={'facecolor': 'lightgray', 'alpha': 0.9, 'pad': 2}
                 )

arr_plt_label_stream(ax[1,1], xy_down=(10,30), xy_up=(140,90), xy_arrow = (10,90))
# fig.tight_layout()
plt.savefig(join(fig_dir, 'example_vka_depth_recharge_mean_inun.tif'), bbox_inches='tight')

# %%
# # took manual adjustments to find 11, 6.5 with h_pad-1 to work
# fig,ax = plt.subplots(3,2, dpi=300,figsize=(10, 5), sharex=True,
#                      )
# im_K, im_d, im_rch = ex_depth_rch(t, 3, fig,ax[:,0], cbar_log=False)
# setback_outer = sfr_setback(grid_sfr, grid_p, 3)
# floodplain_annotate(ax[:,0], setback_outer, grid_sfr, title=['A', 'C','E'])
    
# im_K, im_d, im_rch = ex_depth_rch(t, 6, fig,ax[:,1], cbar_log=False)
# setback_outer = sfr_setback(grid_sfr, grid_p, 6)
# floodplain_annotate(ax[:,1], setback_outer, grid_sfr, title=['B', 'D','F'])

# arr_plt_label_stream(ax=ax[2,0])
# arr_plt_label_stream(ax=ax[2,1])

# # fig.tight_layout(h_pad=-1.5, w_pad=1)
# fig.tight_layout(h_pad=.25) #pad=0.4, w_pad=0.5,
# # fig.tight_layout(h_pad=0, w_pad=0)
# ex_cbar(fig,ax, im_K, im_d, im_rch )

# # fig.tight_layout(h_pad=-.3) #pad=0.4, w_pad=0.5,

# # plt.savefig(join(fig_dir, 'example_vka_depth_recharge_side_by_side_alt.tif'), bbox_inches='tight')
