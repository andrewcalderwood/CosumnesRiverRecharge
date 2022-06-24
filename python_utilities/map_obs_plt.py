"""
map_obs_plt module. 
Different functions for plotting map and output head data from MODFLOW in python based on flopy
First iteration as a Module June 2022
Author: Andrew Calderwood
"""
import numpy as np
import geopandas as gpd

import matplotlib.pyplot as plt
import flopy

def plt_bc_hk(model, ax, hk=False):
    mapview = flopy.plot.PlotMapView(model=model,ax=ax)

    # plot the horizontal hydraulic conductivities
    if hk == True:
        a = model.lpf.hk.array
        csa = mapview.plot_array(a, norm=mpl.colors.LogNorm())
        cb = plt.colorbar(csa, shrink=0.75,ax=ax)
        cb.set_label('Horiz. Cond. (m/d)')
    # Plots boundary condtiions
#     quadmesh = mapview.plot_ibound(ax=ax)
    mapview.plot_bc("GHB", plotAll=True,ax=ax)
    mapview.plot_bc("SFR", plotAll=True,ax=ax)
    mapview.plot_bc("CHD", plotAll=True,ax=ax)
    mapview.plot_bc("WEL", plotAll=True,ax=ax,alpha=0.1, color='red')

    ax.ticklabel_format(style='plain')
#     linecollection = mapview.plot_grid(linewidths = 0.3,ax=ax)
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')

def plot_head_simple(model, hdobj, cbb, spd_stp,time, name, ax,units):
    hd_step = 5
    if units =='ft':
        fact = 0.3048
    elif units== 'm':
        fact=1
    head = hdobj.get_data(kstpkper = spd_stp)/fact
    levels_hmax = np.round(np.max(head[head!=1E+30/fact]),-1)
    levels_hmin = np.round(np.min(head[head>=-200/fact]),-1)
    levels = np.arange(levels_hmin, levels_hmax, int(hd_step/fact))

    ax.set_title(name+' Head Contours at '+ str(time))
    mapview = flopy.plot.PlotMapView(model=model,ax=ax)

    contour_set = mapview.contour_array(head, masked_values=[-999.99], levels=levels, ax=ax)
    hcb = plt.colorbar(contour_set, shrink = 0.5,ax=ax)
    hcb.set_label('Head ('+units+')')
    ax.clabel(contour_set, contour_set.levels[0::2], inline=True, fontsize=8)

#     quiver = mapview.plot_discharge(frf, fff, istep=10, jstep=10)  # no head array for volumetric discharge
    ax.ticklabel_format(style='plain')
#     linecollection = mapview.plot_grid(linewidths = 0.3,ax=ax)
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
#     fig.tight_layout()

def plot_dtw_simple(model, hdobj, cbb, spd_stp,time, name, ax,units):
    hd_step = 10
    if units =='ft':
        fact = 0.3048
    elif units== 'm':
        fact=1
    head = hdobj.get_data(kstpkper = spd_stp)/fact
    dtw = model.dis.top.array-head
    levels_dtw_max = np.round(np.max(dtw[head!=1E+30/fact]),-1)
    # filter out dtw greater than 300 because of mountain ranges
#     levels_dtw_max = 300
    levels_dtw_hmin = np.round(np.min(dtw[head>=-200/fact]),-1)
#     levels_dtw_hmin = -300
        #     frf = cbb.get_data(text='FLOW RIGHT FACE')[0]
    #     fff = cbb.get_data(text='FLOW FRONT FACE')[0]

    levels = np.arange(levels_dtw_hmin, levels_dtw_max, int(hd_step/fact))

    ax.set_title(name+' DTW Contours at '+ str(time))
    mapview = flopy.plot.PlotMapView(model=model,ax=ax)
    
    contour_set = mapview.contour_array(dtw, masked_values=[-999.99], levels=levels, ax=ax)
    hcb = plt.colorbar(contour_set, shrink = 0.5,ax=ax)
    hcb.set_label('DTW ('+units+')')
    ax.clabel(contour_set, contour_set.levels[0::2], inline=True, fontsize=8)
#     ax.clabel(contour_set, contour_set.levels[0::2], inline=True, fontsize=8)

#     quiver = mapview.plot_discharge(frf, fff, istep=10, jstep=10)  # no head array for volumetric discharge
    ax.ticklabel_format(style='plain')
#     linecollection = mapview.plot_grid(linewidths = 0.3,ax=ax)
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
