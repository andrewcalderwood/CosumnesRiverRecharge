"""
map_cln module. 
Different functions for refining plotting spatial data with general python functions
First iteration as a Module June 2022
Author: Andrew Calderwood
"""

import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.geometry import  box

import geopandas as gpd
import pandas as pd
import numpy as np

def gdf_bnds(gdf, ax=None, buf=None):
    """ Use the x-y bounds of a geodataframe to set the axes limits on a plot"""
    gdf_bnd = gpd.GeoDataFrame([0], geometry=[box(*gdf.total_bounds)], crs = gdf.crs)
    if buf != None:
        gdf_bnd.geometry = gdf_bnd.buffer(buf)
    if ax !=None:
        xmin, ymin, xmax, ymax = gdf_bnd.total_bounds
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin, ymax)
    return(gdf_bnd)

def pnt_2_tup(pnt):
    """simple function to get tuple from geopandas/shapely Point"""
    return(list(pnt[0].coords)[0])

def lab_pnt(geom, scale = None, xscale = None, yscale = None):
    """create an x,y coordinate tuple  needed for labeling a plot"""
    xy = pnt_2_tup(geom.geometry.representative_point())
    if xscale != None:
        xy = tuple([p * s for p, s in zip(xy, (xscale, yscale))])
    elif scale != None:
        xy = tuple(map((scale).__mul__, xy ))
    return(xy)

def arr_lab(gdf, text, offset = (0,0)):
    xy = gdf.geometry.unary_union.centroid.coords[0]
    ax.annotate(text=text, xy=xy, ha='center', va = 'bottom', xytext = offset, textcoords='offset pixels',
                arrowprops = {'shrinkA':1,'arrowstyle':'simple', 'color':'black'},
                bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2))


def plt_cln(ax):
    """ Set basic xy axis labels, reduce axis clutter"""
    ax.ticklabel_format(style='plain')
    plt.xlabel('Easting (m)')
    plt.ylabel('Northing (m)')
    # adjust major locators
    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    plt.yticks(rotation=90, verticalalignment = "center")

    
from shapely.geometry import box

def make_multi_scale(ax, xoff,yoff, dist = 1E3, scales = [4,2,1], units = 'km'):
    """ Plot legend with multi distances (3) at fractional axes offset
    A typical scalebar has distance and two time distances and half the distance
    With labels for just the distance and two times
    """
    minx, maxx = ax.get_xlim()
    miny, maxy = ax.get_ylim()
    # translate to actual coordinates
    lx = (maxx-minx)*xoff + minx
    ly = (maxy-miny)*yoff + miny
    height = (maxy-miny)*0.02 # height of bar
    if units=='km':
        km = 1E3
    for n, adj in enumerate(scales):
        color='black'
        if n%2==1:
            color='white'
        rect = box(lx, ly, lx +dist*adj, ly+height)
        gpd.GeoSeries(rect).plot(color=color, edgecolor='black', ax=ax)
        # round to 1 decimal, and drop decimal if it is only a .0
        dist_lab = str(np.round(dist*adj/km, 1)).replace('.0','')
        # dist_lab = str(int(dist*adj/km)) # old version
        # ax.annotate(dist_lab, (lx+dist*adj-dist*0.2,ly-2*height), xycoords='data') # original, manual adjust to make clean
        if len(dist_lab)>1:
            dec_xoff = 3
        else:
            dec_xoff = 1
        ax.annotate(dist_lab, (lx+dist*adj-dist*0.2*dec_xoff,ly-2*height), xycoords='data')
    adj = scales[0]
    ax.annotate(units, (lx+dist*adj+dist*0.4,ly-2*height), xycoords='data')