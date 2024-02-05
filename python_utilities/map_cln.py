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

# def arr_lab(gdf, text, offset = (0,0)):
#     xy = gdf.geometry.unary_union.centroid.coords[0]
#     ax.annotate(text=text, xy=xy, ha='center', va = 'bottom', xytext = offset, textcoords='offset pixels',
#                 arrowprops = {'shrinkA':1,'arrowstyle':'simple', 'color':'black'},
#                 bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2))
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
    return None

def xy_lab(xy, text, ax, offset = (0,0), lw=1, fontsize=10, bbox=True, fc='white', ec='black'):
    if bbox:
        ax.annotate(text=text, xy=xy, ha='center', va = 'bottom', xytext = offset, textcoords='offset pixels', fontsize = fontsize, 
                    bbox=dict(boxstyle="square,pad=0.3", fc=fc, ec=ec, lw=lw))
    else:
        ax.annotate(text=text, xy=xy, ha='center', va = 'bottom', xytext = offset, textcoords='offset pixels', fontsize = fontsize)

def dir_arrow(ax, x, y, dx, dy, arrow_length, text, fontsize=10):
    lw = 1
    ax.annotate(text, xy=(x, y), xytext=(x, y-arrow_length),
                ha='center', va='center', fontsize=fontsize, #rotation =45,
                xycoords=ax.transAxes,
               bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=lw))

    ### gw flow direction arrow approximation
    ax.annotate('', xy=(x+dx, y+dy), xytext=(x, y),
                arrowprops=dict(facecolor='black', width=1.5, alpha=1, headwidth=5),
                ha='center', va='center', fontsize=10,xycoords=ax.transAxes, 
               )
    return None

def plt_cln(ax, nbins=5, label=True):
    """ Set basic xy axis labels, reduce axis clutter"""
    ax.ticklabel_format(style='plain')
    if label:
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
    # adjust major locators
    ax.set_yticklabels(labels=ax.get_yticklabels(), 
                       rotation=90, verticalalignment = "center");
    ax.locator_params(axis='x', nbins=nbins);
    ax.locator_params(axis='y', nbins=nbins);

    return None

def plt_arrow(ax, xoff = 0.7, yoff=0.15):
    x, y, arrow_length = xoff, yoff, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20, 
                xycoords=ax.transAxes)
    return None

def make_multi_scale(ax, xoff,yoff, dist = 1E3, scales = [4,2,1], units = 'km'):
    """ Plot legend with multi distances (3) at fractional axes offset
    A typical scalebar has distance and two time distances and half the distance
    With labels for just the distance and two times
    Input:
    ax: matplotlib subplot axis
    xoff, yoff: = x and y off fractional offset to ax
    dist: distance in 
    """
    minx, maxx = ax.get_xlim()
    miny, maxy = ax.get_ylim()
    # translate to actual coordinates
    lx = (maxx-minx)*xoff + minx
    ly = (maxy-miny)*yoff + miny
    height = (maxy-miny)*0.02 # height of bar
    if units=='km':
        km = 1E3
    elif units=='m':
        km = 1
    for n, adj in enumerate(scales):
        # alternate colors black, white, black
        color='black'
        if n%2==1:
            color='white'
        rect = box(lx, ly, lx +dist*adj, ly+height)
        gpd.GeoSeries(rect).plot(color=color, edgecolor='black', ax=ax)
        # round to 1 decimal, and drop decimal if it is only a .0
        dist_lab = str(np.round(dist*adj/km, 1)).replace('.0','')
        if n ==0:
            dist_lab += ''+units

        # use ha = 'center' to center under lines
        ax.annotate(dist_lab, (lx + dist*adj , ly-2*height), xycoords='data', ha='center')
    return None