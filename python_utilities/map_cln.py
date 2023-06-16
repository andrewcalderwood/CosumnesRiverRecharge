"""
map_cln module. 
Different functions for refining plotting spatial data with general python functions
First iteration as a Module June 2022
Author: Andrew Calderwood
"""

import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import geopandas as gpd
import pandas as pd

def gdf_bnds(gdf, ax, buf=1):
    """ Use the x-y bounds of a geodataframe to set the axes limits on a plot"""
    gdf_bnd = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [gdf.unary_union.buffer(buf)], crs=gdf.crs)
    minx, miny, maxx, maxy = gdf_bnd.bounds.values[0]
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)
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

def plt_cln(ax):
    """ Set basic xy axis labels, reduce axis clutter"""
    ax.ticklabel_format(style='plain')
    plt.xlabel('Easting (m)')
    plt.ylabel('Northing (m)')
    # adjust major locators
    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    plt.yticks(rotation=90, verticalalignment = "center")
