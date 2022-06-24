"""
map_cln module. 
Different functions for refining plotting spatial data with general python functions
First iteration as a Module June 2022
Author: Andrew Calderwood
"""

import matplotlib.pyplot as plt
import geopandas as gpd


def gdf_2_lims(gdf,ax):
	""" Use the x-y bounds of a geodataframe to set the axes limits on a plot"""
	xmin, ymin = gdf.geometry.bounds.loc[:,['minx','miny']].min()
	xmax, ymax = gdf.geometry.bounds.loc[:,['maxx','maxy']].max()
	ax.set_xlim(xmin,xmax)
	ax.set_ylim(ymin,ymax)
