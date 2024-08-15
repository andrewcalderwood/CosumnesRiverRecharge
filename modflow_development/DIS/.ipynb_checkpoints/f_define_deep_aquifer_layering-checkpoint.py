''' function to define the deep aquifer layering in the Cosumnes model
using an elevation threshold believed to define the outcrop of the less permeable geology
and a slope of the dip
First iteraton: 2024-8-2
'''
import numpy as np
import matplotlib.pyplot as plt

from shapely.ops import LineString, polygonize, unary_union
from shapely.geometry import Polygon

import rasterio
from rasterio.features import shapes, rasterize

def get_rast_poly(alldata, maxl, n):
    '''Given object coords return the polygon object
    '''
    maxl2 = np.where(maxl==np.sort(maxl)[-n])[0][0]
    # verify the correct vectorized object is selected
    temp = alldata[maxl2][0].get('coordinates')[0]
    tl = Polygon(temp)
    return(tl)

def define_deep_aquifer_layering(top_botm, dem_data, delr=200, 
                                 cutoff_elev=56, slope=0.02, plt=False):
    '''Function to assign cells to have the low conductivity geology
    botm: 3D array of nlay, nrow, ncol
    dem_data: 2D array of nrow,ncol with ground surface elevations
    cutoff_elev: ground surface elevation associated with outcropping 
        of lower conductivity geology
    slope: The bottom boundary has a dip of 1-2 degrees which is essentially
        a slope of 0.015 based on given cross section data

    '''
    nlay, nrow,ncol = top_botm.shape
    nlay -= 1 # subtract 1 to account for top layer
    # Simplified ibound, only no flow cell if it is below the bottom of the Mehrten Formation
    # Specify no flow boundary based on rough approx of geology (upper basin volcanics)
    ibound = np.ones([nlay, nrow,ncol])
    
    ibound = ibound*(dem_data<cutoff_elev)
    # plot surface definition of boundary
    if plt:
        plt.imshow(ibound[0,:,:])
        plt.show()    
    
    # The function shapes from rasterio requires uint8 format
    ibound_line = ibound.astype(rasterio.uint8)
    out = shapes(ibound_line,connectivity = 8)
    alldata = list(out)
    
    maxl = np.zeros(len(alldata))
    for i in np.arange(0,len(alldata)):
        maxl[i] = len(alldata[i][0].get('coordinates')[0])
    
    # the two largest objects will be the two primary sides of the split
    maxl1 = np.where(maxl==np.sort(maxl)[-1])[0][0]
    
    # the line string object of
    tl1 = get_rast_poly(alldata, maxl, 1)
    tl2 = get_rast_poly(alldata, maxl, 2)
    # determine the smaller of the polygon objects to identify the mountain side
    mtn_id = np.argmin((tl1.area, tl2.area))+1
    
    # # now create a linestring object of the smaller polygon
    maxl2 = np.where(maxl==np.sort(maxl)[-mtn_id])[0][0]
    temp = alldata[maxl2][0].get('coordinates')[0]
    tl = LineString(temp)
    
    
    # Get the constant head or general head boundary after the no flow cells
    linerast = rasterio.features.rasterize([tl], out_shape = np.array((nrow,ncol)))
    # remove far east bound line
    linerast[:,ncol-1] = 0
    fix_bound = np.min(np.argwhere(linerast[0,:]==1))
    linerast[0,:] = 0
    linerast[0,fix_bound]
    
    # Set the polygon/raster for the top layer, no buffer needed
    poly0 = Polygon(temp)
    polyrast0 = rasterio.features.rasterize([poly0], out_shape = np.array((nrow,ncol)))
    
    ibound = np.ones([nlay, nrow,ncol])
    
    # The bottom boundary has a dip of 1-2 degrees which is essentially a slope of 0.015 based on given cross section data
    # calculate average layer thickness to inform step of deep geology
    avg_thick = abs(np.diff(top_botm[:,polyrast0==1],axis=0).mean(axis=1))
    run_total = 0
    for i in np.arange(1,nlay):
        # slope of bottom boundary
        run = (avg_thick[i-1]/slope)/delr
        run_total += run
        polyi = Polygon(temp).buffer(distance = run_total)
        polyrast = rasterio.features.rasterize([polyi], out_shape = np.array((nrow,ncol)))
        # Need to decide whether all layers or just the top layer are affected by ibound from elevation
        ibound[i,polyrast==1] = 0
    
    # wherever the constant head/specified head bound is the cells need to be active
    ibound[0,dem_data>cutoff_elev] = 0
    ibound[0,linerast==1] = 1
    # plot new refined surface boundary
    if plt:
        plt.imshow(ibound[0,:,:])
        plt.colorbar(shrink=0.5)
        plt.show()
    
    return(ibound)