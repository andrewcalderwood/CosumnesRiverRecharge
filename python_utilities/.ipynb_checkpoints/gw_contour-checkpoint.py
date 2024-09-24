import numpy as np
from osgeo import osr
from osgeo import ogr
from osgeo import gdal

import rasterio
from rasterio.transform import Affine


def export_kde_raster(Z, XX, YY, min_x, max_x, min_y, max_y, proj, filename):
    '''Export and save a kernel density raster.'''
    # Get resolution
    xres = (max_x - min_x) / len(XX)
    yres = (max_y - min_y) / len(YY)
    # Set transform
    transform = Affine.translation(min_x - xres / 2, min_y - yres / 2) * Affine.scale(xres, yres)
    # Export array as raster
    with rasterio.open(filename, mode = "w", driver = "GTiff",
            height = Z.shape[0], width = Z.shape[1],
            count = 1,dtype = Z.dtype,
            crs = proj,transform = transform,
    ) as new_dataset:
            new_dataset.write(Z, 1)


def raster2contours(raster_path, contour_path, contour_intervals):
    """Given a raster create a shapefile with contours 
    contour_intervals is list with elevations for contours"""
    #Open tif file as select band
    rasterDs = gdal.Open(raster_path)
    rasterBand = rasterDs.GetRasterBand(1)
    proj = osr.SpatialReference(wkt=rasterDs.GetProjection())

    #Get elevation as numpy array
    elevArray = rasterBand.ReadAsArray()
    #get dem max and min
    # demMax = base_round(elevArray.max(), base=contour_interval)
    # demMin = base_round(elevArray.min(), base=contour_interval)

    contourDs = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(contour_path)

    #define layer name and spatial 
    contourShp = contourDs.CreateLayer('contour', proj)

    #define fields of id and elev
    fieldDef = ogr.FieldDefn("ID", ogr.OFTInteger)
    contourShp.CreateField(fieldDef)
    fieldDef = ogr.FieldDefn("elev", ogr.OFTReal)
    contourShp.CreateField(fieldDef)

    #Write shapefile using noDataValue
      #Generate Contourlines - ContourGenerate from https://gdal.org/python/osgeo.gdal-module.html#ContourGenerate
    gdal.ContourGenerate(
        rasterBand,    #Band srcBand
        0,      #double contourInterval - This defines contour intervals
        0,      #double contourBase
        contour_intervals,      #int fixedLevelCount, defined contours rather than constant interval
        0,      #int useNoData
        0,      #double noDataValue
        contourShp, #Layer dstLayer
        0,      #int idField
        1       #int elevField
        )

    contourDs.Destroy()
