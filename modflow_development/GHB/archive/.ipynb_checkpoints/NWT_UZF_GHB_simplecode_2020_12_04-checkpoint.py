# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Cosumnes Model 
# @author: Andrew

# %%
# standard python utilities
import os
import sys
import glob
import pandas as pd
import numpy as np
import calendar
import time

# standard python plotting utilities
import matplotlib as mpl
import matplotlib.pyplot as plt

# standard geospatial python utilities
import pyproj # for converting proj4string
import shapely
import geopandas as gpd
import rasterio

# mapping utilities
import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %%

# run installed version of flopy or add local path
try:
    import flopy
    from flopy.discretization.structuredgrid import StructuredGrid
    from flopy.utils.reference import SpatialReference
    from flopy.utils import Raster
except:
    import flopy
    fpth = os.path.abspath(os.path.join('..', '..'))
    sys.path.append(fpth)
    from flopy.discretization.structuredgrid import StructuredGrid
    from flopy.utils.reference import SpatialReference
    from flopy.utils import Raster
from flopy.utils.gridgen import Gridgen
from flopy.utils import OptionBlock
import flopy.utils.binaryfile as bf


print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('matplotlib version: {}'.format(mpl.__version__))
print('flopy version: {}'.format(flopy.__version__))

# %%

# Transient -> might want to think about making SP1 steady
end_date = '2019-12-31'
# end_date = '2018-01-02'
strt_date = '2018-01-01'

dates = pd.date_range(strt_date, end_date)

# The number of periods is the number of dates 
nper = len(dates)

# Each period has a length of one because the timestep is one day
perlen = np.ones(nper)
# Steady or transient periods
steady = np.zeros(nper)
steady[0] = 1 # first period is steady state, rest are transient
steady = steady.astype('bool').tolist()
# Reduce the number of timesteps to decrease run time
nstp = np.ones(nper)*6

# %%
#Maribeth's model parameters, had to switch nrow and ncol due to her issue in xul, yul
nrow=100
ncol=230
delr=200
delc=200
rotation=52.9

# The number of layers should be 1 for the Mehrten formation, 1 for the laguna plus the number of TPROGS layers,
# where the Laguna formation will be clipped by the TPROGS layers
num_tprogs = 10
nlay = 2 + num_tprogs
tprog_thick = 100/num_tprogs

# There is essentially no difference bewtween WGS84 and NAD83 for UTM Zone 10N
# proj4_str='EPSG:26910'
proj4_str='+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '

# %%
## Set up directory referencing
# Package data
gwfm_dir = os.path.dirname(os.path.dirname(os.getcwd()))

# %%
from flopy.utils.geometry import Polygon, LineString, Point
# Original model domain, 44.7 deg angle
# m_domain = gpd.read_file(gwfm_dir+'\\GWModelDomain_UTM10N\\GWModelDomain_Rec_UTM10N.shp')
# New model domain 52.9 deg
m_domain = gpd.read_file(gwfm_dir+'\\NewModelDomain\\GWModelDomain_52_9deg_UTM10N_WGS84.shp')

# Need to check this when changing model domains
xul, yul = list(m_domain.geometry.values[0].exterior.coords)[1]
list(m_domain.geometry.values[0].exterior.coords)
# m_domain.geometry.values[0].exterior

# %% [markdown]
# According to Flopy GitHub "Technically you need to create both a SpatialReference object and a ModelGrid object, but in practice the code looks very similar and can easily be implemented in one line."
# WGS84 Zone 10N has EPSG: 32610  
# Lower left corner of model is   
# Zone 10 N  
# Easting: 661211.18 m E  
# Northing: 4249696.50 m N  
# angle is approximate 53 degrees  

# %%
m = flopy.modflow.Modflow(modelname = 'MF', exe_name = 'MODFLOW-NWT', 
                          version = 'mfnwt', model_ws='data')
# m = flopy.modflow.Modflow(modelname = 'MF', exe_name = 'mf2005', 
#                           version = 'mf2005', model_ws='data')
#lenuni = 1 is in ft, lenuni = 2 is in meters
# itmuni is time unit 5 = years, 4=days, 3 =hours, 2=minutes, 1=seconds
dis = flopy.modflow.ModflowDis(nrow=nrow, ncol=ncol, 
                               nlay=nlay, delr=delr, delc=delc,
                               model=m, lenuni = 2, itmuni = 4,
                               xul = xul, yul = yul,rotation=rotation, proj4_str=proj4_str,
                              nper = nper, perlen=perlen, nstp=nstp, steady = steady,
                              start_datetime = strt_date)


# %%
# m.modelgrid.set_coord_info(xoff=xoff, yoff=yoff, proj4='EPSG:32610', angrot=angrot)
mg = m.modelgrid
# Write model grid to shapefile for later use
# mg.write_shapefile(gwfm_dir+'/DIS_data/grid/grid.shp', epsg = '32610')
# mg.write_shapefile(gwfm_dir+'/DIS_data/44_7_grid/44_7_grid.shp', epsg = '32610')


# %%
sensors = gpd.read_file(gwfm_dir+"/DEM_data/allsensor_latlong.csv")
sensors.Latitude = sensors.Latitude.astype(np.float64)
sensors.Longitude = sensors.Longitude.astype(np.float64)
sensors.geometry = gpd.points_from_xy(sensors.Longitude, sensors.Latitude)
sensors.crs = 'epsg:4326'
sensors = sensors.to_crs('epsg:32610')

# %%
# import contextily as ctx
# fig, ax = plt.subplots(figsize = (10,10))
# mg.plot(ax=ax, alpha = 0.3)
# sensors.plot(ax =ax)

# ctx.add_basemap(ax, source = ctx.providers.Esri.WorldImagery, crs='epsg:26910', alpha = 0.6)
# # ctx.add_basemap(ax, source = ctx.providers.Esri.WorldStreetMap, crs=tg.crs.to_string())

# plt.xlabel('Easting (m)')
# plt.ylabel('Northing (m)')
# plt.savefig('Plots/Model_SFR_UZF_progress/Model_grid_overlay_satellite.png', dpi = 300)

# %%
# Load model grid as geopandas object
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
# grid_p = gpd.read_file(gwfm_dir+'/DIS_data/44_7_grid/44_7_grid.shp')


# Find Michigan Bar location
mb_gpd = sensors[sensors.Sensor_id == "MI_Bar"]
mb_grid = gpd.sjoin(mb_gpd, grid_p, how = 'left', op = 'intersects')

# %%
# Get vertexes of model domain
# ll = mg.get_coords(0, 0) #lower left
# lr = mg.get_coords(nrow*delr, 0) #lower right
# ur = mg.get_coords(nrow*delr, ncol*delc) #upper right
# ul = mg.get_coords(0, ncol*delc) #upper left
ll = mg.get_coords(0, 0) #lower left
lr = mg.get_coords(0, nrow*delr) #lower right
ur = mg.get_coords(ncol*delc, nrow*delr) #upper right
ul = mg.get_coords(ncol*delc, 0) #upper left
print(ll, lr, ur, ul)

# Shapefile of model bounds
from shapely.geometry import Polygon
vertices = np.stack(np.asarray((ll,lr, ur, ul)))
vertices
geoms = Polygon(vertices)

# %% [markdown]
# ### Raster files can be loaded using the `Raster.load` method

# %%
# Full size dem of northern sac valley
# raster_name = gwfm_dir+"/DEM_data/USGS_ten_meter_dem/transformed.tif"

# rio10_utm = Raster.load(raster_name)

# %%
# rio10_utm.plot()
# t0 = time.time()
# rio10_utm.crop(vertices, invert=False)
# crop_time = time.time() - t0
# rio10_utm.plot()


# %%
# Nearest neighbor determines the nearest pixel and assumes its value
# linear is as it sounds, cubic is the smoothed version of linear essentially by using a cubic function
# the linear method takes a very, very long time - an hour plus??, just stick with nearest
# nearest takes 170.209 seconds - 220 seconds
# the linear interpolation causes the program to crash
# t0 = time.time()
# dem_data = rio10_utm.resample_to_grid(m.modelgrid.xcellcenters,
#                                 m.modelgrid.ycellcenters,
#                                 band=rio10_utm.bands[0],
#                                 method="nearest")
# resample_time = time.time() - t0
# print("Resample time, nearest neighbor: {:.3f} sec".format(time.time() - t0))

# %%
# np.savetxt(gwfm_dir+'\DIS_data\dem_44_7_200m_nearest.tsv', dem_data, delimiter = '\t')

# Based on Maribeth's grid aligned with Alisha's TPROGS model
# dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_52_9_200m_nearest.tsv', delimiter = '\t')
dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_52_9_200m_linear.tsv', delimiter = '\t')
# dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_44_7_200m_linear_missing_right_corner.tsv', delimiter = '\t')

# dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_44_7_200m_nearest.tsv', delimiter = '\t')

import seaborn as sns
sns.heatmap(dem_data, cmap = 'viridis', vmin = 0)

# %% [markdown]
# # Capture cross section of deeper geology

# %%
# Pre-cretaceous metamorphic rocks - (variable thickness 200-500 ft thick)
# Ione formation (200 ft thick)
# Valley Springs formation (300 ft thick)
# Mehrten Formation (100 ft thick to 300 ft thick) (1-2 deg dip)
# Laguna Formation (less than 100 ft to between 200-300 ft thick) (less than 1 deg dip)
# upper formation (informed by well logs) (100 ft)
# ibound < 0 is constant head
# ibound = 0 is inactive cell
# ibound > 0 is active cell

# %% [markdown]
# ## Need to adjust for change in model grid, based on Michigan Bar previously, maybe also look at effect of model domain angle vs cross section angle

# %%
mb_grid
# The stream gage at michigan bar is now 13 columns in from the boundary
# mehrtenbound

# Cross section E appears to have an angle of 0 compared to the model domain,

# %%
# # columns are xtop_miles, ytop_ft_amsl, xbot_miles, ytop_ft_amsl
# # XS upper bound should be at Michigan bar which is between Jackson road and Sacramento-Amador county line split
# # Mile 36 is approximately where Michigan bar aligns with the cross section
MB_XS_mile = 36
mehrtenbound = pd.read_csv(gwfm_dir+'/DIS_data/Mehrten_boundary_x_y.csv', parse_dates = False, 
                index_col = False, sep = ',', header = 'infer')
# Convert miles to feet and sets x value based on location of Michigan bar
# 0 is michigan bar and each major change in geologic dip is based on distance from Michigan Bar
mehrtenbound.xtop_miles = -5280*(MB_XS_mile - mehrtenbound.xtop_miles)
mehrtenbound.xbot_miles = -5280*(MB_XS_mile - mehrtenbound.xbot_miles)
# No flod boundary based on the original coordinates of the bottom of the Mehrten formation
mehrtenbound.noflow_x_miles = -5280*(MB_XS_mile - mehrtenbound.noflow_x_miles)

# East of mile 32 the entire vertical cross section, including up to the near entire surface
# is composed of old geologic formations that are not water bearing
volcanic_bound = (MB_XS_mile - 32)*-5280
# noflow_ind = int((1-(volcanic_bound/sumx))*ncol)

# Plot the x and y values
fig, ax = plt.subplots(1, 1)
mehrtenbound.plot(x = 'xtop_miles', y = 'ytop_ft_amsl', ax = ax)
mehrtenbound.plot(x = 'xbot_miles', y = 'ybot_ft_amsl', ax = ax)
plt.plot(-100*3.28*np.arange(0,len(dis.top[40,:])), np.flip(3.28*dis.top[40,:]))
# print(mehrtenbound)

# %%

def xy_to_botm(xbound_ft, ybound_ft, nspace,ntransversespace):
    laybotm = np.zeros((ntransversespace, nspace))
    # Nspace will be either nrow or ncol depending model direction
    # ntransversespace is the opposite of nspace (ie nrow if nspace is ncol)
    # Calculate the distance between each major change in dip
    dx = np.diff(xbound_ft)
    # Scale by the total distance across the coordinates to get percentages
    sumx = np.sum(dx)
    dx /= sumx
    # Multiply the number of columns by the percent of columns in each section of constant dip
    dx *= nspace
    # Round the number of columns to allow proper use for indexing
    nx = np.round(dx).astype(int)
    # Fix any discrepancy in number of columns due to issues with rouding the percentages of columns
    # Add a column to the last set of columns because there is already uncertainty at the deeper end
    while(np.sum(nx)-nspace !=0):
        if np.sum(nx)-nspace <0:
            nx[-1] += 1
        elif np.sum(nx)-nspace >0:
            nx[-1] -= 1
    sum(nx)

    # Now split the coordinates into downsized coordinates in between each major change in dip
    k = 0
    for i in np.arange(0,len(nx)):
        for j in np.arange(0,ntransversespace):
            laybotm[j, k:k+nx[i]] = np.arange(ybound_ft[i],ybound_ft[i+1], -(ybound_ft[i]-ybound_ft[i+1])/nx[i])
        k += nx[i]
    return(laybotm)


# %%
# X (east -west) and y (up-down vertical) of major dip changes for Mehrten Formation top boundary

mehrten_top = xy_to_botm(mehrtenbound.xtop_miles,mehrtenbound.ytop_ft_amsl, ncol, nrow)
# X (east -west) and y (up-down vertical) of major dip changes for Mehrten Formation bottom boundary
# drop na is because there are less values to mark changes in the bottom than top boundary
mehrten_bottom = xy_to_botm(mehrtenbound.xbot_miles.dropna(),mehrtenbound.ybot_ft_amsl.dropna(), ncol, nrow)

# Original x,y data for Mehrten bottom boundary to represent the noflow bounds
no_flow_bound = xy_to_botm(mehrtenbound.noflow_x_miles.dropna(), mehrtenbound.noflow_y_ft_amsl.dropna(),ncol,nrow)

# %%
botm = np.zeros((nlay, nrow, ncol))
# Convert units from ft to meters and flip to match direction
botm[-2,:,:] = np.flip(mehrten_top/3.28)
botm[-1,:,:] = np.flip(mehrten_bottom/3.28)
no_flow_bound = np.flip(no_flow_bound/3.28)
dis.botm = botm
# dis.plot()

# %%
botm.shape, dem_data.shape

# %%
num_tprogs, nlay

# %% [markdown]
# ## Adjustment to bottom boundary to ensure sufficient top layer thickness for the TPROGS model
# Although the bottom boundaries are being artifically lowered to allow for sufficient layer thickness, this will be corrected when ibound is implemented based on where the actual bottom boundary is and where there is high elevations based on likelihood to be volcanics geology.

# %%
# The TPROGS model is 100m thick with some of it above the land surface
# to be safe, there should be at least 50 meters below ground surface to the bottom boundary

# Create TPROGS layers from bottom up
# tprog_thick = 100/num_tprogs
# botm[-3,:,:] = -80
# for i in np.arange(-4,-3-num_tprogs,-1):
#     botm[i,:,:] = botm[i+1,:,:] + tprog_thick
    
    
# Create TPROGS layers from top down
# tprog_thick = 100/num_tprogs
botm[0,:,:] = dem_data - tprog_thick
for i in np.arange(1,num_tprogs):
    botm[i,:,:] = botm[i-1,:,:] -tprog_thick
    
# Thickness to give to bottom layers below the TPROGS layers just to provide adequate spacing,
# this will be corrected by changing the geology in the layers above to account for what is actually in
# the Mehrten and what is in the Laguna formations, thickness of 5 also prevents any messy overlap
thickness_to_skip =10
# # Find where top boundary of Mehrten Formation rises within 10 meters of the top layer (10m for sufficient layer thickness)
bot3ind = np.min(np.where(botm[-2,:,:]>botm[-3,:,:]- thickness_to_skip)[1])

# # Where the top boundary of Mehrten was within 10 meters of the top layer 
# # set it equal to top layer elevation minus 10 for sufficient layer thickness
botm[-2,:,bot3ind:] = botm[-3,0,bot3ind]- thickness_to_skip
# # Repeat steps above for bottom of Mehrten formation with the top of the Mehrten formation
bot3ind = np.min(np.where(botm[-1,0,:]>botm[-2,0,:]- thickness_to_skip))
botm[-1,:,bot3ind:] = botm[-2,0,bot3ind]-thickness_to_skip

# %%

fig,ax = plt.subplots(figsize = (12,12))
plt.plot(dem_data[8,:])

for i in np.arange(0,nlay):
    plt.plot(botm[i,8,:])

# %%
# Set the elevation of the top layer based on the DEM
m.dis.top = dem_data
# Bottom of model based on geology
m.dis.botm = botm
chk = dis.check()
chk.summary_array

# %% [markdown]
# ## Import soil data for Lake Package, UZF Package, SFR Package hydraulic parameters

# %%
mb_name = gwfm_dir+"/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp"

mb = gpd.read_file(mb_name)
mb = mb.to_crs('epsg:32610')

# %%
uzf_path = gwfm_dir+'\\UZF_data'
soil_path = uzf_path+'\\wss_gsmsoil_CA'
# # Read in the soil map spatial data
# soil_gpd = gpd.read_file(uzf_path+'\\wss_gsmsoil_CA\\spatial\\gsmsoilmu_a_ca.shp')
# soil_gpd = soil_gpd.to_crs('EPSG:32610')
# # soil_gpd.plot()

# %%
# # Intersect state soil map with model domain
# gpd_mb = gpd.overlay(soil_gpd, mb, how = 'intersection')
# # gpd_mb.plot(column = 'MUKEY')

# %%
# # Read in files relevant to water content, sat K, and epsilon that correspond to the spatial data
# por = pd.read_csv(uzf_path+'\\wss_gsmsoil_CA\\soildb_US_2003_Porosity.txt',sep=',', header = 0)
# ksat = pd.read_csv(uzf_path+'\\wss_gsmsoil_CA\\soildb_US_2003_Ksat.txt',sep=',', header = 0)
# psdi = pd.read_csv(uzf_path+'\\wss_gsmsoil_CA\\soildb_US_2003_PSDI.txt',sep=',', header = 0)
# # Pull out representative values for porosity, saturated hydraulic conductivity and pore size distribution index
# por = por.drop(labels = ['Porosity_Low', 'Porosity_High'], axis = 1)
# ksat = ksat.drop(labels = ['Ksat_Low', 'Ksat_High'], axis = 1)
# psdi = psdi.drop(labels = ['PSDI_Low', 'PSDI_High'], axis = 1)

# %%
# # Prepare the column name for joining spatial data with the reference data
# gpd_mb = gpd_mb.rename(columns={'MUKEY':'mukey'})
# # gpd_mb, por
# # gpd_mb.join(por, on = 'mukey')
# soil = por.merge(ksat, on = 'mukey', how = 'inner')
# soil = soil.merge(psdi, on = 'mukey', how = 'inner')
# # Calculate epsilon using the pore size distribution index equation from Brooks and Corey 1960
# soil['EPS'] = (2 + 3*soil.PSDI_Rep)/soil.PSDI_Rep
# gpd_mb.mukey = gpd_mb.mukey.values.astype('int64')
# gpd_mb = gpd_mb.merge(soil, on = 'mukey')


# %%
# # Samples the points every 10 meters to match with the 100 meter grid
# grid_uzf = gpd.sjoin(gpd_mb, grid_p,   how = "inner", op= "intersects")


# %%
# Dissolves the points every 10 meters to the 100 meter spacing
# grid_uzf = grid_uzf.dissolve(by = 'node', aggfunc = 'mean')
# length of 70,073 before dissolve, 66,000 after dissolve

# %%
# # Convert Ksat from micrometers/sec to m/d
# grid_uzf.Ksat_Rep *= (86400*(1E-6))

# %% [markdown]
# ### Write grid_uzf to shapefile to avoid having to repeat analysis

# %%
# grid_uzf.to_file(uzf_path+'/final_grid_uzf/griduzf.shp')
grid_uzf = gpd.read_file(uzf_path+'/final_grid_uzf/griduzf.shp')


# %%
def fill_uzf(uzfvalues, grid_uzf):
    # convert geopandas object to regular np array for soil data
    temp = np.zeros((nrow,ncol))
    temp[(grid_uzf.row.values-1).astype(int),(grid_uzf.column.values-1).astype(int)] = uzfvalues
    return(temp)


# %%
soilKs_array = np.loadtxt(uzf_path+'/final_soilKs.tsv', delimiter = '\t')
soiln_array = np.loadtxt(uzf_path+'/final_soiln.tsv', delimiter = '\t')
soileps_array = np.loadtxt(uzf_path+'/final_soileps.tsv', delimiter = '\t')

# soilKs_array = fill_uzf(grid_uzf.Ksat_Rep, grid_uzf)
# soiln_array = fill_uzf(grid_uzf.Porosity_R, grid_uzf)
# soileps_array = fill_uzf(grid_uzf.EPS, grid_uzf)

# np.savetxt(uzf_path+'/final_soilKs.tsv', soilKs_array, delimiter = '\t')
# np.savetxt(uzf_path+'/final_soiln.tsv', soiln_array, delimiter = '\t')
# np.savetxt(uzf_path+'/final_soileps.tsv', soileps_array, delimiter = '\t')

# %% [markdown]
# # SFR

# %%
sfr_dir = gwfm_dir+'/SFR_data/'

# %%
# mb_name = gwfm_dir+"/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp"
# # mb_name = gwfm_dir+"/GWModelDomain_UTM10N/GWModelDomain_Rec_UTM10N.shp"

# # Rivers and creeks in the larger area encompassing Cosumnes River in both South American and Cosumnes Subbasins
# rivers = gpd.read_file(gwfm_dir+"/SFR_data/Sac_valley_rivers/Sac_valley_rivers.shp")
# mb = gpd.read_file(mb_name)
# mb = mb.to_crs('epsg:32610')
# rivers = rivers.to_crs('EPSG:32610')
# rivers_clip = gpd.clip(rivers, mb)
# rivers_clip.plot()
# # rivers_clip.GNIS_Name.unique()

# %%
# print(rivers_clip.columns)
# # Split into individual streams/creeks
# cr_ind = rivers_clip[rivers_clip.GNIS_Name == 'Cosumnes River']
# dc_ind = rivers_clip[rivers_clip.GNIS_Name == 'Deer Creek']
# cc_ind = rivers_clip[rivers_clip.GNIS_Name ==  'Coyote Creek']
# # Pull out data for each river/creek
# cr = rivers_clip.loc[cr_ind.index,]
# dc = rivers_clip.loc[dc_ind.index,]
# cc = rivers_clip.loc[cc_ind.index,]


# %%
# import fiona
# from shapely.geometry import shape, mapping
# from shapely.ops import linemerge

# cr.geometry.values.crs = "epsg:32610"
# geom = linemerge(cr.geometry.values)
# # how often to interpolate a point
# dline = 10
# # # length of the LineString
# length = int(geom.length)
# point = np.zeros((int(length/dline)+1,3))
# for i, distance in enumerate(range(0, int(length), dline)):
#          point[i,:] = geom.interpolate(distance).coords[:][0]
# point = point[:,[0,1]]
# plt.plot(point[:,0],point[:,1])


# %%
# raster_name = gwfm_dir+'/DEM_data/USGS_ten_meter_dem/regional_10m.tif'

# # dem10 = rasterio.open(raster_name)

# pnts = pd.DataFrame()
# with rasterio.open(raster_name) as src:
#     pnts['z'] = [sample[0] for sample in src.sample(point)]
# pnts
# pnts['easting'] = point[:,0]
# pnts['northing'] = point[:,1]
# pnts = pnts[pnts.z > -1E4]



# %%
# plt.figure(figsize = (10,10))
# pnts.z.plot()
# pnts['slope'] = 0.002
# for i in np.arange(1,len(pnts)):
#     if pnts.z.values[i] >= pnts.z.values[i-1]:
#         # if strtop is greater than previous strtop use previous elevation minus the average slope
#         slope = ((np.max(pnts.z) - np.min(pnts.z))/geom.length)*dline
#         if pnts.index[i] < 800:
#             slope = 0.002
#         elif pnts.index[i] < 2700:
#             slope = 0.0003
#         elif pnts.index[i] < 3200:
#             slope = 0.001
#         else:
#             slope = 0.0003
#         pnts.z.values[i] = pnts.z.values[i-1] - slope*dline
#         pnts.slope.values[i] = slope
# pnts.z.plot()

# %%
#  pnts['Point_order'] = pnts.index

# %%
# pnts_gpd = gpd.GeoDataFrame(pnts, geometry = gpd.points_from_xy(pnts.easting, pnts.northing))
# pnts_gpd.crs = 'epsg:32610'

# %%
# plot_lines(m.modelgrid.grid_lines)
# Write model grid to shapefile of polygons
# m.modelgrid.write_shapefile('grid/grid.shp', epsg = '32610')
# grid_p = gpd.read_file('grid/grid.shp')

# %%
# # Samples the points every 10 meters to match with the 100 meter grid
# grid_sfr = gpd.sjoin(grid_p, pnts_gpd, how = "inner", op= "intersects")
# grid_sfr.plot()

# %%
# # Set reach length for each reach based on the separation used to create the points from the line object
# # dline is 10 meters
# grid_sfr['length_m'] = dline

# # Dissolve the points again but using sum this time to get the total length of each reach
# temp = grid_sfr.dissolve(by = 'node', aggfunc = 'sum').length_m.values

# # Dissolves the points every 10 meters to the 200 meter spacing, using mean because the interested component is elevation
# grid_sfr = grid_sfr.dissolve(by = 'node', aggfunc = 'mean')
# grid_sfr.length_m = temp

# grid_sfr = grid_sfr.sort_values(by = 'Point_order')
# grid_sfr['reach'] = np.arange(1,len(grid_sfr)+1)

# %%
# grid_sfr.to_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')
grid_sfr = gpd.read_file(sfr_dir+'/final_grid_sfr/grid_sfr.shp')


# %%
# floodplain sensor info
fp_sensors = gpd.read_file(gwfm_dir+"/LAK_data/floodplain_logger_metadata.csv", header = True)
fp_sensors.Northing = fp_sensors.Northing.astype(np.float64)
fp_sensors.Easting = fp_sensors.Easting.astype(np.float64)
fp_sensors.geometry = gpd.points_from_xy(fp_sensors.Easting, fp_sensors.Northing)
fp_sensors.crs = 'epsg:32610'

# %%
# od_breach is the sensor location where the breach was made in the levees for flow to leave the river
od_breach = fp_sensors[fp_sensors['Logger Location']=='OD_Excavation']
od_breach


# %%
# first need to try with just adding extra cross section
# append the extra reach to the list of reaches and resort by reach
# grid_sfr = grid_sfr.append(grid_breach).sort_values(by = 'reach')
# # next, need to relabel reaches to account for the added reach witha duplicate number
# grid_sfr.reach = np.arange(1,len(grid_sfr)+1)
# grid_sfr

# %%
# # Buffer the location of the breach sensor to have overlap with the river streamline
# just sjoin the geometry because the extra info is unnecessary
# spatial join breach sensor polygon with sfr grid locations to find match
grid_breach = gpd.sjoin(grid_sfr, 
                        gpd.GeoDataFrame(geometry = od_breach.geometry.buffer(25), crs = 'epsg:32610'), how = "inner", op= "intersects")
# add a reach to the overlap cell that will be used to divert flow (there will be two reaches in one cell)
grid_breach

# %%
XS8pt = pd.read_csv(sfr_dir+'8pointXS.csv')
XSlocs = gpd.read_file(sfr_dir+'8pointXS_locs\\8pointXS_locs.shp')
XSlocs.crs = 32610

XSg  = gpd.sjoin(grid_sfr, XSlocs, how = "inner", op= "contains", lsuffix = 'sfr',rsuffix = 'xs')
# Append the grid_breach location to the list of cross sections to split the segment
XSg = XSg.append(grid_breach).sort_values('reach')
# Copy the XS site name from the previous last site to the breach site to keep same XS
XSg.Site.iloc[-1] = XSg.Site.iloc[-2]
len(XSg), len(XS8pt.loc[0,:])/2

# %% [markdown]
# ### Read in 8 pt XS, revised by simplifying from Constantine 2001

# %%
# There is one reach for each cell that a river crosses
NSTRM = -len(grid_sfr)
# There should a be a stream segment if there are major changes
# in variables in Item 4 or Item 6
# 1st segment is for the usgs Michigan Bar rating curve, one for each XS, plus 2 for the floodplain diversion
NSS = 1 + len(XSg) 
# nparseg (int) number of stream-segment definition with all parameters, must be zero when nstrm is negative
NPARSEG = 0
CONST = 86400 # mannings constant for SI units, 1.0 for seconds, 86400 for days
# real value equal to the tolerance of stream depth used in
# computing leakage between each stream reach and active model cell
DLEAK = 0.0001 # unit in lengths, 0.0001 is sufficient for units of meters
IPAKCB = 55
# writes out stream depth, width, conductance, gradient when cell by cell
# budget is specified and istcb2 is the unit folder
ISTCB2 = 54
# specifies whether unsat flow beneath stream or not, isfropt 2 has properties read for each reach, isfropt 3 also has UHC
# read for each reach, isfropt 4 has properties read for each segment (no UHC), 5 reads for each segment with UHC
ISFROPT = 1
# nstrail (int), number of trailing weave increments used to represent a trailing wave, used to represent a decrease 
# in the surface infiltration rate. Can be increased to improve mass balance, values between 10-20 work well with error 
# beneath streams ranging between 0.001 and 0.01 percent, default is 10 (only when isfropt >1)
NSTRAIL = 20
# isuzn (int) tells max number of vertical cells used to define the unsaturated zone beneath a stream reach (default is 1)
ISUZN = 1
#nsfrsets (int) is max number of different sets of trailing waves (used to allocate arrays), a value of 30 is sufficient for problems
# where stream depth varies often, value doesn't effect run time (default is 30)
NSFRSETS = 30
# IRTFLG (int) indicates whether transient streamflow routing is active, must be specified if NSTRM <0. If IRTFLG >0 then
# flow will be routed with the kinematic-wave equations, otherwise it should be 0 (only for MF2005), default is 1
IRTFLG = 1
# numtim (int) is number of sub time steps used to route streamflow. Streamflow time step = MF Time step / NUMTIM. 
# Default is 2, only when IRTFLG >0
NUMTIM = 4
# weight (float) is a weighting factor used to calculate change in channel storage 0.5 - 1 (default of 0.75) 
WEIGHT = 0.75
# flwtol (float), flow tolerance, a value of 0.00003 m3/s has been used successfully (default of 0.0001)
FLWTOL = 0.0001

sfr = flopy.modflow.ModflowSfr2(model = m, nstrm = NSTRM, nss = NSS, nparseg = NPARSEG, 
                           const = CONST, dleak = DLEAK, ipakcb = IPAKCB, istcb2 = ISTCB2, 
                          isfropt = ISFROPT, nstrail = NSTRAIL, isuzn = ISUZN, irtflg = IRTFLG, 
                          numtim = NUMTIM, weight = WEIGHT, flwtol = FLWTOL,
                                reachinput=True, transroute=True, tabfiles=True,
                                tabfiles_dict={1: {'numval': nper, 'inuit': 56}})

# %%
# Add option block at the top of the sfr input file for tabfiles
tab_option = flopy.utils.OptionBlock(options_line = ' reachinput transroute tabfiles 1 ' + str(nper), package = sfr, block = True)
sfr.options = tab_option
# sfr.options

# %%
# Modflow NWT additions to SFR package set up
# sfr.transroute = True
# sfr.reachinput = True
# sfr.tabfiles = True
# # numval is the number of values in the flow tab files, inuit is the corresponding unit file
# sfr.tabfiles_dict = {1: {'numval': nper, 'inuit': 56}}

# %%
xs_sfr = grid_sfr.set_index('reach')
# set all reaches to start as segment 1 which will be changed iteratively based on the number of cross-sections
xs_sfr['iseg'] = 1
# add a column reach_new that will be changed iteratively as the segment number is changed
xs_sfr['reach_new'] = xs_sfr.index
# xs_sfr

# %% [markdown]
# ## Define reach data based on ISFROPT

# %%
## Given the reach number of each XS, the 718 reaches will be broken down into each segment
## create a new reach column based on XS reach number and 

segcount = 2
for i in np.arange(0,len(XSg)):
    temp_reach = XSg.reach.values[i]
    rchnum = xs_sfr.index[-1] - temp_reach+1
    xs_sfr.reach_new.loc[temp_reach:] = np.linspace(1,rchnum, rchnum)
    xs_sfr.iseg.loc[temp_reach:] = segcount
    segcount +=1

# %%
xs_sfr.reach_new = xs_sfr.reach_new.astype(int)


# %%
# Determine which layer the streamcell is in
# since the if statement only checks whether the first layer is greater than the streambed elevation, 
# otherwise it would be less than and zero (most should be in layer 0)
sfr_lay = np.zeros(len(grid_sfr))

for i in np.arange(0,nlay-1):
    # pull out elevation of layer bottom
    lay_elev = botm[i, (grid_sfr.row.values-1).astype(int), (grid_sfr.column.values-1).astype(int)]
    for j in np.arange(0,len(grid_sfr)):
        # want to compare if streambed is lower than the layer bottom
        # 1 will be subtracted from each z value to make sure it is lower than the model top in the upper reaches
        if lay_elev[j] < (grid_sfr.z.values-1)[j]:
            sfr_lay[j] = i 
        


# %%
# KRCH, IRCH, JRCH, ISEG, IREACH, RCHLEN, STRTOP, SLOPE, STRTHICK, STRHC1, THTS, THTI, EPS, UHC

columns = ['KRCH', 'IRCH', 'JRCH', 'ISEG', 'IREACH', 'RCHLEN', 'STRTOP', 
               'SLOPE', 'STRTHICK', 'STRHC1', 'THTS', 'THTI', 'EPS', 'UHC']

sfr.reach_data.node = grid_sfr.index
sfr.reach_data.k = sfr_lay.astype(int)
sfr.reach_data.i = grid_sfr.row.values-1
sfr.reach_data.j = grid_sfr.column.values-1
sfr.reach_data.iseg = xs_sfr.iseg
sfr.reach_data.ireach = xs_sfr.reach_new
sfr.reach_data.rchlen = xs_sfr.length_m.values
sfr.reach_data.strtop = grid_sfr.z.values-1
sfr.reach_data.slope = grid_sfr.slope.values
sfr.reach_data.strthick = 2 # guess 2 meters thick streambed
# 0.004125 m/s average which is 356.4266 m/d
sfr.reach_data.strhc1 = soilKs_array[sfr.reach_data.i, sfr.reach_data.j]
sfr.reach_data.thts = soiln_array[sfr.reach_data.i, sfr.reach_data.j]
sfr.reach_data.thti = sfr.reach_data.thts
sfr.reach_data.eps = soileps_array[sfr.reach_data.i, sfr.reach_data.j]


# %%
mb4rl = pd.read_csv(sfr_dir+'michigan_bar_icalc4_data.csv', skiprows = 0, sep = ',')
# mb4rl.plot(x='gage_height_va',y='discharge_va', legend = False)
# plt.xlabel('Gage height (m)')
# plt.ylabel('Discharge $(m^3/d$)')
# plt.ticklabel_format(style='scientific') # plain to show all zeros
# plt.title('Simplified USGS Michigan Bar Rating Curve')
# plt.savefig('Plots/Model_SFR_UZF_Progress/MB_ratingcurve', dpi = 300, bbox_inches='tight')

# %% [markdown]
# ## Define segment data

# %% [markdown]
# Median grain size (mm) ranges from 1 mm – 30 mm along surveyed sites, which gives a range of 0.026-0.035 for a stable channel
# Moderate channel irregularity due to channel scouring and pools alternating, range of 0.006-0.010
# Gradual cross section change: 0.000 adjustment
# Effect of obstructions: minor due to occasional downed logs and debris in river, 0.005-0.015
# Amount of vegetation: large on banks due to willows and cottonwood trees, 0.025-0.050, and negligible in the channel
# Degree of meandering: minor due to levees, m = 1.0
#
# n = (nb+n1+n2+n3+n4)*m (b=base,1=surface irregularity, 2 = XS variation, 3 = obstructions, 4 = vegetation, m = correction for meandering)
# n = (0.03+0.08+0.01) = 0.048 in channel
# n = (0.048 +0.03) = 0.078 on banks
#

# %%
# There is one dictionary key for each stress period (starting at 0) and in each dictionary key there is a 
# rec array holding an entry for each stream segment where nseg shows which segment it is (ie no dictionary key for segment)
# If ITMP = 0 Item 4 is used, if ITMP >0 Item 6 is used, 
# if ITMP <0 the stream segment data not defined in Item 4 will be reused form the last stress period
nss = NSS
if sfr.dataset_5[0][0] > 0 :
    # For the first stress period
    t = 0
    for i in np.arange(0, nss):
        # Section 6a
        sfr.segment_data[t][i].nseg = int(i)+1
        if sfr.segment_data[t][i].nseg ==1:
            sfr.segment_data[t][i].icalc = 4
        else:
            sfr.segment_data[t][i].icalc = 2 # Mannings and 8 point channel XS is 2 with plain MF, 5 with SAFE
        if sfr.segment_data[t][i].nseg != nss:
            sfr.segment_data[t][i].outseg =sfr.segment_data[t][i].nseg +1
        elif sfr.segment_data[t][i].nseg == nss:
            sfr.segment_data[t][i].outseg = 0
        if sfr.segment_data[t][i].nseg != 1:
            sfr.segment_data[t][i].iupseg =sfr.segment_data[t][i].nseg - 1
            # Integer value that defines priority for diversion
            sfr.segment_data[t][i].iprior = -2 #diversion made will max out at flow in channel leaving no flow in channel
        elif sfr.segment_data[t][i].nseg ==1:
            sfr.segment_data[t][i].iupseg = 0
        if sfr.segment_data[t][i].icalc ==4:
            sfr.segment_data[t][i].nstrpts = len(mb4rl)
        # Defined 
        if sfr.segment_data[t][i].nseg != 1:
            sfr.segment_data[t][i].flow= 1
        elif sfr.segment_data[t][i].nseg ==1:
            sfr.segment_data[t][i].flow= 2.834*86400. # m3/day, originally 15 m3/s
        sfr.segment_data[t][i].runoff = 0.0
        sfr.segment_data[t][i].etsw = 0.01
        sfr.segment_data[t][i].pptsw = 0.01
        # Manning's n data comes from Barnes 1967 UGSS Paper 1849 and USGS 1989 report on selecting manning's n
        # RoughCH is only specified for icalc = 1 or 2
        if sfr.segment_data[t][i].icalc == 1 or sfr.segment_data[t][i].icalc ==2:
            sfr.segment_data[t][i].roughch = 0.048
        # ROUGHBK is only specified for icalc = 2
        if sfr.segment_data[t][i].icalc == 2 or sfr.segment_data[t][i].icalc == 5:
            sfr.segment_data[t][i].roughbk = 0.083 # higher due to vegetation
            
            
# Define stress period data need one for each stress period
# Dataset 5 will be built automatically from segment_data unless specified
# ITMP (int) for reusing or reading stream seg data that can change each stress period
#IRDFLG, 0 is input data printed, greater than 0 input data is not printed
# doesn't seem to change the value
# IPTFLG, 0 is streamflow-routing results printed, greater than 0 not

# %%
# Pull out data for upstream and downstream reach of each segment
up_data = xs_sfr.drop_duplicates('iseg')
dn_data = xs_sfr.sort_values('reach_new',ascending = False).drop_duplicates('iseg').sort_values('iseg')


# Need to return to later and remove hard coding
# These are getting used for initial guesses
# Read in first stress period when ICALC = 1 or 2 and ISFROPT is 5
# Dataset 6b
sfr.segment_data[0].hcond1 = sfr.reach_data.strhc1[0]
sfr.segment_data[0].thickm1 = 2
sfr.segment_data[0].elevup = up_data.z.values
sfr.segment_data[0].width1 = 20
sfr.segment_data[0].depth1 = 1
sfr.segment_data[0].thts1 = 0.4
sfr.segment_data[0].thti1 = 0.15
sfr.segment_data[0].eps1 = 4
sfr.segment_data[0].uhc1 = sfr.reach_data.strhc1[0]

# Dataset 6c
sfr.segment_data[0].hcond2 = sfr.reach_data.strhc1[-1]
sfr.segment_data[0].thickm2 = 2
sfr.segment_data[0].elevdn = dn_data.z.values
sfr.segment_data[0].width2 = 20
sfr.segment_data[0].depth2 = 1
sfr.segment_data[0].thts2 = 0.4
sfr.segment_data[0].thti2 = 0.15
sfr.segment_data[0].eps2 = 4
sfr.segment_data[0].uhc2 = sfr.reach_data.strhc1[-1]


# %%
# Change column name to float type for easier referencing in iteration
XS8pt.columns = XS8pt.columns.astype('float')
# Pre-create dictionary to be filled in loop
sfr.channel_geometry_data = {0:{j:[] for j in np.arange(2,len(XSg)+2)}  }

xsnum = 2
for k in XSg.Site.values:
        pos = int(XS8pt.columns.get_loc(k))
        XCPT = XS8pt.iloc[:,pos].values
        ZCPT = XS8pt.iloc[:,pos+1].values
        ZCPT_min = np.min(ZCPT)
        ZCPT-= ZCPT_min
        sfr.channel_geometry_data[0][xsnum] = [XCPT, ZCPT]
        xsnum += 1


# %%
FLOWTAB = mb4rl.discharge_va.values
DPTHTAB = mb4rl.gage_height_va.values
WDTHTAB = mb4rl.chan_width.values
sfr.channel_flow_data = {0: {1: [FLOWTAB, DPTHTAB, WDTHTAB]}}

# %%
# sfr.plot_path(start_seg=1, end_seg=0, plot_segment_lines=True)
# plt.savefig('Plots/Model_SFR_UZF_Progress/sfr_elev_vs_model_top.png', dpi = 600, bbox_inches='tight')

# %% [markdown]
# ## Complex ibound

# %% [markdown]
# ### Define no flow cells based on elevation, informed by DWR cross sections and geologic maps of volcanic geology fingers leaving the mountains
# In general, the location of Michigan Bar is near the boundary where there is total volcanics to majority alluvium. However there is a major finger North and South of the Cosumnes River of andesititc conglomerate, sandstone, breccia

# %%
# Simplified ibound, only no flow cell if it is below the bottom of the Mehrten Formation
# Specify no flow boundary based on rough approx of geology (upper basin volcanics)
# ibound = np.ones([nlay, nrow,ncol])
strt = np.ones((nlay, nrow, ncol), dtype = np.float32)
# The model should start in hydraulic connection
strt[:,:,:] = np.mean(m.dis.top[:,:], axis = 0)
# strt[:,:,:] = m.dis.top[:,:]

# # no_flow_bound[0,:]
# for k in np.arange(0,nlay):
#     for j in np.arange(0,ncol):
#         for i in np.arange(0,nrow):
#             # Check whether the bottom elevation is actual below where there shouldn't be flow
#             if botm[k,i,j] < no_flow_bound[i,j]:
#                 ibound[k,i,j] = 0


# %%
plt.imshow(strt[0,:,:])
plt.colorbar()
plt.show()
plt.imshow(dem_data-strt[0,:,:])
plt.colorbar()
plt.show()

# %%
ibound = np.ones([nlay, nrow,ncol])
strt = np.ones((nlay, nrow, ncol), dtype = np.float32)
# The model should start in hydraulic connection
strt[:,:,:] = np.mean(m.dis.top[:,:], axis = 0)

cutoff_elev = 56
ibound = ibound*(dem_data<cutoff_elev)

chd_locs = np.where(((dem_data<cutoff_elev)&(dem_data>cutoff_elev*0.985)))

plt.imshow(ibound[0,:,:])

# %%
# the heads at the constant head cells should be at or below the ground surface
strt[0,chd_locs[0],chd_locs[1]] = dem_data[chd_locs[0],chd_locs[1]]-10

# %%
strt_heads_sp6 = np.loadtxt('model_heads_after_sp6_use_for_strt_heads.tsv', delimiter = '\t')
# strt[:,:,:] = strt_heads_sp6
plt.imshow(dem_data- strt_heads_sp6, vmin = 0, vmax = 20)
plt.colorbar()

# %% [markdown]
# ### Create a line bounding the noflow region to set the specified head boundary

# %%
import pprint
from rasterio.features import shapes, rasterize

# The function shapes from rasterio requires uint8 format
ibound_line = ibound.astype(rasterio.uint8)
out = shapes(ibound_line,connectivity = 8)
alldata = list(out)

# maxl = 0
maxl = np.zeros(len(alldata))
for i in np.arange(0,len(alldata)):
    maxl[i] = len(alldata[i][0].get('coordinates')[0])
#     if len(alldata[i][0].get('coordinates')[0])>maxl:
#         maxl = len(alldata[i][0].get('coordinates')[0])
#         ind = i
# select the two longest linestring indexes (1st will be chunk down of divide (lower elevation) 2nd will chunk above (high elev))
maxl1, maxl2 = np.where(maxl>np.mean(maxl))[0]
print(maxl[maxl>np.mean(maxl)])

# %%
temp = alldata[maxl2][0].get('coordinates')[0]
tl = LineString(temp)
tl

# %%
from shapely.ops import LineString, linemerge, polygonize, unary_union
tl = LineString(temp)

# Get the constant head or general head boundary after the no flow cells
linerast = rasterio.features.rasterize([tl], out_shape = np.array((nrow,ncol)))
# remove far east bound line
linerast[:,ncol-1] = 0
fix_bound = np.min(np.argwhere(linerast[0,:]==1))
linerast[0,:] = 0
linerast[0,fix_bound]
np.shape(linerast)

# ibound[0,linerast==1] = -1

# %%
from shapely.ops import LineString, linemerge, polygonize, unary_union
tl = LineString(temp)
tu = unary_union(tl)
poly = list(polygonize(tu))
# Set the polygon/raster for the top layer, no buffer needed
poly0 = poly[0].buffer(distance = 0)
polyrast0 = rasterio.features.rasterize([poly0], out_shape = np.array((nrow,ncol)))
# Set the polygon/raster for the top layer, slight buffer needed to expand geologic formation outward with depth as 
# naturally occurs
poly1 = poly[0].buffer(distance = 13)
polyrast1 = rasterio.features.rasterize([poly1], out_shape = np.array((nrow,ncol)))
# Set the polygon/raster for the bottom layer, largest buffer needed
poly2 = poly[0].buffer(distance = 17)
polyrast2 = rasterio.features.rasterize([poly2], out_shape = np.array((nrow,ncol)))

ibound = np.ones([nlay, nrow,ncol])
# Need to decide whether all layers or just the top layer are affected by ibound from elevation
# it is better to define the top layer with a simple dem>elevation check than the rasterize functins that isn't perfect
# ibound[0,polyrast0==1] = 0
# Need to decide whether all layers or just the top layer are affected by ibound from elevation
ibound[-2,polyrast1==1] = 0
# Need to decide whether all layers or just the top layer are affected by ibound from elevation
ibound[-1,polyrast2==1] = 0

# The bottom boundary has a dip of 1-2 degrees which is essentially a slope of 0.015 based on given cross section data
# The layer thickness for TPROGS
laythk = tprog_thick
# It appeared shapely buffer is on the scale of kilometers
run = (laythk/0.015)/1000
run_const = run
for i in np.arange(1,nlay-2):
    # error saying poly[i] is not subscriptable
    polyi = poly[0].buffer(distance = run)
    polyrast = rasterio.features.rasterize([polyi], out_shape = np.array((nrow,ncol)))
    # Need to decide whether all layers or just the top layer are affected by ibound from elevation
    ibound[i,polyrast==1] = 0
    run += run_const

# %%
# wherever the constant head/specified head bound is the cells need to be active
ibound[0,dem_data>cutoff_elev] = 0
ibound[0,linerast==1] = 1
plt.imshow(ibound[0,:,:])
plt.colorbar()

# %%
# # copy the ibound array to alter the geology array to set these cells as low permeability formations
# either marine or volcanic based
deep_geology = ibound[:,:,:]

# reset ibound to all active cells to reduce non-linearity
# still need to take account of no flow cells for lake package
ibound = np.ones([nlay, nrow,ncol])

# %%
# plt.plot(dem_data[:,end_no_flow])

# plt.plot((dem_data[:,end_no_flow]+np.min(dis.top.array[:,end_no_flow]))/2)
# plt.plot([0,80],[np.mean(dis.top.array[:,end_no_flow])]*2)
# plt.plot([0,80],[np.min(dis.top.array[:,end_no_flow])]*2)

head = dem_data[linerast==1]
# head = dem_data[:,mean_col]
plt.plot(head)
headmin = np.min(dem_data[linerast==1])
ch_weight = 0.6
head = head*(1-ch_weight)+headmin*ch_weight
plt.plot(head)
plt.plot(headmin)


# %% [markdown]
# ### Remove no flow cells in the first layer where there are stream cells

# %%
# Get indices of SFR cells
sfr_x = grid_sfr.row.values.astype(int)-1
sfr_y = grid_sfr.column.values.astype(int)-1
# Check where SFR cells overlap the no flow cell array 
overlap = ibound[0, sfr_x, sfr_y]==0
# Convert the cells that overlap from inactive to active cells
ibound[0,sfr_x[overlap], sfr_y[overlap]] = 1
# Check where SFR cells overlap the constant head cell array 
overlap = ibound[0, sfr_x, sfr_y]==-1
# Convert the cells that overlap from inactive to active cells
ibound[0,sfr_x[overlap], sfr_y[overlap]] = 1

# %%
plt.imshow(ibound[0,:,:])
plt.colorbar()

# %%
lak_inactive = (lakarr[0,:,:]-1)*-1
# Make the cells in the first layer that underly the lake inactive
ibound[0,:,:] = lak_inactive*ibound[0,:,:]

# %%
# where the specified head boundary is the cells must be active
# ibound[0,chd_locs[0], chd_locs[1]] = -1

# %%
fig, ax = plt.subplots(figsize = (10,10))
plt.imshow(ibound[0,:,:])

cbar = plt.colorbar(shrink = 0.4, ticks = np.arange(-1,2), values = np.arange(-1,2), spacing = 'proportional' )
cbar.ax.set_yticklabels(np.array(('Constant head', 'No flow cell', 'Active Cell'), dtype = 'object' ))
plt.savefig('Plots/Model_SFR_UZF_Progress/Ibound_plot.png', dpi = 300)

# %%
# Basic package, BAS

# ibound < 0 is constant head
# ibound = 0 is inactive cell
# ibound > 0 is active cell
# strt is array of starting heads
bas = flopy.modflow.ModflowBas(model = m, ibound=ibound, strt = strt)

# %%
bas.check()


# %% [markdown]
# "When subsurface recharge (MBR2) is negligible,
# stream runoff at the mountain front (runoff measured at
# point B in Figure 1, or RO) may be considered the total
# contribution to MFR [Anderholm, 2000]." (Wilson and Guan 2004)

# %% [markdown]
# ## CHD Package Time variant head

# %%
chd = flopy.modflow.ModflowChd(model=m)

# %%

# %%
# Determine which layer the specified head cell is in
# since the if statement only checks whether the first layer is greater than the streambed elevation, 
# otherwise it would be less than and zero (most should be in layer 0)
chd_lay = np.zeros(nrow)

head = dem_data[:,ncol-1]
headmin = np.min(head)
ch_weight = 0.9
chd_vals = head*(1-ch_weight)+headmin*ch_weight



for k in np.arange(0,nlay-1):
    # pull out elevation of layer bottom
    lay_elev = botm[k, :, ncol-1]
    for i in np.arange(0,nrow):
        # want to compare if streambed is lower than the layer bottom
        # 1 will be subtracted from each z value to make sure it is lower than the model top in the upper reaches
        if lay_elev[i] > chd_vals[i]:
            chd_lay[i] = k



# %%
# layer (int), row (int), column (int), shead (float), ehead (float) shead is the head at the
# start of the stress period, and ehead is the head at the end of the stress period
# nlay_ghb = 1

# constant head boundary for mountain front recharge
# assume that near the mountains the head should be at the surface becaues the aquifer is thin

# new specified head boundary will be linear at the uppermost column to reduce nonlinearity
# as the no flow cells will be removed and replaced with low hydraulic conductivity cells
head = dem_data[chd_locs]
# chd_spd = np.zeros((len(chd_locs[0]),5))
chd_spd = np.zeros((int(np.sum((nlay-chd_lay))),5))

# # head for mountain front recharge
shead = chd_vals
ehead = chd_vals
p=0
for i in np.arange(0,nrow):
    for k in np.arange(chd_lay[i],nlay):
        chd_spd[p] = [int(k), i, ncol-1, shead[i], ehead[i]]
        p+=1
print('Number of CHD cells for upland bound', p)

# p = 0
# # head for mountain front recharge
# shead = head
# ehead = head

# for i, j in zip(chd_locs[0], chd_locs[1]):
#     chd_spd[p] = [0, i, j, shead[p], ehead[p]]
#     p+=1
# print('Number of CHD cells for upland bound', p)


# %%
plt.plot(chd_vals)

# %%
chd.stress_period_data =  {0: chd_spd}

# %%
chd.check()

# %% [markdown]
# ## Read in TPROGS data

# %%
import scipy.io
upw_dir = gwfm_dir + '/UPW_data'
tprogs_dir = upw_dir + '/TPROGS_upscale_50_mat'
# Units of K are m/d for both
khmat = scipy.io.loadmat(tprogs_dir + '\\K_xy_Upscale_50')['Upscale_50']
kvmat = scipy.io.loadmat(tprogs_dir + '\\Kz_Upscale_50')['Upscale_50']
ssmat = scipy.io.loadmat(tprogs_dir + '\\Ss_Upscale_50')['Upscale_50']
symat = scipy.io.loadmat(tprogs_dir + '\\Sy_Upscale_50')['Upscale_50']
# Column 1 is nodes, number1 is the bottom left hand corner, x,y coords are the 2nd/3rd columns, the 4th columns
# is the elevation of the center coordinate, the 5th is the water table elevation from Steve's model
# CRS is albers 3310
tcoords = pd.read_csv(tprogs_dir+'\\Center_coordinates_for_Andrew_July_2020.csv', 
                      header = None, names = ['Node', 'X','Y','elev','DTW'])
gcoords = gpd.GeoDataFrame(tcoords, geometry = gpd.points_from_xy(tcoords.X, tcoords.Y))
gcoords.crs = 'epsg:3310'
# Switch from Albers crs to UTM Zone 10 with WGS84
gcoords = gcoords.to_crs('epsg:32610')

# %% [markdown]
# ## General head boundary representing Delta/Sea Level

# %%
ghb = flopy.modflow.ModflowGhb(model = m)

# %%
# head = dem_data[linerast==1]
# find where the model boundary is actually less than 60 meters, not needed for a simple ibound
# linerast[linerast==1]= linerast[linerast==1]*(head<60)

# %%
# head = dem_data[linerast==1]
# plt.plot(head)
# headmin = np.min(dem_data[linerast==1])
# # the smaller the weight, the less impact the minimum head has
# ch_weight = 0.4
# head = head*(1-ch_weight)+headmin*ch_weight
# plt.plot(head)
# plt.plot(headmin)

# %% [markdown]
# ## Northeast GHB boundary (specified head for mountain front recharge)

# %%
# All of the rows and layers in the first column
# Head of zero for Delta at sea level, distance based on approx distance to Sac/SJ Delta
# Reduced from 5000 to 1000 meters to increase the forcing
# nlay*nrow for the number of cells on the left face of the model, 5 for the parameters associated with the GHB
# package: lay, row, col, stage, cond

# nlay_ghb = 1

# # constant head boundary for mountain front recharge
# ch_ind = np.where(linerast==1)

# ghbup_spd = np.zeros((len(ch_ind[0]),5))

# p = 0
# distance = 750

# # head for mountain front recharge
# ghbup_head = head


# for i, j in zip(ch_ind[0], ch_ind[1]):
#     cond = -np.mean(khmat)*(dis.top.array[i,j]-dis.botm.array[0,i,j])*delc/distance
#     ghbup_spd[p] = [0, i, j, ghbup_head[p], cond]
#     p+=1
# print('Number of GHB cells for upland bound', p)

# %% [markdown]
# ## Northwest and Southeast GHB boundaries based on historical WSEL

# %%
# raster cropping will be done in outside script so the only part read in will be the final array
ghb_dir = gwfm_dir+'/GHB_data'

# %%
strtyear = 2013
endyear = 2019

for year in np.arange(strtyear,endyear):
    filename = glob.glob(ghb_dir+'/final_fallWSEL_arrays/'+str(year)+'_kriged_WSEL.tsv')[0]
    df_grid = np.loadtxt(filename)

    # northwest boundary is first row
    NW = df_grid[0,:]
    # southeast boundary is last row
    SE = df_grid[-1,:]
    plt.plot(NW, label = str(year)+' NW', linestyle = '--')
    plt.plot(SE, label = str(year)+' SE')
    # 2018 looks like a pretty average year for water elevations so will use that for now
    if year ==2018:
        nwhead = df_grid[0,:]
        sehead = df_grid[-1,:]
plt.plot(dem_data[0,:],label = 'NW dem')  
plt.plot(dem_data[-1,:],label = 'SE dem')   

plt.legend(loc = [1, 0])

#  need to correct aquifer levels in foothills to be below the dem_data

# %%
# Set kriged water table elevations that are above land surface to land surface
nwhead = np.where(nwhead>dem_data[0,:], dem_data[0,:], nwhead)
sehead = np.where(sehead>dem_data[0,:], dem_data[0,:], sehead)


# %%
# Determine which layer the specified head cell is in
# since the if statement only checks whether the first layer is greater than the streambed elevation, 
# otherwise it would be less than and zero (most should be in layer 0)
nw_lay = np.zeros(ncol)
se_lay = np.zeros(ncol)


for k in np.arange(0,nlay-1):
    # pull out elevation of layer bottom
    lay_elevnw = botm[k, 0, :]
    lay_elevse = botm[k, -1, :]
    for j in np.arange(0,ncol):
        if lay_elevnw[j] > nwhead[j]:
            nw_lay[j] = k+2
        if lay_elevse[j] > sehead[j]:
            se_lay[j] = k+2
# correct from float to integer for cell referencing
nw_lay = nw_lay.astype(int)
se_lay = se_lay.astype(int)

# %%
se_lay

# %%
plt.subplots(figsize=(8,8))
for k in np.arange(0,nlay):
    plt.plot(botm[k,-1,:], label = k,alpha =0.6)

plt.plot(sehead)
plt.plot(nwhead)
plt.legend(loc=(1,.05))

# %%
# package: lay, row, col, stage, cond
# the northwest boundary is in the first row (0) and the southeast boundary is in the last row (nrow, -1)
# number of layers to go below the first active GHB cell
nlay_senw = 4
ghbnw_spd = np.zeros((np.sum(nlay - nw_lay),5))
ghbse_spd = np.zeros((np.sum(nlay - se_lay),5))
p = 0
distance = 10
for j in np.arange(0,ncol):
    
    # iterate through northwest layers first to avoid setting dry cells
    for k in np.arange(nw_lay[j], nlay):
        if k == 0:
            condnw = -np.mean(khmat)*(dis.top.array[0,j]-dis.botm.array[k,0,j])*delr/distance
        else:
            condnw = -np.mean(khmat)*(dis.botm.array[k-1, 0,j]-dis.botm.array[k,0,j])*delr/distance
        # take off 0.1 meter in head for each layer drop to support a slight vertical gradient
        ghbnw_spd[p] = [k, 0, j, nwhead[j]-k*0.1, condnw]
        
    # iterate through southeast layers to avoid setting dry cells
    for k in np.arange(se_lay[j], nlay):
        if k == 0:
            condse = -np.mean(khmat)*(dis.top.array[-1,j]-dis.botm.array[k,-1,j])*delr/distance
        else:
            condse = -np.mean(khmat)*(dis.botm.array[k-1, -1,j]-dis.botm.array[k,-1,j])*delr/distance
        # take off 0.1 meter in head for each layer drop to support a slight vertical gradient
        ghbse_spd[p] = [k, nrow-1, j, sehead[j]-k*0.1, condse]
        p+=1
print('Number of GHB cells for NW and SE bound', p*2)

# %%

# %% [markdown]
# ## Southwest GHB boundary (specified head for outflow to the Delta)

# %% [markdown]
# How much would the expected head gradient be near the delta, how fast would head decrease with depth.
# Perhaps it would only go down a few meters for every layer

# %%
# package: lay, row, col, stage, cond

# nlay_ghb = nlay
nlay_ghb = int(num_tprogs/4)

ghbdn_spd = np.zeros((nlay_ghb*nrow,5))
p = 0
distance = 1000
for i in np.arange(0,nrow):
    for k in np.arange(0,nlay_ghb):
        # temporary change of column from 1 to 0, maybe a mistake?
        if k == 0:
            cond = -np.mean(khmat)*(dis.top.array[i,0]-dis.botm.array[k,i,0])*delc/distance
        else:
            cond = -np.mean(khmat)*(dis.botm.array[k-1, i,0]-dis.botm.array[k,i,0])*delc/distance
        # take off 0.1 meter in head for each layer drop to support a slight vertical gradient
        ghbdn_spd[p] = [k, i, 0, 0-k*0.1, cond]
        p+=1
print('Number of GHB cells for Delta bound', p)

# %%
# (botm[0:nlay_ghb,:,0]<0).shape
# botm[ghbdn_spd[:,:3].astype(int)].shape

# %%
# lay, row, col for delta ghb
zxy = ghbdn_spd[:,:3].astype(int)
# (botm[zxy[:,0],zxy[:,1],zxy[:,2]]<0).shape
# drop any delta ghb cells where cell bottom is below sea level
ghbdn_spd =  ghbdn_spd[botm[zxy[:,0],zxy[:,1],zxy[:,2]]<0]



# %%
# ghb joining for delta, SE, NW boundary
# ghb_spd = np.vstack((ghbdn_spd, ghbup_spd))
ghb_spd = np.vstack((ghbdn_spd, ghbse_spd, ghbnw_spd))
# ghb_spd = np.vstack(( ghbse_spd, ghbnw_spd))

# GHB for west, north, and south model boundaries
ghb.stress_period_data =  {0: ghb_spd}
# GHB for only Delta, west side of model
# ghb.stress_period_data =  {0: ghbdn_spd}

# %%
ghb.check()

# %% [markdown]
# ## Aquifer properties with LPF

# %%
from scipy.stats import hmean


# Model homogeneous values for hydraulic properties based on averages from TPROGS
hk = np.mean(khmat)
layvka = 0 # 0 means vka is vert K, not 0 means vka is ratio of horiz to vert K
vka = hmean(hmean(hmean(kvmat)))
sy = np.mean(symat)
ss = np.mean(ssmat)

# LAYTYP MUST BE GREATER THAN ZERO WHEN IUZFOPT IS 2
laytyp = 1 # 0 is confined, >0 convertible, <0 convertible unless the THICKSTRT option is in effect
# Laywet must be 0 if laytyp is confined laywet = [1,1,1,1,1]
#ipakcb = 53 means cell-by-cell budget is saved because it is non zero (default is 53)

upw = flopy.modflow.ModflowUpw(model = m, hk =hk, layvka = layvka, vka = vka, sy=sy, ss=ss,
                               laytyp=laytyp, ipakcb=53)


# %%

def tprog2mf(arr, gcoords, lpf_arr):
    # Reshape TPROGS format so it firts in the gcoords dataframe format
    new = np.reshape(arr,[len(gcoords),arr.shape[2]])
    # Put reshaped array into dataframe so it will merge
    new = pd.DataFrame(new)
    # Convert column names to strings so they're no errors when calling them
    new.columns = new.columns.astype(str)
    # Join the TPROGS data and their spatial references
    gcoords_new = gcoords.join(new)
    
# Spatial join the grid and coordinates of TPROGS data because grid intersect is no longer working
    merged = gpd.sjoin(grid_p, gcoords_new, op = 'intersects', how = 'inner')
    
    # there are two layers on the bottom without TPROGS data
    for i in np.arange(0,nlay -2):
        index = i
        index = index.astype(str)
        lpf_arr[i, merged.row.values-1 ,merged.column.values-1] = merged[index].values
        
    # Output the array that will then be set directly in lpf. hk or vka or sy or ss
    return(lpf_arr)



# %%
1e-5*86400


# %% [markdown]
# The tuff breccia is very dense, hard and low water yielding. It is supposedly responsible for the many "haystack" hills in the eastern part of the county
#
# DWR report has a few final well pumping rates, drawdowns and specific capacities but limited.
#
# Fleckenstein et al. 2006 found the Mehrten had
# Kh = 1 to 1.8 x10^-5 m/s
# Kv = 1 to 1.8 x10^-7 m/s
# Sy = 0.15 to 0.2
# Ss = 1e-4 to 1e-3 m^-1

# %%

print(upw.hk.array.shape)
def mehrten_hp(upw_array, hp_value):
    temp = upw_array
    temp[-1,:,:] = hp_value
    return(temp)
# convert hydr. cond. from m/s to m/d
upw.hk = mehrten_hp(upw.hk.array, 1e-5*86400)
upw.vka = mehrten_hp(upw.vka.array, 1e-7*86400)
upw.sy = mehrten_hp(upw.sy.array, 0.15)
upw.ss = mehrten_hp(upw.ss.array, 1e-4)


# %%
file_name = upw_dir+'/RGM_001A_Sacramento_1981/RGM_001A_Sacramento_1981.tif'
src = rasterio.open(file_name)
band1 = src.read(1)
band1 = band1.astype('int16')
band1 = band1.astype(float)
band3 = band1.reshape(1, band1.shape[0], band1.shape[1])[:,:,:]

# turn the rasterio format into a flopy format
# for some reason flopy will not directly read in the raster .tif
sac_geology = Raster(band3, bands = (1,), crs = src.crs, transform = src.transform, 
             nodataval = None)


# %%
sac_geology.bounds

# %%
# # crop the raster to the model domain
# sac_geology.crop(vertices, invert=False)

# %%
sac_geology.plot()

# %% [markdown]
# ## Define the volcanic/marine sediment geology as very low permeability

# %%
# plt.imshow(-1*deep_geology[0,:,:]+1)
# plt.colorbar()
# convert to boolean values for referencing upw array
deep_geology = deep_geology.astype(bool)
# invert selection so previously set inactive cells will be selected to set deep geology
deep_geology = np.invert(deep_geology)

def deep_geo_hp(upw_array, hp_value, deep_geology):
    temp = upw_array
    temp[deep_geology] = hp_value
    return(temp)

# 1e-7 hydr. cond. is in the range of
# sandstone: 10^-10 to 10^-6 and
upw.hk = deep_geo_hp(upw.hk.array, 1e-7*86400, deep_geology)
upw.vka = deep_geo_hp(upw.vka.array, 1e-9*86400, deep_geology)
upw.sy = deep_geo_hp(upw.sy.array, 0.1, deep_geology)
upw.ss = deep_geo_hp(upw.ss.array, 1e-6, deep_geology)

# %%
# upw.hk = tprog2mf(khmat, gcoords, upw.hk.array)
# upw.vka = tprog2mf(kvmat, gcoords, upw.vka.array)
# upw.ss = tprog2mf(ssmat, gcoords, upw.ss.array)
# upw.sy = tprog2mf(symat, gcoords, upw.sy.array)


# %% [markdown]
# ## Well Package
#

# %% [markdown]
# each township range section is 1 square mile (2.6 square km) split into 36 sections
# so each is roughly 0.07222 km^2 each, (0.269 km)^2, which is basically the size of one model cell so there is no need to perform a spatial distribution of wells because there fluxes will add to the flux for each individual cell

# %%
wells = pd.read_csv(gwfm_dir+'/WEL_data/all_wells_type.csv')
wells_grid = gpd.GeoDataFrame(wells, geometry = gpd.points_from_xy(wells.easting,wells.northing), crs = 'epsg:32610')

# %%
# Location of wells in lat lon
# the wells are numbered from 0 to 3, 0 = domestic, 1 = irrigation, 2 = municipal, 3 = industrial
# wells = pd.read_csv(gwfm_dir+'/WEL_data/final_wcr_list_welluse.csv')
# wells

# %%
# make a geopandas dataframe with lat long
from geopandas import GeoDataFrame
from shapely.geometry import Point

geometry = [Point(xy) for xy in zip(wells.Longitude, wells.Latitude)]
df = wells.drop(['Longitude', 'Latitude'], axis=1)
wells_gdf = GeoDataFrame(df, crs="EPSG:4326", geometry=geometry)


# %%
fig, ax = plt.subplots()
# make well locations into correct coordinate system
wells_crs = wells_gdf.to_crs('epsg:32610')
wells_grid.plot(color = 'red',alpha = 0.5)
wells_crs.plot(ax=ax)
grid_p.plot(ax=ax)
#wells_crs

# %%
# Conduct spatial join using the intersect of two datasets
# ISSUE here with losing some of the wells
# wells_grid = gpd.sjoin(wells_crs, grid_p, how = "left", op= "intersects")
# wells_grid.row.isna().sum()

# %%
# wells_grid['depth'] = 0
wells_grid['flux'] = 0
wells_grid['layer'] = 0

# %%
# gpm to m^3/day, standard pumping rate is 500-1500 gpm for agricultural wells avg of 1000, maybe should be a bit lower
# that's assuming it's run 24 hours a day, they are probably usually only run for a few hours each day?
# 1000 gpm * (1 ft^3/ 7.48 gal) (ft^3/m^3)*60 min hour* 6 hours in a day

irrig_flux = 1000*(1/7.48)*(0.3048**3)*60*6

public_flux = 1500*(1/7.48)*(0.3048**3)*60*12

# averge pumping rate of domestic wells is 10-100 gpm, typically on the lower end and 
# should end up being around 50 gal/person a day and average 3 people is 150 gal/day
# 10 gpm * (1 ft^3/ 7.48 gal) (ft^3/m^3)*60 min hour* 3 hours in a day

dom_flux = 20*(1/7.48)*(0.3048**3)*60*24


# %%
# define pumping rate based on well use, average pumping rate in m^3/day
# wells_grid.loc[wells_grid.PlannedUseFormerUse == 'irrigation', 'flux'] = -irrig_flux
# wells_grid.loc[wells_grid.PlannedUseFormerUse == 'domestic', 'flux'] = -dom_flux
# wells_grid.loc[wells_grid.PlannedUseFormerUse == 'public', 'flux'] = -irrig_flux


wells_grid.loc[wells_grid.Simple_type == 'irrigation', 'flux'] = -irrig_flux
wells_grid.loc[wells_grid.Simple_type == 'domestic', 'flux'] = -dom_flux
wells_grid.loc[wells_grid.Simple_type == 'public', 'flux'] = -public_flux

# %%
if wells_grid.row.min()==1:
    wells_grid.row = (wells_grid.row-1).astype(int)
    wells_grid.column = (wells_grid.column -1).astype(int)

# %%
wells_grid.depth_m = wells_grid.depth*0.3048
# for k in np.arange(0,nlay-1):
#     # pull out elevation of layer bottom
#     lay_elev = botm[k, :, ncol-1]
#     for i in np.arange(0,nrow):
#         # want to compare if streambed is lower than the layer bottom
#         # 1 will be subtracted from each z value to make sure it is lower than the model top in the upper reaches
#         if lay_elev[i] > wells_grid.depth_m.iloc[i]:
#             wells_grid.layer.iloc[i] = k

for k in np.arange(0,nlay-1):
    # pull out elevation of layer bottom
    lay_elev = botm[k, :, :]
    for i in np.arange(0,len(wells_grid)):
        # want to compare if streambed is lower than the layer bottom
        # 1 will be subtracted from each z value to make sure it is lower than the model top in the upper reaches
        if lay_elev[wells_grid.row.values[i],wells_grid.column.values[i]] > dem_data[wells_grid.row.values[i],wells_grid.column.values[i]]-wells_grid.depth_m.iloc[i]:
            wells_grid.layer.iloc[i] = k     

# %%
# Create a list of the items needed for the stress period data from the geopandas dataset
spd_gdf = wells_grid.loc[:,['layer','row','column', 'flux']]
spd_gdf = spd_gdf.dropna()
spd_gdf

# %%
# Create a dictionary for the stress periods

spd = { j: spd_gdf.values for j in np.arange(0,nper)}
# spd = { 0: spd_gdf.values}

# %%
# Create well flopy object
wel = flopy.modflow.ModflowWel(m, stress_period_data=spd)

# %% [markdown]
# ## Unsaturated Zone

# %%
# required to do this method to insure correct formatting and then change parameters in code below, these are 
# mostly default values
uzf = flopy.modflow.ModflowUzf1(model = m,
                                nuztop=3, iuzfopt=2, irunflg=0, ietflg=1,
                                ipakcb=0, # unformatted output file of  ground-water recharge, ET, and ground-water discharge to land surface rates
                                iuzfcb2=61,# binary output of recharge and groundwater discharge
                                ntrail2=25, nsets=20,
                                surfdep=np.mean(np.abs(np.diff(m.dis.top.array, axis = 0))),
                                iuzfbnd=ibound[0,:,:], 
                                finf=0.01, #  temporary value infiltration rates
                                specifythtr = True,
                                thtr = 0.1, # temporary value
                                specifythti = False,
                                eps=3.5, #  temporary value Brooks-Corey relation of water content to hydraulic conductivity (epsilon)
                                thts = 0.35, #  temporary valuesaturated water content of the uz in units of volume of water to total volume
                                pet=5.000000E-08, #  temporary value potential ET
                                extdp=1, #  temporary valueET extinction depth(s)
                                extwc=0.1, #  temporary valueextinction water content below which ET cannot be removed from the unsaturated zone
                                unitnumber=19)

# %%
uzf_option = flopy.utils.OptionBlock(options_line = 'specifythtr etsquare', package = uzf, block = True)
uzf.options = uzf_option

# %% [markdown]
# ### For each simulation

# %%
# nuztop define which cell in column that recharge and discharge is simulated
# 1 = recharge/discharge for top model layer only, 2 = rch/dis to layer specified in iuzfbnd
# 3 = rch/dis to highest active cell in each column, top layer is iuzfbnd, const head node prevents deeper perc
uzf.nuztop = 3
# iuzfopt (int) 0 = no unsat zone routing, 1 = kvert specified by UZF in VKS, 2 = kvert taken from LPF or BCF input
uzf.iuzfopt = 2
# irunflg (int), whether gw discharge to land surface will be routed to stream (+1) segments specified by irunbnd or removed
# from the model and accounted as a loss (0)
# start with removed for simplicitly and this would most likely occur on the lower river where it is near model bounds
uzf.irunflg = 0
# ietflg (int) specifies whether ET will be simulated, 0 = no ET, other = ET
# Start with no ET to get simple UZF running
uzf.ietflg = 1
# flag for writing rch, ET, dis to land surface to separate unformatted file
uzf.iuzfcb2 = 0
# nstrail2, number of trailing waves used to define he water content profile following a decrease in the infiltration rate
# range between 10 and 20 is usually adequate
uzf.ntrail2 = 20
# nsets2 (in) number of wave sets used to simulate multiple infiltration periods, should be set to 20 for most problems
uzf.nsets = 20
# nuzgag = num of cells that will be specified for printing detailed information
# uzf.nuzgag = 0
# surfdep = avg undulation depth in the land surface within a finite-difference cell
# average difference between elevation along rows to estimate undulation
uzf.surfdep = np.mean(np.abs(np.diff(m.dis.top.array, axis = 0)))

# %%
# areal extent in which UZF will be simulated
iuzfbnd = ibound[0,:,:]
# removes -1 from const head boundary and convert to just active vs inactive cells
iuzfbnd[iuzfbnd[:,:] == -1] = 0
uzf.iuzfbnd = iuzfbnd[:,:]
# if irunflg>0 then specify irunbnd
# array of (int) to define stream segments that overland runoff from excess infil and gw discharge
# to land surface will be added

# uzf.irunbnd = iuzfbnd[:,:]



# %% [markdown]
# ### For each stress period

# %% [markdown]
# The code that originally imported the state soil map data was moved up above to be readily applicable to the LAK package and SFR package. gpd_mb and grid_uzf are created above

# %%
uzf.eps = soileps_array
# # Convert from percent into decimal
uzf.thts = soiln_array/100

# %%
# residual water content is equal to the porosity minus specific yield
uzf.thtr = uzf.thts.array - upw.sy.array[0,:,:]

# %%
plt.imshow(uzf.thts.array)
plt.colorbar()
plt.show()
plt.imshow(uzf.thtr.array)
plt.colorbar()
plt.show()
plt.imshow(upw.sy.array[0,:,:])
plt.colorbar()
plt.show()

# %%

fig, ax = plt.subplots(figsize = (10,10))
plt.imshow(soilKs_array)
plt.colorbar(label = 'Ksat m/day', shrink = 0.5)
# plt.savefig('Plots/Model_SFR_UZF_Progress/Ksat_from_soil_map.png', dpi = 300)

# %% [markdown]
# ### Crop coefficients and ETo

# %%
# crop_path = uzf_path+'\\Oversized_rectangle_2007_2019'
crop_path = uzf_path+'\\Newmodeldomain_cropdata_2007_2019'

import pathlib

crop_raster_list = list(pathlib.Path(crop_path).glob('*.tif'))

crop_dbf = gpd.read_file(crop_path+'\\CDL_2018_clip_20201020153149_1205888229.tif.vat.dbf')
crop_dbf1 = gpd.read_file(crop_path+'\\CDL_2007_clip_20201020153149_1205888229.tif.vat.dbf')


# %%
# A-B, A-C, A-D are the percent of the season for each Kc
# The dates are the dates of the growing season
Kc = pd.read_csv(uzf_path+'\\Kc\\Kc_Current.csv',skiprows = 1)
Kc = Kc.rename(columns={'Unnamed: 0' : 'Match_kc'})


# %%
def crop_raster2array(file_name, Kc):
    src = rasterio.open(file_name)
    band1 = src.read(1)
    band1 = band1.astype('int16')
    band1 = band1.astype(float)
    band3 = band1.reshape(1, band1.shape[0], band1.shape[1])[:,:,:]
    # turn the rasterio format into a flopy format
    # for some reason flopy will not directly read in the raster .tif
    croprio = Raster(band3, bands = (1,), crs = src.crs, transform = src.transform, 
                 nodataval = 255)
    # crop the raster to the model domain
    croprio.crop(vertices, invert=False)
    
    # The original crop raster has a cell size of 56 by 56 m so if there are less than 4 cells of one crop
    # then for certain they will not fill one cell and most likely have minimal impact considering there are 
    # 6300 model cells in one layer
    crop_hist = np.histogram(band3, bins = np.arange(0,257))
    crops_in_domain = crop_hist[1][:-1][crop_hist[0]>4]

    domain_dbf = crop_dbf.iloc[crops_in_domain]

    domain_dbf['CLASS_NAME'] = domain_dbf.CLASS_NAME.str.replace('Dbl Crop ','')
    domain_dbf['crop_hist'] = crop_hist[0][crops_in_domain]
    Kcmatch = pd.read_csv(uzf_path+'\\Kc\\Cosumnes_crops.csv', index_col = 0)

    # domain_dbf['crop1'] = domain_dbf.CLASS_NAME.str.split('/', expand = True)[0]
    # domain_dbf['crop2'] = domain_dbf.CLASS_NAME.str.split('/', expand = True)[1]

    domain_dbf = domain_dbf.merge(Kcmatch, on = 'CLASS_NAME')
    domain_dbf = domain_dbf.merge(Kc, left_on = 'Match_kc', right_on = 'Match_kc', how = 'left')
    return(croprio, domain_dbf)



# %%
croprio, domain_dbf = crop_raster2array(crop_raster_list[-2], Kc)

# %%
## Before resampling to the grid, I need to convert from a crop label with integer format to
## a crop coefficient format that is a float that can be weighted by the number of cells

# %%
## Potential ETo spatial interpolation from CIMIS
# daily_data = pd.read_csv(uzf_path+'\\Cosumnes_daily_multistations.csv',index_col = ['Date'], parse_dates = True)
daily_data = pd.read_csv(uzf_path+'\\Cosumnes_dailyET_precip_2010_01_01_to_2019_12_31.csv',index_col = ['Date'], parse_dates = True)

coords = pd.read_csv(uzf_path+'\\CIMIS_station_lat_long.csv', index_col = 0)
coords = gpd.GeoDataFrame(coords,geometry = gpd.points_from_xy(coords['Long'] ,coords['Lat']))
# Convert WGS Lat and long to Easting and Northing in Zone 10N
coords.crs = 'epsg:4326'
coords = coords.to_crs('epsg:32610')

# clean up data so columns are by location, units of ETo are in mm
ET = daily_data.pivot_table(index = 'Date', columns = ['Stn Id', 'Stn Name', 'CIMIS Region'], values = 'ETo (mm)')
ET.iloc[-1:-10:-1].values.shape

# clean up data so columns are by location, units of Precip are in mm
rain = daily_data.pivot_table(index = 'Date', columns = ['Stn Id', 'Stn Name', 'CIMIS Region'], values = 'Precip (mm)')
rain.iloc[-1:-10:-1].values.shape
# there are quite a few days missing rainfall at one station, could generally assume those are zero rain days

# %%
# The upper cosumnes water shed covers 2,459.98 km^2 from i5 to 
# roughly cut in half for 1250 for upper watershed


# %%
# mtft = rain.resample('A').sum().filter(like = 'Sierra Foothill')/1000
# # multiply rainfall by upper watershed area in sq. m
# mtft*1250*(1000**2)/365

# %% [markdown]
# # Need to fix missing rainfall/ET data
# temporary fix for now is to set any missing data to zero

# %%
# rain.loc[rain.index > '2018-01-01'].iloc[34:40]
# Check on missing rain data, Ryde, Holt, and Staten Island aren't installed until 2015... need to find more stations
fig, ax = plt.subplots(figsize = (8,6))
rain.groupby(rain.index.year).count().plot(ax=ax)

plt.ylabel('Number of days with CIMIS data')
# plt.savefig('Plots/Model_SFR_UZF_progress/cimis_stations_data_availability.png', dpi = 600, bbox_inches='tight')

# Same situation for ET
# ET.groupby(rain.index.year).count().plot()

# %%
# code will need to be removed and the data actually fixed
rain = rain.fillna(0)
ET = ET.fillna(0)

# %%
ncell = nrow*ncol

# Filter ET, precip data for the strt_date and end_date of the model stress periods
ET_spd = ET.loc[strt_date:end_date]
rain_spd = rain.loc[strt_date:end_date]

# Get the xy cell centers for the model grid
xy = np.append(mg.xyzcellcenters[0],mg.xyzcellcenters[1], axis = 0)
# reshape the xy data
out_xy = np.transpose(np.vstack([np.reshape(xy[0:nrow], ncell),np.reshape(xy[nrow:],ncell)]))
out_xy.shape, xy.shape
# out_xy

from scipy.interpolate import griddata
in_xy = np.transpose(np.vstack([coords.geometry.x.values, coords.geometry.y.values]))
# Final reshaped array of interpolated ET
ET_final = np.zeros((len(ET_spd), nrow, ncol))
rain_final = np.zeros((len(rain_spd), nrow, ncol))

for i in np.arange(0,len(ET_spd)):
    ET_grid = griddata(in_xy, ET_spd.iloc[i].values, xi = out_xy, method = 'cubic')
    ET_final[i,:,:] = np.reshape(ET_grid, (nrow,ncol))
    rain_grid = griddata(in_xy, rain_spd.iloc[i].values, xi = out_xy, method = 'cubic')
    rain_final[i,:,:] = np.reshape(rain_grid, (nrow,ncol))

# plt.imshow(test)
# plt.colorbar()

# plt.savefig('Plots/Model_SFR_UZF_Progress/ETo distribution in mm_day.png', dpi = 300)

# %%
# convert from unis of mm to meters
rain_final = rain_final/1000
rain_final.shape
# correct for any interpolated values that set rainfall to less than 0
rain_final[rain_final<0] = 0

# %%
i = 1
# i += 1
plt.imshow(rain_final[i,:,:])
plt.colorbar()
np.nanmean(rain_final)

# %%
np.histogram(np.isnan(rain_final))
# there are still 920,000 cells with nans
daily_domain_avg = np.mean(np.mean(rain_final, axis = 2),axis = 1)
plt.plot(daily_domain_avg)
np.sum(np.isnan(daily_domain_avg))
# there are 40 days without data, but are mostly late summer so it should be okay to count as zeros for now

# %%
# t0 = time.time()
crop_data = croprio.resample_to_grid(m.modelgrid.xcellcenters,
                                m.modelgrid.ycellcenters,
                                band=croprio.bands[0],
                                method="nearest")
# resample_time = time.time() - t0
# print("Resample time, nearest neighbor: {:.3f} sec".format(time.time() - t0))

# %%
domain_dbf.columns
# def calc_kc_dates(i):
# The year for each crop for each set of dates needs to change iteratively for each crop individually because
# some crops have dates that extend into the next year that must not change until the final date of the 
# season is reached (e.g. 2018-11-01 to 2019-09-17 must stay 2018 and 2019 until 2019-09-17 is reached)
i = 2018
dates = domain_dbf.loc[:,['Beg Month','Beg Day', 'End Month', 'End Day', 'A-B', 'A-C', 'A-D']]

# Set the pandas datetime from the start and end dates of crops
dates['A'] = pd.to_datetime({'year': i, 'month':dates['Beg Month'].values, 'day': dates['Beg Day'].values})
dates['E'] = pd.to_datetime({'year': i, 'month':dates['End Month'].values, 'day': dates['End Day'].values})
# Make correction for any end dates that are in the next year
dates.E.loc[dates.E < dates.A] = dates.E.loc[dates.E < dates.A] + pd.offsets.DateOffset(years=1)

# Get the length of the growing periods
dates['num_days'] = dates.E-dates.A

# set the end date of growing period A/ start of period B
dates['B'] = dates.A + dates.num_days*(dates['A-B']/100)
# Round the dates, as we will be one a daily time step
dates.B = pd.to_datetime(dates.B.dt.date)

# # set the end date of growing period B/ start of period C
dates['C'] = dates.B + dates.num_days*((dates['A-C']-dates['A-B'])/100)

# # set the end date of growing period C/ start of period D
dates['D'] = dates.C + dates.num_days*((dates['A-D']-dates['A-C'])/100)


# %%
# ET with crop coefficients will be the same dimension as the ETo
ETc = ET_final.copy()
# set to zero initially
ETc[:,:,:] = 0
time = 0

domain_dbf['Kc'] = 0
for y in np.arange(pd.to_datetime(strt_date).year, pd.to_datetime(end_date).year+1):
    # set start and end date for range for the year to be iterated over
    yr_strt = pd.to_datetime(str(y)+'-01-01')
    yr_end = pd.to_datetime(str(y)+'-12-31')
    if yr_strt < pd.to_datetime(strt_date):
        yr_strt = pd.to_datetime(strt_date)
    if yr_end > pd.to_datetime(end_date):
        yr_end = pd.to_datetime(end_date)
        
    for dt in pd.date_range(yr_strt, yr_end):
        # First step is to get the current Kc for each crop for the time step
        domain_dbf.Kc.loc[dt > dates.A] = domain_dbf.loc[dt > dates.A, 'Kc1']
        domain_dbf.Kc.loc[dt > dates.B] = domain_dbf.loc[dt > dates.B, 'Kc2']
        domain_dbf.Kc.loc[dt > dates.C] = domain_dbf.loc[dt > dates.C, 'Kc3']
        domain_dbf.Kc.loc[dt > dates.D] = domain_dbf.loc[dt > dates.D, 'Kc4']
        domain_dbf.Kc.loc[dt > dates.E] = domain_dbf.loc[dt > dates.E, 'Kc4']
        for i,j in zip(domain_dbf.VALUE.values, domain_dbf.Kc.values):
            ETc[time,crop_data==i] = j
        ETc[time,:,:] = ETc[time,:,:]*ET_final[time,:,:]
        time += 1


# %%
# Keys correspond to the stress periods
keys = np.arange(0,nper-1).astype(int)
# Remove any na values and replace with zeros 
rain_final[rain_final<0] = 0
# Input each array for each corresponding stress period
ETc_dict = {k: ETc[k,:,:] for k in keys}
# Convert any rain values less than zero due to interpolation to zero
rain_final[rain_final<0] = 0
rain_dict = {k: rain_final[k,:,:] for k in keys}

# %%
plt.imshow(ETc[1,:,:])
plt.colorbar()

# %% [markdown]
# ## Adjusting ETc data
# For now set any ETc value that is NA to simply 0 ETc. It is mainly the developed land that doesn't have crop coefficients anyway and they normally wouldn't transpire. Rather than try to estimate Kc for developed land, one can simply remove any infiltration from rain and any ETc, assuming that all rainfall runs off to sewer systems that discharge out of the model domain or is quickly evaporated.

# %%
# nuzf1>=0 infiltration rates are specified, nuzf1 <0 infil rates from
# previous stress period are used
# infiltration rates (L/T) at land surface
uzf.finf = rain_dict
# Potential ET (L/T)
# average ETo for CA zone 14 is 57.0 in/month, below is a list from Jan-Dec of in/month
# ETo = np.array([1.55, 2.24, 3.72, 5.10, 6.82, 7.80, 8.68, 7.75, 5.70, 4.03, 2.10, 1.55])
# ETo = ETo*(2.54/30) # convert to cm/day roughly
# print(ETo, np.mean(ETo))

uzf.pet = ETc_dict

# Keys correspond to the stress periods
keys = np.arange(0,nper-1).astype(int)
rain_dict = {k: rain_final[k,:,:] for k in keys}
# ET extinction depths (L)
# relatively deep for vines, but what about the rest of ag in the basin?
# Need to get crop maps to estimate PET, EXTDEP and soil maps for FINF
uzf.extdp = 1.5
# Extinction water content, (THTS-Sy) < extwc < THTS
uzf.extwc = 0.15

# %% [markdown]
# ## Tabfile set up for SFR
#

# %%
# For the tab files the left column is time (in model units) and the right column is flow (model units)
# Time is days, flow is cubic meters per day
import numpy as np
import pandas as pd
# USGS presents flow in cfs (cubic feet per second)
inflow = pd.read_csv(sfr_dir+'MB_daily_flow_cfs_2011_2019.csv', index_col = 'datetime', parse_dates = True)

# filter out data between the stress period dates
inflow = inflow.loc[strt_date:end_date]
# covnert flow from cubic feet per second to cubic meters per day
inflow['flow_cmd'] = inflow.flow_cfs * (86400/(3.28**3))

# # np.arange(0,len(flow_cmd))
time_flow = np.vstack((np.arange(0,len(inflow.flow_cmd)),inflow.flow_cmd))
time_flow = np.transpose(time_flow)
np.savetxt('data/MF.tab',time_flow, delimiter = '\t')

# %%
# the threshold is 23m^3/s
23*(86400)/1e6, inflow.flow_cmd.min(), inflow.flow_cmd.mean(),inflow.flow_cmd.median()

# %%
plt.plot(time_flow[:,0],time_flow[:,1])
plt.xlabel('Days since '+strt_date)
plt.ylabel('Discharge ($m^3/d)$')
plt.title('Measured Flow at Michigan Bar')
plt.ticklabel_format(style='scientific') # or plain for all zeros

# plt.savefig('Plots/Model_SFR_UZF_Progress/dailyflow_MB.png', dpi = 300, bbox_inches='tight')

# %% [markdown]
# ### Add the outside TAB file for SFR Package

# %%
# Only needs to be run one time
flopy.modflow.mfaddoutsidefile(model = m, name = 'DATA',extension = 'tab',unitnumber = 56)

# %% [markdown]
# ### Lake Bathymetry file set up

# %%
# Exactly 151 lines must be included within each lake bathymetry input file and each line must contain 1 value 
#  of lake stage (elevation), volume, and area (3 numbers per line) if the keyword “TABLEINPUT” is specified in item 1a.
# A separate file is required for each lake. 
# For oneto-denier, the levees are relatively vertical and on the very most exterior such that the lake area changes
# very litle for any change in stage
# It may be of interest to subset the lakes to have varying elevaion

# the initial stage of each lake at the beginning of the run
lak_elev = dem_data[gplak2D.row-1, gplak2D.column-1]
# minimum elevation based on resampled model dem because lake stage shouldn't be below the lowest cell
lak_elev_min = np.min(lak_elev)
# maximum stage based on data from the 2m DEM because max stage shouldn't be above the levees and
# the model dem doesn't capture the effect of the levees
lak_elev_max = lak2D.MaxElev.values[0]
stages = lak_elev_min+0.1
# (ssmn, ssmx) max and min stage of each lake for steady state solution, there is a stage range for each lake
# so double array is necessary
stage_range = [[lak_elev_min, lak_elev_max]]
print(lak2D)
# gplak2D
plt.plot(dem_data[gplak2D.row-1, gplak2D.column-1])
print('The minimum lake stage is', np.min(dem_data[gplak2D.row-1, gplak2D.column-1]), 
      'm at which there is no water ponding')



# %%
# for i in np.arange(0,151):
np.sort(lak_elev)
len(lak_elev)
lak_stage = np.append(np.sort(lak_elev), np.linspace(np.max(lak_elev),lak_elev_max, num = 151-len(lak_elev)))
plt.plot(lak_stage)

# %%
# nonzero = lakarr > 0
# bdlknc = np.zeros((nrow,ncol))
temp = np.zeros((nper, nlay, nrow,ncol))
# Calculate bed leakance based on soil maps K representative values
# temp[:,:, grid_uzf.row.values-1,grid_uzf.column.values-1] = grid_uzf.Ksat_Rep.values
temp[:,:,:] = soilKs_array
bdlknc = temp[0,:,:,:]

bdlknc.shape
plt.imshow(bdlknc[0,:,:])
plt.colorbar()
# lak.bdlknc.array.shape

# %% [markdown]
# ## Create the bathymetry file for the LAK package

# %% [markdown]
# Becase the grid resolution is so coarse, there are both levee and floodplain cells in each grid cell that lead to the nearest cell being sampled to be of a much higher elevation than expected

# %%
# lake stage (elevation), volume, and area (3 numbers per line)
# for i in np.arange(0,151):
np.sort(lak_elev)
len(lak_elev)
lak_stage = np.append(np.sort(lak_elev), np.linspace(np.max(lak_elev),lak_elev_max, num = 151-len(lak_elev)))
lak_depth = lak_stage - lak_stage[0]
lak_area = np.arange(0,len(lak_elev))*(200*200)
lak_area = np.append(lak_area, lak_area[-1]*np.ones(len(lak_stage) - len(lak_area)))
lak_volume = lak_depth*lak_area
bathtxt = np.column_stack((lak_stage, lak_volume, lak_area))
# lak_area[-1]/1e6, lak2D
# plt.plot(lak_stage)
# for i in np.arange(0,len(lak_stage)):
np.savetxt('data/MF.txt', bathtxt, delimiter = '\t')

# %%
## Need to specify flux data
# Dict of lists keyed by stress period. The list for each stress period is a list of lists,
# with each list containing the variables PRCPLK EVAPLK RNF WTHDRW [SSMN] [SSMX] from the documentation.
# flux_data = np.zeros((nrow,ncol))

flux_data = {0:{0:[0,0,0,0]}}

# %%
# filler value for bdlknc until soil map data is loaded by uzf
lak = flopy.modflow.ModflowLak(model = m, lakarr = lakarr, bdlknc = bdlknc,  stages=stages, 
                               stage_range=stage_range, flux_data = 0,tabdata= True, 
                               tab_files='MF.txt', tab_units=57)

# %%
lak.tabdata = True
lak.iunit_tab


# %%
# the lak package doesn't specify the tab file unit number when the files are written
# example:      110.0     100.0     170.0   22   Item 3:  STAGES,SSMN,SSMX,IUNITLAKTAB

lak.options = ['TABLEINPUT']
# option block is not yet available for the lake package
# lak_option = flopy.utils.OptionBlock(options_line = 'TABLEINPUT ', package = lak, block = True)


# %% [markdown]
# ### Add the outside bathymetry text file for LAK package

# %%
flopy.modflow.mfaddoutsidefile(model = m, name = 'DATA',extension = 'txt',unitnumber = 57)

# %% [markdown]
# ## Output Control

# %%
# Output control
# default unit number for heads is 51, cell by cell is 53 and drawdown is 52
# (0,0) is (stress period, time step)

# For later model runs when all the data is needed to be saved
# spd = { (j,0): ['save head', 'save budget'] for j in np.arange(0,nper)}
# # Print the budget in the list file for every 50 time steps to check
# for i in np.arange(0,730,50):
#     spd[(i,0)] = ['print budget', 'save head', 'save budget']

# spd = { (j,0): ['print budget', 'save head', 'save budget'] for j in np.arange(0,nper,int(nper/30))}
spd = { (j,0): ['print budget', 'save head', 'save budget'] for j in np.arange(0,nper,2)}

# test oc 
# spd = { (0,0): ['print budget', 'save head', 'save budget'],
#       (1,0): ['print budget', 'save head', 'save budget']}

    
oc = flopy.modflow.ModflowOc(model = m, stress_period_data = spd, compact = True)

# %% [markdown]
# ## The model fails to converge on the first time step due to issues with the upper and lower specified head/ghb boundary conditions
# One solution is to increase the maximum outer iterations so that they model is able to converge by the first time step

# %%
# pcg = flopy.modflow.ModflowPcg(model = m)
# nwt = flopy.modflow.ModflowNwt(model= m)
# thickfact: portion of cell thickness used for smoothly adjusting storage and conductance coefficients to zero (default is 1e-5)
# linmeth (linear method): 1 for GMRES and 2 for XMD (1 is default)
# iprnwt: flag whether additional info about solver convergence will be printed to the main listing file (default is 0)
# ibotav: flag whether corretion will be made to gw head relative to cell-bottom if surrounded by dry cells.
# 1 = corrections and  0 = no correction (default is 0)
# options: specify comlexity of solver. SIMPLE : default solver for linear models, MODERATE for moderately nonlinear models,
# COMPLEX for highly nonlinear models (default is COMPLEX)
# Continue: if model fails to converge during a time step it will continue to solve the next time step (default is False) 
# epsrn (XMD) is the drop tolerance for preconditioning (default is 1E-4)
# hclosexmd (XMD) head closure criteria for inner (linear) iterations (default 1e-4)
nwt = flopy.modflow.ModflowNwt(model = m, headtol=0.0001, fluxtol=500, maxiterout=200, thickfact=1e-05, 
                               linmeth=1, iprnwt=1, ibotav=0, options='COMPLEX', Continue=False,
                               maxbackiter=50, backtol=1.1, maxitinner=50, ilumethod=2, 
                               levfill=5, stoptol=1e-10, msdr=15, iacl=2, norder=1, level=5, north=7, 
                               iredsys=0, rrctols=0.0, idroptol=1, epsrn=0.0001, hclosexmd=0.0001, 
                               mxiterxmd=50, extension='nwt', unitnumber=None, filenames=None)
# GMG is more successful than pcg which is fine for steady state model
# mxiter, max outer, iiter = max inner, hclose = head change criterion for convergence, 
# rclose = residual criterion for convergence

# gmg = flopy.modflow.ModflowGmg(model = m, mxiter=50, iiter=30, hclose = 1e-5, rclose = 1e-5)

# %%
m.get_package_list()
# m.remove_package('DATA')
# m.remove_package('UZF')
# m.remove_package('LAK')
# m.remove_package('WEL')


# %%
# dis.nstp = dis.nstp.array*2.4
# sfr.plot_path(start_seg = 1, end_seg = 0)
# sfr.plot(key='ireach');
m.check()
# m.check()
# lak.check()
# upw.check()

# %% [markdown]
# ## Write the input files

# %%
# Writing the MODFLOW data files
m.write_input()


# %% [markdown]
# # Run the model

# %%
dis.get_lrc(6011)
dis.get_lrc(420)

# %%
success, buff = m.run_model()

# %%
