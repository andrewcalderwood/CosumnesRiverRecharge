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
# Translatation of DEM to SFR usable data
#
# - loads HEC-RAS dem to look at auto-transects, this code is improved in the Oneto-Denier project
# - code uses Constantine cross-sections and profile to identify the location of the cross-sections spatially
# - code uses Michigan bar field measurements to create a rating curve table

# %%
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyproj # for converting proj4string
import pandas as pd
import shapely 
import shapefile
import geopandas as gpd
from flopy.utils import Raster


# %%
usr_dir = os.getcwd()
while basename(usr_dir) != 'Users':
    temp = basename(usr_dir)
    usr_dir = dirname(usr_dir)
usr_dir += '/' + temp

gwfm_dir = usr_dir+'/Box/research_cosumnes/GWFlowModel/'


# %%
mb_name = os.path.join(os.getcwd(),'polygon\polygon.shp')
mb_name

# %%

# Rivers and creeks in the larger area encompassing Cosumnes River in both South American and Cosumnes Subbasins
rivers = gpd.read_file("Sac_valley_rivers/Sac_valley_rivers.shp")
mb = gpd.read_file(mb_name)

# %%
rivers.plot()
mb.plot()
rivers = rivers.to_crs('EPSG:32610')
rivers.plot()

# %%
rivers_clip = gpd.clip(rivers, mb)

# %%
rivers_clip.plot()

# %%
rivers_clip

# %%
# Stream segments, there is a new segment each time two rivers/creeks join one another
cr_line = shapely.geometry.MultiLineString(cr.geometry.values)
dc_line = shapely.geometry.MultiLineString(dc.geometry.values)
cc_line = shapely.geometry.MultiLineString(cc.geometry.values)


# %%
import fiona
from shapely.geometry import shape, mapping
cr_line.crs = "epsg:32610"
line = cr_line
m.modelgrid.epsg

crs = "epsg:32610"

# line = line.next()

geom = shape(line)

# # length of the LineString
length = geom.length

# # creation of the resulting shapefile
schema = {'geometry': 'Point','properties': {'id': 'int'}}

with fiona.open('sfr_points/sfr_points.shp', 'w', 'ESRI Shapefile', schema, crs=crs) as output:
#     # create points every 100 meters along the line
    for i, distance in enumerate(range(0, int(length), 1)):
         point = geom.interpolate(distance)   
         output.write({'geometry':mapping(point),'properties': {'id':i}}) 

# %%
# Reaches are determined by where the model grid intersects a stream segment
sfr_p = gpd.read_file('sfr_points/sfr_points.shp')
sfr_p.plot()


# %%

# Gets the corresponding polygons that contain the points
grid_sfr = gpd.sjoin(grid_p, sfr_p, how = "inner", op= "intersects")
# Essentially returns all the points as they all exist within the grid
# sfr_grid = gpd.sjoin(sfr_p, grid_p, how = "inner", op= "intersects")

# %% [markdown]
# ## HEC RAS model XS tif view

# %%
import rasterio
# src = rasterio.open("C:\\Users\\ajcalder\\Box\\Thalweg_data\\XS.tif")
src = rasterio.open("C:\\Users\\ajcalder\\Box\\Research_Calderwood\\Cos\\Terrain\\Terrain.Cos.tif")

# %%
band1 = src.read(1)

# %%
src.meta

# %%
import matplotlib.pyplot as plt
fig,ax = plt.subplots(figsize = (12,12))
t = band1[:,:]
t[t==-9999.]= float('NaN')
plt.imshow(t)
plt.colorbar()

# %%
band1

# %%
import seaborn as sns
sns.heatmap(band1, vmin = 0)

# %% [markdown]
# ## Automating cross section calculation 

# %%
cr_ind = rivers_clip[rivers_clip.GNIS_Name == 'Cosumnes River']
cr = rivers_clip.loc[cr_ind.index,]


# %%
import time
import shapely
try:
    import flopy
    from flopy.utils import Raster
except:
    fpth = os.path.abspath(os.path.join('..', '..'))
    sys.path.append(fpth)
    import flopy
    from flopy.utils import Raster

# %%
spath = "C://Users/ajcalder/Box/Research_Calderwood/dem"

raster_name = spath+'/USGS_ten_meter_dem/regional_10m.tif'

rio10_utm = Raster.load(raster_name)

# %%
from shapely.ops import linemerge
geom = linemerge(cr.geometry.values)

geom = geom.simplify(10)

# %%
# geom = linemerge(cr.geometry.values)
geom_up = geom.parallel_offset(100,'left', resolution = 32, join_style = 2)
# right hand side offsets are returned in the reverse direction
geom_down = geom.parallel_offset(100,'right', resolution = 32, join_style = 2)
# simplifying the line increased the closeness of the total length of the lines

# %%
geom_down.length, geom_up.length, geom.length

# %%
temp = pd.DataFrame(np.linspace(1,3,3))
temp['geometry'] = geom
temp.geometry[1] = geom_up
temp.geometry[2] = geom_down
tempg = gpd.GeoDataFrame(temp)
print(tempg)
# tempg.plot(markersize = 0.01, column =0, legend = True)
# shows that the parallel offset is doing what is should

# %%
import fiona
from shapely.geometry import shape, mapping
from shapely.ops import linemerge


# how often to interpolate a point
dline = 1000
# # length of the LineString
length = int(geom.length)

pointup = np.zeros((int(length/dline)+1,2))
pointdown = np.zeros((int(length/dline)+1,2))

for i, distance in enumerate(range(0, int(length), dline)):
    pointup[i,:] = geom_up.interpolate(distance).coords[:][0]
    # Making it negative because it is flipped
    pointdown[i,:] = geom_down.interpolate(-distance).coords[:][0]



# %%
from shapely.geometry import LineString
transects = pd.DataFrame(np.zeros((len(pointup),1)), columns = ['line'])
transect = LineString([pointup[0], pointdown[0]])
# transects.loc[0]= transect
transects['geometry'] = transect
for i in np.arange(0,len(pointup)):
    transects.geometry[i] = LineString([pointup[i],pointdown[i]])
    
transects.line  = np.linspace(1,len(pointup),len(pointup))
# transects


# %%
transg = gpd.GeoDataFrame(transects)
transg = transg.drop_duplicates('geometry')
transg = transg.drop(index = np.max(transg.index))
transg = transg.drop(index = 0)
transg.plot(column = 'line')

# %%
transp = transg.geometry.buffer(5)[:]
# buffer of 1 didn't capture enough pixels, buffer of 10 captured so it was 2-3 wide, buffer of 5 is just about 1 wide

# %%
import rasterio
r10 = rasterio.open(raster_name)
dem10 = r10.read(1)


# %%
r10 = Raster.load(raster_name)
# transp.values[10]
# r10.plot()
# test = r10.crop(transp.values[10])

# %%

r10.sample_polygon(transp.values[0], band = 1, invert = False)

# %%
r10.crop(transp.values[0], invert = False)

# %%

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

ax = r10.plot(ax=ax)
plt.colorbar(ax.images[0], shrink=0.7);

# %%
cos = Raster.load("C:\\Users\\ajcalder\\Box\\Research_Calderwood\\Cos\\Terrain\\Terrain.Cos.tif")

# %%
transp3 = transg.geometry.buffer(1.5)[:]
# buffer of 1 didn't capture enough pixels, buffer of 10 captured so it was 2-3 wide, buffer of 5 is just about 1 wide

# %%
cos.sample_polygon(transp3.values[0], band = 1, invert = False)

# %%
cos.plot()

# %%
cos.sample_polygon(transp.values[0], band = 1)

# %% [markdown]
# ## Constantine 2001 cross sections

# %%
import pandas as pd
# Left value is site, right is distance downstream in kilometers
xsm = pd.read_csv('setup_data\\XS_Constantine_manual.csv', sep = ',', skiprows = 1)
xsm_refh = xsm.columns.values.astype(float)
temp = np.vstack((xsm_refh[0::2], xsm_refh[1::2]))

xsm_ref = pd.DataFrame(np.transpose(temp), columns = ['Site','Distance_m'])
xsm_ref.Distance_m *=1000

# %%
import contextily as ctx

print(ctx.providers.keys())
ctx.providers.Esri.keys()

# %%

# %%
tpath = "C:\\Users\\ajcalder\\Box\\Thalweg_data\\thalweg_profile\\thalweg_profile.shp"
import geopandas as gpd
import contextily as ctx



tg = gpd.read_file(tpath)
tg.crs = 'epsg:26910'
tg = tg.rename(columns = {'Distance (':'Distance_m'})
tg.Distance_m = np.round(tg.Distance_m, -1)
# tg.loc[tg.Distance_m <= 840
match = tg.merge(xsm_ref, on = 'Distance_m')

# Plotting
fig, ax = plt.subplots(figsize = (10,10))
# extent for ctx axes comes from minx, maxx, miny, maxy
minx, miny, maxx, maxy = tg.geometry.total_bounds
ax.axis([minx, maxx, miny, maxy])

match.plot(ax=ax)
ctx.add_basemap(ax, source = ctx.providers.Esri.WorldImagery, crs=tg.crs.to_string(), alpha = 0.6)
# ctx.add_basemap(ax, source = ctx.providers.Esri.WorldStreetMap, crs=tg.crs.to_string())

plt.xlabel('Easting (m)')
plt.ylabel('Northing (m)')

plt.savefig('Plots/Model_SFR_UZF_Progress/Constantine XS Locations.png', dpi = 300)

# %%
fig, ax = plt.subplots(figsize = (12,12))
# for i in np.arange(0,len(xsm)):
#     plt.plot(xsm.iloc[:,i], xsm.iloc[:,i+1], ax=ax)
plt.xlabel('Distance from left bank (m)')
plt.ylabel('Elevation (m)')
plt.gca().set_aspect('equal', adjustable='box')

plt.plot(xsm.iloc[:,0::2], xsm.iloc[:,1::2])
plt.savefig('Plots/Model_SFR_UZF_Progress/Constantine_original_XS.png', dpi = 300, pad_inches=0.01)

# %%
# from shapely.geometry import LineString
# lp.geometry.iloc[0] = LineString([(0,0),(1,0)])

# %%
from shapely.geometry import LineString
i = 0
# Number of cross sections
numxs = int(len(xsm.columns))
# i is the cross-section number
lp = pd.DataFrame(np.linspace(1,int(numxs/2),int(numxs/2)))
lp['geometry'] = LineString([(0,0),(0,1)])

for i in np.arange(0,numxs,2):
    # Number of points in each cross section
    numl = np.sum(pd.notna(xsm.iloc[:,i]))
    # Create empty array to fill with coordinates
    lines = np.zeros((numl,2))
    # j is the number of points in each individual cross-section
    for j in np.arange(0,numl):
        lines[j,:] = (xsm.iloc[j,i], xsm.iloc[j,i+1])
    lm = LineString(lines)
    tol = 0.0
    deltol = 0.0001
    count = 0
    lms = LineString(lm).simplify(tolerance = tol)
    while len(list(lms.coords))>10:
        lms = LineString(lm).simplify(tolerance = tol)
        tol += deltol
        count += 1
    lp.geometry.iloc[int(i/2)] = LineString(lms)

# %%
print(len(list(lm.coords)), len(list(lms.coords)))
len(list(LineString(lm).simplify(tolerance = 0.4).coords))

# %%

lpg = gpd.GeoDataFrame(lp[:])
# print(lpg)
lpg.plot(markersize = 0.01, column =0, figsize = (10,10))
# len(list(lms1.coords))
plt.xlabel('Distance from left bank (m)')
plt.ylabel('Elevation (m)')

plt.savefig('Plots/Model_SFR_UZF_Progress/Constantine_simplified_XS_plot.png', dpi = 300, pad_inches=0.1)

# %%
# Performing this kind of plotting messes with the scale of the XS so do not do if printing output
# lpg.geometry = lpg.scale(xfact=1,yfact=10)
# lpg.plot(column = 0)

# %%
# for i in np.arange(0, len(lp)):
#     print(len(list(lp.geometry.iloc[i].coords)))
# lp

# %%
xscoords = np.zeros((10,numxs))
filler = np.zeros(2)
filler[:] = np.nan
for i in np.arange(0, numxs,2):
    coordtemp = np.array(list(lp.geometry.iloc[int(i/2)].coords))
    while len(coordtemp)<10:
        coordtemp = np.vstack((coordtemp,filler))
    xscoords[:,i:i+2] = coordtemp


# %%
# xs_pd, xs_pd.loc[:,29]
xs_pd.iloc[1,30]

# %%

# %%
xs_pd = pd.DataFrame(xscoords, columns = xsm_refh)
xs_pd = xs_pd.drop(index = [0,9])
for i in np.arange(0,len(xs_pd.columns),2):
    xs_pd.iloc[:,i] -= xs_pd.iloc[0,i]
xs_pd.loc[:,8:]
xs_pd.to_csv('8pointXS.csv', index = False)


# %%
match.to_file('8pointXS_locs\\8pointXS_locs.shp')

# %% [markdown]
# ## Michigan Bar USGS Flow-Depth-Width Cross Section from manual measurement data

# %%
# units are ft, ft2, ft3/s or ft/s
mb4 = pd.read_csv('setup_data\\michigan_bar_field_measurements.tsv', skiprows = 14, sep = '\t', 
                  usecols = ['measurement_nu', 'measurement_dt','gage_height_va', 'discharge_va', 'chan_width',
                             'chan_area','chan_velocity'], 
                  parse_dates = ['measurement_dt'])
mb4 = mb4.drop(index = 0)
mb4 = mb4.dropna(axis = 0)
mb4 = mb4.loc[mb4.measurement_dt > '2000-01-01 00:00:00']
mb4.gage_height_va = mb4.gage_height_va.astype('float')
mb4.discharge_va = mb4.discharge_va.astype('float')
mb4.chan_width = mb4.chan_width.astype('float')

# mb4.plot.scatter('gage_height_va', 'discharge_va')
mb4 = mb4.sort_values(by = 'gage_height_va')
mblowest = mb4.iloc[0][['gage_height_va','discharge_va','chan_width']]

# mb4.loc[mb4.gage_height_va.diff()>0.03].plot.scatter('gage_height_va', 'discharge_va')
# Perform a rolling mean to create a more continuous data set that won't cause issues due to sharpness
mb4r = mb4.rolling(window = 20).mean()
# Remove NAs create by rolling mean
mb4r = mb4r.dropna()
# Works better to consecutively remove close points in case the first clearing gives sufficient room around some points
mb4rl = mb4r.loc[mb4.gage_height_va.diff()>0.02]
mb4rl = mb4rl.loc[mb4.gage_height_va.diff()>0.04]
mb4rl = mb4rl.loc[mb4.gage_height_va.diff()>0.06]

mb4rl = mb4rl.append(mblowest)
mb4rl = mb4rl.sort_values(by = 'gage_height_va')

mb4rl.gage_height_va = mb4rl.gage_height_va.values/3.28
mb4rl.discharge_va = mb4rl.discharge_va.values/((3.28^3)/86400)
mb4rl.chan_wdith = mb4rl.chan_width.values/3.28

# mb4rl.plot.scatter('gage_height_va', 'discharge_va')
len(mb4rl)
mb4rl.to_csv('michigan_bar_icalc4_data.csv', index = False)

# %% [markdown]
# ## Tabfile set up for SFR
#

# %%
# For the tab files the left column is time (in model units) and the right column is flow (model units)
# Time is days, flow is cubic meters per day
import numpy as np
import pandas as pd
# USGS presents flow in cfs (cubic feet per second)
inflow = pd.read_csv('USGS_MB_2017_oct.tsv', delimiter = '\t')
# inflow = pd.read_csv('USGS_MB_2018_01_01_to_2019_12_31_daily.tsv', delimiter = '\t')

inflow.columns
flow_cfs = inflow['9996_00060_00003'].values.astype('float')
flow_cmd = flow_cfs * (86400/(3.28**3))

# np.arange(0,len(flow_cmd))
time_flow = np.vstack((np.arange(0,len(flow_cmd)),flow_cmd))
time_flow = np.transpose(time_flow)
np.savetxt('data/MF.tab',time_flow, delimiter = '\t')

# %%
