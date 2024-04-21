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

# %% [markdown]
# ## External code for model set up
# This notebook will hold code used in the initial set up of the model but is not needed for every run (ie transforming the large scale 10 meter dem from epsg 4326 to epsg 32610 only needs to be done once)

# %%
spath = "C://Users/ajcalder/Box/Research_Calderwood/dem"

# 1 meter dem
raster_name = spath+"/model_dem.tif"

rio = Raster.load(raster_name)

# %%
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

ax = rio.plot(ax=ax)
plt.colorbar(ax.images[0], shrink=0.7);

# %%
# 10 meter dem
raster_name = spath+"/USGS_ten_meter_dem/USGS_13_n39w122_10meterdem.tif"

rio10 = Raster.load(raster_name)

# %%
# Convert 10 meter dem crs from lat long to utm zone 10n
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling

dst_crs = 'EPSG:32610'

raster_name = spath+"/USGS_ten_meter_dem/USGS_13_n39w122_10meterdem.tif"
with rasterio.open(raster_name) as src:
    transform, width, height = calculate_default_transform(
        src.crs, dst_crs, src.width, src.height, *src.bounds)
    kwargs = src.meta.copy()
    kwargs.update({
        'crs': dst_crs,
        'transform': transform,
        'width': width,
        'height': height
    })

    with rasterio.open(spath+'/USGS_ten_meter_dem/transformed.tif', 'w', **kwargs) as dst:
        for i in range(1, src.count + 1):
            reproject(
                source=rasterio.band(src, i),
                destination=rasterio.band(dst, i),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest)

# %% [markdown]
# # Extras

# %% [markdown]
# # Model vertices, shapefile

# %%
# Get vertexes of model domain
ll = mg.get_coords(0, 0) #lower left
lr = mg.get_coords(nrow*delr, 0) #lower right
ur = mg.get_coords(nrow*delr, ncol*delc) #upper right
ul = mg.get_coords(0, ncol*delc) #upper left
print(ll, lr, ur, ul)

# Shapefile of model bounds
vertices = np.stack(np.asarray((ll,lr, ur, ul)))
vertices

# %%
geoms = Polygon(vertices)
geoms.plot() # this feature requires descartes
geoms.type

# %% [markdown]
# # Saving a polygon to a shapefile

# %%
# How to save a polygon to shapefile
import shapely
import shapefile
w = shapefile.Writer('polygon')
w.field('name', 'C')
w.poly([vertices])
w.record('polygon1')
w.close()

# %% [markdown]
# # Raster cropping

# %%
t0 = time.time()
rio10_utm.crop(vertices, invert=False)
crop_time = time.time() - t0

# %%
rio10_utm.plot()

# %%
rio10_utm.write()

# %% [markdown]
# # Capture cross section of deeper geology
# Not the best method to create the cells, its better to just plot the change points and then apply those to a grid

# %%
# vertical resolution is -400 ft to 500 ft divided by 5m (16.4 ft) 55 layers or 500/16.4 is 30 layers
# horizontal resolution is 0-38 miles, 0 to 200,640 ft
# assuming 100 meter grid would be 328 feet and need 611 rows
gelnlay = 55
gelnrow = 611
gel = np.ones((gelnlay,gelnrow))

# Pre-cretaceous, Ione and Valley Springs will be set to inactive due to low water bearing
# Slope of bottom of Mehrten is run= 32.25 mi- 26.125 and rise =(200ft--300ft)
mslp = (200-(-300))/((32.25-26.125)*5280)
# assuming the bottom of the model is the last row of the array

xlb = int((6/38)*gelnrow)
ylb = int((50/900)*gelnlay)
xlt = int((20/38)*gelnrow)
ylt = int((150/900)*gelnlay)
aslp =  (ylt-ylb)/(xlt-xlb)
ratio = int(1/aslp) # for every 1 up three overs
j = xlb
gel[-1:-ylb:-1,:] = 2
for i in np.arange(-ylb,-ylt,-1):
    gel[i,j:] = 2 # represents mehrten unit
    j = j+ratio
    
xmb = int((26/38)*gelnrow)
ymb = int((100/900)*gelnlay)
xmt = int((32/38)*gelnrow)
ymt = int((600/900)*gelnlay)
aslp =  (ymt-ymb)/(xmt-xmb)
ratio = int(1/aslp) # for every 1 up three overs
j = xmb
for i in np.arange(-1,-ymt,-1):
    gel[i,j:] = 0 # represents deep geology
    j = j+ratio
