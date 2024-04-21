
library(readr)
library(dplyr)
library(stringr)
library(mapview)
library(mapedit)
library(sp)
library(splot)
library(rgdal)


## locations of wells already coded by Alisha
alisha_wcr <- read_csv("QGIS_CSV_4_15_19_Andrew.csv")

## Creat spatial data points from csv of lat long, epsg: 4326 for WGS84 lat long
alisha_sp <-SpatialPointsDataFrame(coords = data.frame(alisha_wcr$`Decimal Longitude`, alisha_wcr$`Decimal Latitude`)
                                   , data = alisha_wcr, proj4string =CRS("+init=epsg:4326"))

coded <- read.csv("WCR_Alisha_OD_coded.csv")
df <- coded[,c('WCR_No','Legacy_Log_No', 'from','to','Simplified_Lithology','Lat','Long')]

## Convert Lat and Long to easting and northing using sf package
library(sf)
lat_long <- na.omit(df[,c('WCR_No', 'Lat','Long')])
cord_sf <- st_as_sf(x = lat_long, coords = c("Long", "Lat"), crs = "+proj=longlat +datum = WGS84")

## Transforming coordinate to UTM using 32610 for WGS=84 
## Sacramento Area is in Zone 10 for UTM
cord.UTM <- st_transform(cord_sf, crs=("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
## Convert coordinates into EPSG:3310 from latlong
cord_3310 <- st_transform(cord_sf, crs=("=proj=utm +init=epsg:3310"))