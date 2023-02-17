# Final run through of WCRs from SGMA data viewer

library(readr)
library(dplyr)
library(stringr)
library(mapview)
library(mapedit)
library(sp)
library(splot)
library(rgdal)
library(sf)

current_wcr <- read_csv("final_wcr_list_welluse.csv")

## Creat spatial data points from csv of lat long, epsg: 4326 for WGS84 lat long
current_sp <-SpatialPointsDataFrame(coords = data.frame(current_wcr$Longitude, current_wcr$Latitude)
                                   , data = current_wcr, proj4string =CRS("+init=epsg:4326"))


w <- read_csv("wellcompletionreports.csv", col_types = cols(
  DecimalLatitude = col_double(),
  DecimalLongitude = col_double(),
  DateWorkEnded = col_datetime(format = ""),
  ReceivedDate = col_datetime(format = ""),
  TotalCompletedDepth = col_double()))
# Get relevant columns 
# 15, 16 are lat, long, 24,25,26 are township, range and section, 11 is PlannedUseFormerUse
wp <- w[c(1,2,9, 10, 11, 15, 16, 24, 25, 26)]
# If there is no way to plot it then we can't clip it out as in the dataframe
# This still includes WCRs with Lat and Long from TRS, removes unnecessary NA that can't be mapped
wp = wp[!is.na(wp$DecimalLatitude),]
wp = wp[!is.na(wp$DecimalLongitude),]
# Prepare to be mapped
wp_sp <-SpatialPointsDataFrame(coords = data.frame(wp$DecimalLongitude, wp$DecimalLatitude)
                              , data = wp, proj4string =CRS("+init=epsg:4326"))



####################################################
## Mapping of WCRs on the lower Cosumnes
# lower = editMap(mapview(coded_loc, color = 'blue'))
# shp_lower = as(lower$finished, 'Spatial')
# shp_lower = spTransform(shp_lower, crs(coded_loc))
# writeOGR(shp_lower, layer = 'polygon', dsn = 'lower_polygon.shp', driver="ESRI Shapefile")
# polygon = readOGR(dsn ='lower_polygon/lower_polygon.shp')
domain_rec = readOGR(dsn = 'GWModelDomain_UTM10N/GWModelDomain_Rec_UTM10N.shp')
# Translate proj4string of domain to mapped wells
domain_rec = spTransform(domain_rec, proj4string(wp_sp))

# Pull out all WCRs based on any lat longs that fall in the domain shape
wcr = wp_sp[domain_rec,]
wcr_data = wcr@data

# mapview(wcr)

new_data <- anti_join(wcr_data, current_wcr, by = 'WCRNumber')

