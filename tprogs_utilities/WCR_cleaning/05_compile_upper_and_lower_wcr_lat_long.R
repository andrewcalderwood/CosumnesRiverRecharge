#05_compile_upper_and_lower_wcr_lat_long

library(readr)
library(dplyr)
library(stringr)
library(mapview)
library(mapedit)
library(sp)
library(splot)
library(rgdal)

## locations of wells already coded by Alisha
alisha_wcr <- read_csv("Alisha_final_north_no_duplicates.csv")

## Creat spatial data points from csv of lat long, epsg: 4326 for WGS84 lat long
alisha_sp <-SpatialPointsDataFrame(coords = data.frame(alisha_wcr$Longitude, alisha_wcr$Latitude)
                                   , data = alisha_wcr, proj4string =CRS("+init=epsg:4326"))

## Coded well log loading and preparation for mapping
coded <- read_csv("WCRLower_Cosumnes_coded_final.csv", 
                  col_types = cols(WCR_No = col_character()))
coded <- coded %>% filter(!is.na(Lat))
coded_compile <- coded[,c('WCR_No', 'Legacy_Log_No','Well_use', 'Lat','Long')]
num_wells <- length(coded$WCR_No)
coded_loc <- SpatialPointsDataFrame(coords = data.frame(coded$Long, coded$Lat),
                                    data = coded, proj4string = CRS("+init=epsg:4326"))

#BRing in well database file of all wells from DWR
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
w_sp <-SpatialPointsDataFrame(coords = data.frame(wp$DecimalLongitude, wp$DecimalLatitude)
                              , data = wp, proj4string =CRS("+init=epsg:4326"))


## Join alisha's completed WCRs with state database to fill in the well use
alisha_welluse <- left_join(alisha_wcr, wp, by = 'WCRNumber')
alisha_welluse <- alisha_welluse %>% select(c('WCRNumber','Latitude','Longitude','LegacyLogNumber',
                                              'PlannedUseFormerUse'))
# rename new data to merge properly with alisha's columns
coded_compile <- coded_compile %>% rename(WCRNumber = 'WCR_No', Latitude = "Lat", Longitude = 'Long',
                                          PlannedUseFormerUse = 'Well_use', LegacyLogNumber = 'Legacy_Log_No')
# Add the north and south wcr data together
final_welluse <- rbind(alisha_welluse, coded_compile)

# Clean up well use naming convention for Maribeth
final_welluse$PlannedUseFormerUse[final_welluse$PlannedUseFormerUse == 'Water Supply Domestic'] = 'domestic'
final_welluse$PlannedUseFormerUse[final_welluse$PlannedUseFormerUse == 'Water Supply Irrigation - Agriculture'] = 'irrigation'
final_welluse$PlannedUseFormerUse[final_welluse$PlannedUseFormerUse == 'Water Supply Irrigation - Agricultural'] = 'irrigation'
final_welluse$PlannedUseFormerUse[final_welluse$PlannedUseFormerUse == 'Water Supply Industrial'] = 'industrial'
final_welluse$PlannedUseFormerUse[final_welluse$PlannedUseFormerUse == 'Water Supply Public'] = 'public'

final_welluse$PlannedUseFormerUse <- str_to_lower(final_welluse$PlannedUseFormerUse)

unique(final_welluse$PlannedUseFormerUse)
write.csv(x = final_welluse, 'final_wcr_list_welluse.csv', row.names = FALSE)

## PLotting to check
domain_rec = readOGR(dsn = 'GWModelDomain_UTM10N/GWModelDomain_Rec_UTM10N.shp')
# Translate proj4string of domain to mapped wells
domain_rec = spTransform(domain_rec, proj4string(w_sp))

mapview(alisha_sp, color = 'yellow')+
  mapview(coded_loc, color = 'blue')+
  mapView(domain_rec)
