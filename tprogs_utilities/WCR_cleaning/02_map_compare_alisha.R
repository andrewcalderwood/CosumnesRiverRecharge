
## setwd("C:/Users/ajcalder/Documents/WCR_Steve_project")

library(readr)
library(dplyr)
library(stringr)
library(mapview)
library(mapedit)
library(sp)
library(splot)
library(rgdal)
library(readr)
## locations of wells already coded by Alisha
alisha_wcr <- read_csv("QGIS_CSV_4_15_19_Andrew.csv")

## Creat spatial data points from csv of lat long, epsg: 4326 for WGS84 lat long
alisha_sp <-SpatialPointsDataFrame(coords = data.frame(alisha_wcr$`Decimal Longitude`, alisha_wcr$`Decimal Latitude`)
                                   , data = alisha_wcr, proj4string =CRS("+init=epsg:4326"))

#BRing in well database file of all wells from DWR
w <- read_csv("wellcompletionreports.csv", col_types = cols(
  DecimalLatitude = col_double(),
  DecimalLongitude = col_double(),
  DateWorkEnded = col_datetime(format = ""),
  ReceivedDate = col_datetime(format = ""),
  TotalCompletedDepth = col_double()))
# Get relevant columns 
# 15, 16 are lat, long, 24,25,26 are township, range and section
wp <- w[c(1,2,9, 10, 15, 16, 24, 25, 26)]
# If there is no way to plot it then we can't clip it out as in the dataframe
# This still includes WCRs with Lat and Long from TRS, removes unnecessary NA that can't be mapped
wp = wp[!is.na(wp$DecimalLatitude),]
wp = wp[!is.na(wp$DecimalLongitude),]
# Prepare to be mapped
wp_sp <-SpatialPointsDataFrame(coords = data.frame(wp$DecimalLongitude, wp$DecimalLatitude)
                               , data = wp, proj4string =CRS("+init=epsg:4326"))

# Add well use to alisha's data

## Coded well log loading and preparation for mapping
coded <- read.csv("WCRLower_Cosumnes_coded2020_04_23_2020.csv")[,c("WCR_No","Legacy_Log_No", "WCRLat", "WCRLong", "Lat","Long")]
coded <- coded %>% filter(!is.na(WCRLat))
num_wells <- length(coded$WCR_No)
coded_loc <- SpatialPointsDataFrame(coords = data.frame(coded$WCRLong, coded$WCRLat),
                                    data = coded, proj4string = CRS("+init=epsg:4326"))

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

# Prepare already coded data to contrast with newly downloaded
coded <- select(coded, WCRNumber = WCR_No, LegacyLogNumber = Legacy_Log_No,
                DecimalLatitude = Lat, DecimalLongtitude = Long)
coded$WCRNumber = as.character(coded$WCRNumber)
coded$LegacyLogNumber = as.character(coded$LegacyLogNumber)

## Pull out data from Alisha again as dataframe, essentially the same but cleaner colnames
alisha_data <- data.frame(alisha_sp@data) 
alisha_data <- select(alisha_data, WCRNumber = WCR.Number, LegacyLogNumber = Legacy.Log.Number,
                      DecimalLatitude = Decimal.Latitude, DecimalLongtitude = Decimal.Longitude)

# Check whether there is any overlap
wcr_data <- anti_join(wcr_data, alisha_data, by = 'WCRNumber')
wcr_data <- anti_join(wcr_data, coded, by = 'WCRNumber')


# To include pdfs that were not digitized (e.g. destruction reports) grab
# all file names that were sorted
# Old version had WCRNumber as the name
f_og = list.files(path = 'Sorted_WCR_pdfs', pattern = '.pdf', recursive = T) %>% basename()
f_daisy = list.files(path = 'Daisy_WCR_Coding', pattern = '.pdf', recursive = T) %>% basename()

f1 = list.files(path = 'Oneto_Denier_WCRs', pattern = '.pdf', recursive = T) %>% basename()
f2 = list.files(path = 'WCR_99ish_to_OD', pattern = '.pdf', recursive = T) %>% basename()
# New version had TRS and Legacy Log or WCR number
f3 = list.files(path = 'OD_to_5', pattern = '.pdf', recursive = T) %>% basename()
f_LLN = str_split_fixed(f3, pattern = '_', n = 2)[,2]
f_LLN = str_remove(f_LLN, pattern = '.pdf') %>% data.frame()
colnames(f_LLN) = 'LegacyLogNumber'
f_LLN$LegacyLogNumber = as.character(f_LLN$LegacyLogNumber)
f_WCR = append(f1, f2) %>% str_remove(pattern = '.pdf') %>% data.frame() 
colnames(f_WCR) = 'WCRNumber'
f_WCR$WCRNumber = as.character(f_WCR$WCRNumber)

wcr_data <- anti_join(wcr_data, f_WCR, by = 'WCRNumber')
# None were removed with f_LLN because they were still included in coded above
wcr_data <- anti_join(wcr_data, f_LLN, by = 'LegacyLogNumber')

# first made 4/10/2020
# write.csv(wcr_data, 'WCR_full_domain.csv')

mapview(alisha_sp, color = 'yellow')+
  mapview(coded_loc, color = 'blue')+
  mapView(domain_rec)




###################################################################
## Legacy code as DWR has changed hyperlink system

# urls = inner_join(wp_data, hp, by = 'WCRNumber')
# keep = urls[,c('WCRNumber','LegacyLogNumber', 'MTRS','WCRLink.y')]
# urls = urls$WCRLink.y
# 
# urlnames <- paste(wp_data$WCRNumber, 'pdf', sep = '.') #Create unique names using WCR Number
# (every well has a WCR Number) for pdfs to be downloaded, argument destfile

# loop to download the pdfs of all WCR that fall inside the drawn polygon
# can't use as of March 2020 as WCR pdf storage at DWR changed
# library(tictoc)

# tic()
# for (i in 1:length(urls)) {
# if(!(is.na(urls[i])))   #some of the WCRs do not have URLs so the !(is.na(urls)) statement skips the file if there is no url
#   # download.file doesn't work for https:// connections
#   download.file(urls[i], urlnames[i]) #download the pdfs to a folder in WCR_project
# }
# toc()