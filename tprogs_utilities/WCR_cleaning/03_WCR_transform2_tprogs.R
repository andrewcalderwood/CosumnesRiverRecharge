#### Transform typed up WCR to the needed EAS format
# Andrew Calderwood
# Last edit: 9/4/2019

library(raster)
library(sf)
library(dplyr)
library(sp)
library(stringr)

# Read in coded wells
# coded <- read.csv("WCRLower_Cosumnes_coded_final.csv") #397 of the coded lat long
coded <- read.csv("WCR_Cosumnes_coded_final.csv") #315 of the coded lat long

df <- coded[,c('WCR_No','Legacy_Log_No', 'from','to','Simplified_Lithology','Lat','Long')]

## Desired format
## x (easting) y (northing) z (elevation) 1 (gravel) 2 (sand) 3 (sandy mud) 4 (mud)
## 1 foot intervals for log classifications
## elevations will come from a DEM


## Convert Lat and Long to easting and northing using sf package
lat_long <- na.omit(df[,c('WCR_No', 'Lat','Long')])
cord_sf <- st_as_sf(x = lat_long, coords = c("Long", "Lat"), crs = "+proj=longlat +datum = WGS84")

## Transforming coordinate to UTM using 32610 for WGS=84 
## Sacramento Area is in Zone 10 for UTM
cord.UTM <- st_transform(cord_sf, crs=("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
## Convert coordinates into EPSG:3310 (CA Albers) from latlong
cord_3310 <- st_transform(cord_sf, crs=("=proj=utm +init=epsg:3310"))

# Take UTM coordinates and place in their individual columns
df_UTM <- cord.UTM %>% st_coordinates()
df_UTM <- data.frame(df_UTM, lat_long$WCR_No) %>% rename("WCR_No" = lat_long.WCR_No, "Easting"=X, "Northing" = Y)

fill <- inner_join(df,df_UTM, by = "WCR_No")
fill[,c("Lat","Long")]<- NULL
# Remove data that has no lithology
fill <- fill[!is.na(fill$from),]
# Factor to numeric for depth intervals
fill$from <- as.numeric(fill$from)
fill$to <- as.numeric(fill$to)


## Need to split depths to 1 ft intervals
all = list()
for(i in 1:dim(fill)[1]){
  x = fill[i,]
  if(!is.na(x$from) & !is.na(x$to) & round(x$to)-round(x$from)>0){
    from = seq(round(x$from), round(x$to)-1, by = 1)
    to = from+1
    y = data.frame(from, to)
    y[3:7] = x[c(1,2,5,6,7)]
    all[[i]] = dplyr::select(y, 'WCR_No', 'Legacy_Log_No', 'from', 'to', 
                             'Simplified_Lithology', 'Easting', 'Northing')
  }
  else{all[[i]] = x}

}

clean = do.call(rbind.data.frame, all)

# write cleaned data to csv
# write.csv(clean,'WCR_OD_99_plusAR_new_cleaned.csv',row.names = F)

# Get elevations for input file
## DEM is already in EPSG:3310 format
dem <- raster("C:/Users/ajcalder/Box/research_cosumnes/Possible 99 to 5 extension/WCR_project/Cosumnes_DEM_3310.tif")
# Extract relevant elevations from the DEM
elev <- extract(dem,cord_3310,method='simple',buffer=NULL)

# join coordinates and elevation
elevcord <- cord_3310
elevcord$elev <- elev

final <- left_join(clean, elevcord, by = "WCR_No")
# Convert from and to in units of feet to units of meters
final$from <- final$from/3.28084
final$to <- final$to/3.28084

# Convert from and to from depths to elevations
final$from <- final$elev - final$from
final$to <- final$elev - final$to

# Convert Lithology from text to numeric values
# Change from factors to strings
final$Simplified_Lithology <- as.character(final$Simplified_Lithology)

# Get the indexes of data corresponding to the strings
mud_ind <- str_detect(final$Simplified_Lithology, regex("mud", ignore_case = TRUE))
sandymud_ind <- str_detect(final$Simplified_Lithology, regex("sandy mud", ignore_case = TRUE))
sand_ind <- str_detect(final$Simplified_Lithology, regex("sand", ignore_case = TRUE))
gravel_ind <- str_detect(final$Simplified_Lithology, regex("gravel", ignore_case = TRUE))

# Set the numeric values, sand must come before sandy mud or it over rides the values
final$Numeric_Lithology <- 0
final$Numeric_Lithology[mud_ind] <- 4
final$Numeric_Lithology[sand_ind] <- 2
final$Numeric_Lithology[sandymud_ind] <- 3
final$Numeric_Lithology[gravel_ind] <- 1

plot(final$Numeric_Lithology)

# Convert UTM Zone 10 easting and northing into California Albers coordinates
ret <- as_tibble(final$geometry)
library(tidyverse)
final$Easting <- unlist(map(final$geometry,1))
final$Northing <- unlist(map(final$geometry,2))
# remove geometry column because the , causes issues when writing to a csv file, could use txt file to avoid...
final_nogeom <- select(final, -c(geometry))

write.csv(final_nogeom, 'WCR_OD_99_plusAR_new_clean_elevation.csv',row.names = F)
