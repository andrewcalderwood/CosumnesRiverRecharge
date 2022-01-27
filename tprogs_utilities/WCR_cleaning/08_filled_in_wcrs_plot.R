#
library(readr)
library(sp)
library(dismo)

final <- read_csv('final_wcr_list_welluse.csv')

final_loc <- SpatialPointsDataFrame(coords = data.frame(final$Longitude, final$Latitude),
                                    data = final, proj4string = CRS("+init=epsg:4326"))


new <- read_csv('Lower_cosumnes_WCR_2020_09_21.csv')
new <- new[!is.na(new$Lat),]

new_loc <- SpatialPointsDataFrame(coords = data.frame(new$Long, new$Lat),
                                  data = new, proj4string = CRS("+init=epsg:4326"))
  

final.layer <- list("sp.points", final_loc, zcol = "Latitude")
new.layer <- list("sp.points", new_loc, zcol = "Lat")

# basemap <- gmap(final_loc@bbox, type = "terrain")

spplot(final_loc, zcol = 'Latitude', 
       sp.layout = list("sp.points",new_loc, zcol = 'Lat'))


spplot(new_loc, zcol = 'Lat', add = TRUE)

