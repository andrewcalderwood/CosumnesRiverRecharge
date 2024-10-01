# helper: get sets of overlapping points

get_set <- function(x, y){zerodist(x)[, y]}



# kriging function

ordinary_krige <- function(d, grid, aoi){
  
  
  
  # fix sets of overlapping points
  
  # index of sets 1 and 2: wells with an overlapping observation
  
  s1 <- get_set(d, 1)
  
  s2 <- get_set(d, 2)
  
  
  
  # get parallel minima of overlapping points
  
  min_list_mt = pmin(d[s1, ]$wse_ft, d[s2, ]$wse_ft)
  
  
  
  # replace DGBS of set 2 wells with average of set 1 and 2
  
  d[s2, "wse_ft"] <- min_list_mt
  
  
  
  # remove set 1 wells
  
  d <- d[-s1, ]
  
  
  
  # rm extreme values
  
  d <- d[d@data$wse_ft > quantile(d@data$wse_ft, 0.025) &
           
           d@data$wse_ft < quantile(d@data$wse_ft, 0.975), ]
  
  
  
  # depth below land surface
  
  gs <- gstat(formula = wse_ft ~ 1, # spatial data -- fit xy as idp vars
              
              locations = d)      # groundwater monitoring well locations
  
  
  
  v <- variogram(gs,              # gstat object
                 
                 width = 0.5)     # lag distance
  
  plot(v)                         # check for range and patial sill and input below
  
  fve <- fit.variogram(v,         # takes `gstatVariogram` object
                       
                       vgm(700,  # partial sill: semivariance at the range
                           
                           "Lin", # linear model type
                           
                           15, # range: distance where model first flattens out
                           
                           200))  # nugget: semivariance at the y-intercept
  
  
  
  # plot variogram and fit
  
  plot(v, fve)
  
  
  
  # ordinary kriging
  
  kp <- krige(wse_ft ~ 1, d, g, model = fve)
  
  spplot(kp)
  
  
  
  # covert to raster brick and crop to area of interest
  
  ok <- brick(kp)                          # spatialgrid df -> raster brick
  
  ok <- mask(ok, aoi)                      # mask to aoi extent
  
  names(ok) <- c('prediction', 'variance') # name raster layers in brick
  
  
  
  # base R
  
  # plot(ok$prediction)
  
  # contour(ok$prediction, add=TRUE, nlevels=25)
  
  # plot(sasb, add=TRUE)
  
  # plot(d, add=TRUE, pch = 16, cex = 1, col="red")
  
  
  
  return(ok)}
