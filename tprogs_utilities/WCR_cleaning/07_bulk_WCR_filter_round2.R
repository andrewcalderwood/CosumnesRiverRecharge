# Get all the file names of downloaded pdfs
library(stringr)

all = list.files('All_WCR_pdfs', recursive = T) %>% basename()

#Bring in relevant well files to pull out
new = new_data
t = list()
# attempt to join file names and corresponding name
for( i in 1:dim(new)[1]){
  t[[i]] = c(all[str_detect(all, paste(as.character(new$LegacyLogNumber[i]),'.pdf', sep = ''))], as.character(new$LegacyLogNumber[i]))
  
  if( is.na(new$LegacyLogNumber[i]) ){
    t[[i]] = c(all[str_detect(all, paste(as.character(new$WCRNumber[i]), '.pdf', sep = ''))], as.character(new$WCRNumber[i]))
  }
  
}
# just keeps the file names
for( i in 1:dim(new)[1]){
  t[[i]] = all[str_detect(all, paste(as.character(new$LegacyLogNumber[i]), '.pdf', sep = ''))]
  
  if( is.na(new$LegacyLogNumber[i]) ){
    t[[i]] = all[str_detect(all, paste(as.character(new$WCRNumber[i]), '.pdf', sep = ''))]
  }
  
}
final = unlist(t)[!is.na(unlist(t))]
pdf_available <- final[!is.na(str_match(final, ".pdf"))]
available <- data.frame(pdf_available)
colnames(available) = 'pdf'
available$pdf = as.character(available$pdf)

pdf_done = list.files(path = 'Daisy_WCR_Coding', pattern = '.pdf', recursive = T) %>% basename()
done = data.frame(pdf_done)
colnames(done) = 'pdf'
done$pdf = as.character(done$pdf)
newpdfs = anti_join(available, done, by = 'pdf')

newpdfs$LLN = 0
for (i in 1:length(newpdfs[,1])){
  newpdfs$LLN[i] = str_sub(newpdfs[i,],start = 10, end = -5)
}

newpdfs = newpdfs %>% rename(LegacyLogNumber = 'LLN')
pdf_plot = inner_join(newpdfs, new, by = 'LegacyLogNumber')
pdf_sp <-SpatialPointsDataFrame(coords = data.frame(pdf_plot$DecimalLongitude, pdf_plot$DecimalLatitude)
                                , data = pdf_plot, proj4string =CRS("+init=epsg:4326"))

pdf_plot %>% distinct()
library(mapedit)
library(sp)
library(rgdal)
library(sf)

upper<- editMap(mapView(pdf_sp)+
  mapView(current_sp, color = 'yellow'))

up = spTransform(as(upper$finished, 'Spatial'), pdf_sp@proj4string)
pdf_up = pdf_sp[up,]
pdf_upper = pdf_up@data

pdf_lower = anti_join(pdf_plot, pdf_upper, by = 'WCRNumber')
low_sp <-SpatialPointsDataFrame(coords = data.frame(pdf_lower$DecimalLongitude, pdf_lower$DecimalLatitude)
                                , data = pdf_lower, proj4string =CRS("+init=epsg:4326"))


mapView(pdf_up, color = 'red')+
  mapView(low_sp,color = 'Yellow')+
  mapView(current_sp, color = 'blue')

write.csv(pdf_upper, 'Upper_cosumnes_WCR_2020_09_01.csv', row.names = F)
write.csv(pdf_lower, 'Lower_cosumnes_WCR_2020_09_01.csv',row.names=F)

# Transfer files for the upper Cosumnes
f = pdf_upper$pdf
f = paste(str_sub(f, start = 1, end = 6), f, sep = '/')
f = paste('All_WCR_pdfs/', f, sep = '')
# Some WCRs may be duplicated
# copy relevant files to folder
file.copy(f, 'All_WCR_pdfs/Allupper', overwrite = FALSE )

## Repeat for the lower Cosumnes
f = pdf_lower$pdf
f = paste(str_sub(f, start = 1, end = 6), f, sep = '/')
f = paste('All_WCR_pdfs/', f, sep = '')
# Some WCRs may be duplicated
# copy relevant files to folder
file.copy(f, 'All_WCR_pdfs/Alllower', overwrite = FALSE )

# final = data.frame(matrix(final, nrow=length(final)/2, byrow=T))

# antijoin(keep, final)