# Get all the file names of downloaded pdfs
library(stringr)

all = list.files('All_WCR_pdfs', recursive = T) %>% basename()

#Bring in relevant well files to pull out
# new = read.csv('new_coded.csv')
new = read.csv('WCR_full_domain.csv')

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


f = paste(str_sub(final, start = 1, end = 6), '/', sep = '')
f = paste('All_WCR_pdfs/', f, final, sep = '')
# Some WCRs may be duplicated
# copy relevant files to folder
file.copy(f, 'All_WCR_pdfs/All', overwrite = TRUE )

# final = data.frame(matrix(final, nrow=length(final)/2, byrow=T))

# antijoin(keep, final)