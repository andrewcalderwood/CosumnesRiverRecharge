#' The following example demonstrates
#' how a search and replace string task 
#' can be performed with R across several files

# dir_nam = 'Projects/'
dir_nam = 'modflow_development/'
# dir_nam = 'projects/sasb_smc_reeval/'
# dir_nam = '/'
filenames <- list.files(dir_nam,recursive = TRUE)
# filter for file names with .R, .py, .ipynb
filenames <- filenames[stringr::str_detect(filenames, "\\.R|\\.py|\\.ipynb")]
f = filenames[1]

## Find the source R script that creates an output
output_name <- 'CIMIS'
# output_name <- 'sasb'

for( f in filenames ){
  # create filepath
  fp <- paste0(dir_nam,f)
  x <- readLines(fp)
  # replace standard here("data") file path with here(proj_dir,"data") path
  # x <- stringr::str_contains(x,output_name)
  
  lines <- (1:length(x))[stringr::str_detect(x,output_name)]
  if(length(lines)>0){
    print(paste0('Filename: ', f))
    print('Line numbers: ')
    print(paste(lines, sep=' '))
  }

}


