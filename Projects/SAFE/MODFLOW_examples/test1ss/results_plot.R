library(tidyverse)
ggl1 <- read_lines('testsfr2.gg1')
ggl2 <- read_lines('testsfr2.gg2')
ggl7 <- read_lines('testsfr2.gg7')
clean_gg <- function(df){
  site = df[1]
  header = df[2]
  header = str_remove(header, '\"DATA:')
  header = str_remove(header, '\"')
  header = str_trim(header)
  header = str_split(header, pattern = '\\s+')
  header = header[[1]]
  data = df[c(-1,-2)]
  data = str_trim(data)
  data = str_split(data, pattern = '\\s+')
  data = lapply(data, as.numeric)
  data = data.frame(matrix(unlist(data), ncol = length(header), byrow = T))
  colnames(data) = header
  return(data)
}
gg1 = clean_gg(ggl1)
gg7 = clean_gg(ggl7)

# elevup is 140, elevdn is 110
# reach 20 out of 100 is 140-0.2*(30)
plot_gg <- function(gg1){
gg1$Depth = strtop + gg1$Depth

ggplot(gg1)+
  geom_line(aes(x=Time, y = `GW-Head`, color = 'GW-Head'))+
  geom_line(aes(x=Time, y=Depth, color = 'Depth'))
}
strtop = 140-0.2*(30)
plot_gg(gg1)
strtop = 140-0.8*(30)
plot_gg(gg7)
