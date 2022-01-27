## Code converts the text file with the data for all the well logs used in the Tprogs model to individual text files
## Created: 1/16/2018
## By: Rich Pauloo with Andrew Calderwood

## setwd("C:/Users/ajcalder/Documents/WCR_Steve_project")

library(dplyr)
# Uncorrected xy coordinates
#Tsim <- read.delim("Tsim_s_input_rot_noheader.eas.txt", header = F,  sep = "\t")

Tsim <- read.delim("Tsim_corrected_xy_011719_noheader.txt", header = F, sep = "\t")
col <- dim(Tsim)[1]; row <- dim(Tsim)[2]

Tsim$EN <- paste(Tsim$V1, Tsim$V2, sep = "_") ## Creating a single identifier for each well based on Easting and Northing to use for split

splitTsim <- split(Tsim, Tsim$EN) ## make sure to use only one identifier when using split

Tsim %>% group_by(EN) %>% summarise(count = n()) ## Check the number of rows in minTsim
r1 <- Tsim %>% count(EN) #same as group_by %>% summarise(count = n())
r2 <- sapply(splitTsim, nrow) ## Check the number of rows in Test to compare that results have equal number of rows
# plot(r1$n, r2) #Check to see that the number of rows adds up (1:1)

#adds a computer readable number at the beginning of each file
names(splitTsim) <- paste(formatC(x=1:length(r2),digits = 4,flag = "0"), names(splitTsim), sep = "_") 

         # Removes the 9th row which contained the unique identifier, EN
         for (i in 1:length(splitTsim)){
splitTsim[[i]][9] <- NULL
}

# Creates individuatl text files for each well log
# row.names = F and col.names = F remove the unncecessary data
#sapply(names(splitTsim), 
#       function (x) write.table(splitTsim[[x]], file=paste(x, "txt", sep="."), row.names = F, col.names = F )   ) 

# Pulling out the given label from each xy log from the names of the dataframes in the datatable
new_xy_label <- names(splitTsim)
label <- substr(new_xy_label, start = 1, stop = 5)

########################################################################
## Pulling out x and y coordinates from data.frames in the datatable ##
x = 1; y = 1; z1 = 1; z2 = 1
for (i in 1:length(r2)) {
x[i] <- splitTsim[[i]][1,1]
y[i] <- splitTsim[[i]][1,2]
z1[i] <- splitTsim[[i]][1,3]
z2[i] <- splitTsim[[i]][length(splitTsim[[i]][,3]),3]
}

xyz = data.frame(x, y, z1, z2)
xyz$depth = xyz$z1-xyz$z2
map_xy <- data.frame(label, x, y)
#write.csv(map_xy, "xy_coord_label.csv")

