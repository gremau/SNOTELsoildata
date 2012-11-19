# Simple R script to count the number of met and soil datafiles available and 
# identify which sites are short of full records

metfiles <- read.csv('../rawdata/allsensors_daily/filelist.txt')
colnames(metfiles)[1] <- 'siteid'
colnames(metfiles)[2] <- 'years'


soilfiles <- read.csv('../rawdata/soilsensors_hourly/filelist.txt')
colnames(soilfiles)[1] <- 'siteid'
colnames(soilfiles)[2] <- 'years'

table(metfiles) # will show a nice table of sites, years and number of files
table(soilfiles)
