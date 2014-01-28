# Load the needed climate, soil temperature, and soil vwc data for doing linear
# models
#

datapath = '../processed_data/'

# Load a few datafiles
climData <- read.csv(paste(datapath, 'wyear_climatesummary.txt', sep=''))
soilTData <- read.csv(paste(datapath, 'wyear_soiltempsummary_hourly.txt', 
			    sep=''))
soilVWCData <- read.csv(paste(datapath,
			      'wyear_soilwatersummary_hourly_smnorm.txt',
			      sep=''))

# Get subset of climData that matches the soil data using merge
climrows <- data.frame(siteID=climData$siteClim,year=climData$year)
climrows$included_clim <- TRUE
soilrows <- data.frame(siteID=soilVWCData$siteVWC,year=soilVWCData$yearVWC)
soilrows$included_soil <- TRUE
finder <- merge(climrows, soilrows, all=TRUE)
remove <- !is.na(finder$included_soil)
climData.sub <- subset(climData, remove)

# First create jasMAT - the sum of jul, aug, sept Air T
climData.sub$jasMAT <- (climData.sub$julTairmean+climData.sub$augTairmean
                        +climData.sub$sepTairmean)
# And JFM MAT - the sum of jan, feb and march MAT
climData.sub$jfmMAT <- (climData.sub$janTairmean+climData.sub$febTairmean
                        +climData.sub$marTairmean)

rm(finder, remove, soilrows, climrows, datapath)
