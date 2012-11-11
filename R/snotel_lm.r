datapath = '../processed_data/'
# Load a few datafiles
climData <- read.csv(paste(datapath, 'wyear_climatesummary.txt', sep=''))
soilTData <- read.csv(paste(datapath, 'wyear_soiltempsummary_hourly.txt',
			    sep=''))
soilVWCData <- read.csv(paste(datapath, 'wyear_soilwatersummary_hourly.txt',
			      sep=''))

# Get subset of climData that matches the soil data using merge
climrows <- data.frame(siteID=climData$siteClim,year=climData$year)
climrows$included_clim <- TRUE
soilrows <- data.frame(siteID=soilTData$siteTsoil,year=soilTData$year)
soilrows$included_soil <- TRUE
finder <- merge(climrows, soilrows, all=TRUE)
remove <- !is.na(finder$included_soil)
climData.sub <- subset(climData, remove)

reg1 <- lm(climData$totaldaysSC ~ climData$maat)
plot(climData$maat, climData$totaldaysSC)
abline(reg1)
# totaldaysSC and maat are highly correlated

reg2 <- lm(soilTData$mast20cm ~ climData.sub$totaldaysSC)
plot(climData.sub$totaldaysSC, soilTData$mast20cm)
abline(reg2)

diff <- (soilTData$mast20cm - climData.sub$maat)
reg3 <-lm (diff ~ climData.sub$totaldaysSC)
plot(climData.sub$totaldaysSC, diff)
abline(reg3)
# This looks good, totaldaysSC is highly significant

reg4 <-lm (diff ~ climData.sub$maat)
plot(climData.sub$maat, diff)
abline(reg4)
# This looks good, maat is highly significant

reg5 <-lm (diff ~ climData.sub$totaldaysSC + climData.sub$maat)
# This looks good, both terms are highly significant
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(reg5)
# However, maat and totalSCdays are highly correlated (see reg1). now what?

reg6 <-lm (diff ~ climData.sub$totaldaysSC + climData.sub$maat + 
	   climData.sub$totaldaysSC*climData.sub$maat)
# Interaction term doesn't look as significant
