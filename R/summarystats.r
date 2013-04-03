# #summarystats.r
#
# ## A script for calculating inter-annual variability and between-site
#    variability in some variables (Tsoil or VWC means). It also calculates
#    some Tsoil-Tair offets in winter and the growing season.
#
#    The saved numbers (end comments) are for mean below-snow Tsoil, but 
#    the script can easily be run with VWC variables

# Get the summary datasets
source('getdata.r')

#==INTER-ANNUAL VARIABILITY===================================================

# How much inter-annual variability (Tsoil or VWC) is there at each site?

# Get below-snow tsoil means
ts <- soilTData$snowcovTs5mean # VWC means can be used here too
# And a list of sites
sitevec <- soilTData$siteTsoil # If using VWC, change this too

# Aggregate by site to get mean snowcov Tsoil and its stddev for each site
tsAgg <- aggregate(ts, by=list(sitevec),FUN=mean, na.rm=TRUE)
tsSdAgg <- aggregate(ts, by=list(sitevec),FUN=sd,na.rm=TRUE)
# Calculate the range in snowcov Tsoil at each site
tsMaxAgg <- aggregate(ts, by=list(sitevec),FUN=max, na.rm=TRUE)
tsMinAgg <- aggregate(ts, by=list(sitevec),FUN=min, na.rm=TRUE)
tsRange <- tsMaxAgg-tsMinAgg
tsRange[tsRange[,2]==0,2] <- NA
tsRange[!is.finite(tsRange[,2]),2] <- NA
# Calculate the mean of the means (+sd) and the mean of the sd (by site)
sapply(tsAgg, mean, na.rm=TRUE) # 0.2867, 0.7066,  1.3298  agg. means
sapply(tsAgg, sd, na.rm=TRUE) #0.6459, 0.7025, 0.7201  SD of the agg. means
sapply(tsSdAgg, mean, na.rm=TRUE) # 0.4212 0.4293 0.4006 mean of agg. SD
# Calculate the mean, max, and min of the interannual ranges
mean(tsRange[,2], na.rm=TRUE) # 1.0779 1.0858 1.0290
max(tsRange[,2], na.rm=TRUE) # 7.5136 7.1708 7.424
min(tsRange[,2], na.rm=TRUE) # 0.0197 0.0528 0.0666

## Inter-annual variability in Tair========== 
tair <- climData.sub$scovmat
sitevec2 <- climData.sub$siteClim
tairAgg <- aggregate(tair, by=list(sitevec2),FUN=mean, na.rm=TRUE)
tairSdAgg <- aggregate(tair, by=list(sitevec2),FUN=sd,na.rm=TRUE)
# Calculate the mean of the means (+sd) and the mean of the sd (by site)
sapply(tairAgg, mean, na.rm=TRUE) # -1.774 snowcovered and 11.2666 snow-free
sapply(tairAgg, sd, na.rm=TRUE) # 1.697 and 1.42334
sapply(tairSdAgg, mean, na.rm=TRUE) # 0.717 and 0.8490


## GRADIENTS with elevation or SWE? ==========
# Is there an elevation gradient in snowcov tsoil standard deviation?
elevAgg <- aggregate(climData.sub$elev, by=list(sitevec),FUN=mean, na.rm=TRUE)
plot(elevAgg[,2], tsSdAgg[,2])
summary(lm(tsSdAgg[,2] ~ elevAgg[,2])) # No
# What about in snowcov Tsoil range?        
plot(elevAgg[,2], tsRange[,2])
summary(lm(tsRange[,2] ~ elevAgg[,2])) # No

# Is there a peak SWE gradient in snowcov tsoil standard deviation?
sweAgg <- aggregate(climData.sub$maxswe, by=list(sitevec),FUN=mean, na.rm=TRUE)
plot(sweAgg[,2], tsSdAgg[,2])
summary(lm(tsSdAgg[,2] ~ sweAgg[,2])) # Yes - high swe = low SD
# What about in snowcov Tsoil range?        
plot(sweAgg[,2], tsRange[,2])
summary(lm(tsRange[,2] ~ sweAgg[,2])) # Yes - high swe = low range

#==BETWEEN-SITE VARIABILITY====================================================

# How much variability is there (Tsoil or VWC) among the sites in each year.

# Get below-snow tsoil means
ts <- soilTData$snowcovTs50mean # VWC means can  be used here too
# And a list of years
yearvec <- soilTData$yearTsoil  # If using VWC, change this too

# First, this is mean snowcov Ts in all years - not sure how useful it is.
mean(ts, na.rm=TRUE) # 0.3396 0.7644 1.390422 - means

# Aggregate by year to get mean snowcov Tsoil and its stddev for each year
# (all sites)
tsAgg <- aggregate(ts, by=list(yearvec),FUN=mean,na.rm=TRUE)
tsSdAgg <- aggregate(ts, by=list(yearvec),FUN=sd,na.rm=TRUE)
# Calculate the range in snowcov Tsoil at each year
tsMaxAgg <- aggregate(ts, by=list(yearvec),FUN=max, na.rm=TRUE)
tsMinAgg <- aggregate(ts, by=list(yearvec),FUN=min, na.rm=TRUE)
tsRange <- tsMaxAgg-tsMinAgg
tsRange[tsRange[,2]==0,2] <- NA
tsRange[!is.finite(tsRange[,2]),2] <- NA
# Calculate the mean of the means (+sd) and the mean of the sd (by year)
sapply(tsAgg, mean, na.rm=TRUE) # 0.7044 0.9079 1.8725
sapply(tsAgg, sd, na.rm=TRUE) # 0.4861 0.6885 0.6648  SD of the agg. means
sapply(tsSdAgg, mean, na.rm=TRUE) # 0.6956 0.8709 0.7802  mean of agg. SD
# Calculate the mean, max, and min of the yearly ranges
mean(tsRange[,2], na.rm=TRUE) # 4.458 4.8852  4.8155
max(tsRange[,2], na.rm=TRUE) # 8.4155 11.0697 7.8216
min(tsRange[,2], na.rm=TRUE) # 0.62885 0.4799 0.452

## GRADIENTS with Tair or SWE? =============
# Is there a mean Tair gradient in snowcov tsoil standard deviation?
scovmatAgg <- aggregate(climData.sub$scovmat,by=list(yearvec),FUN=mean,
                        na.rm=TRUE)
plot(scovmatAgg[,2], tsSdAgg[,2])
summary(lm(tsSdAgg[,2] ~ scovmatAgg[,2])) # No
# What about in snowcov Tsoil range?        
plot(scovmatAgg[,2], tsRange[,2])
summary(lm(tsRange[,2] ~ scovmatAgg[,2])) # No

# Is there a peak SWE gradient in snowcov tsoil standard deviation?
sweAgg <- aggregate(climData.sub$maxswe, by=list(yearvec),FUN=mean, na.rm=TRUE)
plot(sweAgg[,2], tsSdAgg[,2])
summary(lm(tsSdAgg[,2] ~ sweAgg[,2])) # No
# What about in snowcov Tsoil range?        
plot(sweAgg[,2], tsRange[,2])
summary(lm(tsRange[,2] ~ sweAgg[,2])) # No


#==WINTER and GS OFFSETS=====================================================

# Get offsets and means for the snow-covered period
offset <- ts - climData.sub$scovmat
offsetAgg <- aggregate(offset, by=list(sitevec), FUN=mean, na.rm=TRUE)
sapply(offsetAgg, mean, na.rm=TRUE) # 2.0603, 2.4505, 3.1041
sapply(offsetAgg, sd, na.rm=TRUE) # 1.4064, 1.3481, 1.2677

# Get snow-free tsoil means
ts <- soilTData$snowfreeTs50mean
sitevec <- soilTData$siteTsoil
tsAgg <- aggregate(ts, by=list(sitevec),FUN=mean, na.rm=TRUE) 
mean(ts, na.rm=TRUE) # 10.4557
sapply(tsAgg, mean, na.rm=TRUE) # 10.3257, 9.7863, 8.8871 
sapply(tsAgg, sd, na.rm=TRUE) # 2.0646, 1.9670, 1.8493

# Get offsets and means for the snow-free period
offset <- ts - climData.sub$freemat
offsetAgg <- aggregate(offset, by=list(sitevec), FUN=mean, na.rm=TRUE)
sapply(offsetAgg, mean, na.rm=TRUE) # -0.9691, -1.5094, -2.4066 
sapply(offsetAgg, sd, na.rm=TRUE) # 1.8950, 1.8813, 1.8772
