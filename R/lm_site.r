source('getdata.r')


# Make some dataframes to analyze with linear models
reg.dat <- cbind(climData.sub[,c('siteClim','yearClim','maxswe','novSWEmean',
				 'decSWEmean','janSWEmean','febSWEmean','meltdoy','onsetdoy',
				 'totaldaysSC','maat','freemat','scovmat',
				 'jfmMAT','jasMAT','preonsetTair','JASprecip')],
		 soilTData[,c('jfmTs5mean','jfmTs20mean','jfmTs50mean',
			      'snowcovTs5mean','snowcovTs20mean',
			      'snowcovTs50mean','preonsetTs5',
			      'preonsetTs20','preonsetTs50')],
		 soilVWCData[,c('jasVWC5mean','jasVWC20mean','jasVWC50mean',
				'jfmVWC5mean','jfmVWC20mean','jfmVWC50mean',
				'preonsetVWC5','preonsetVWC20',
				'preonsetVWC50')])
# List of sites
sites <- unique(reg.dat$siteClim)

# BELOW-SNOW Tsoil - Create a data frame to store regression values
scovTs.lm <- data.frame(site=sites,
			       nyrs=as.vector(table(reg.dat$siteClim)),
			       x1Beta=0,x1Pval=0,x2Beta=0,x2Pval=0,
			       x3Beta=0,x3Pval=0,x4Beta=0,x4Pval=0,
			       x5Beta=0,x5Pval=0,x6Beta=0,x6Pval=0,
			       x7Beta=0,x7Pval=0,x8Beta=0,x8Pval=0,
                   x9Beta=0,x9Pval=0)
for (i in 1:length(sites)) {
	# Subset for the site
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	# Assign the independent variable
	y <- tmp$snowcovTs5mean
	# Create a list of independent variables
	# peak SWE, onsetdoy, preonset Tair, nov & dec SWE
	xlist <- data.frame(tmp$maxswe,tmp$onsetdoy,tmp$preonsetTair,
			    tmp$scovmat,tmp$meltdoy,tmp$novSWEmean,tmp$decSWEmean,
                tmp$janSWEmean,tmp$febSWEmean)
	cols <- seq(3, 19, 2)
	for (j in 1:9) {
		x <- xlist[,j] # Get the x variable from the list
		# Only do the linear model if there are 3 or more non-NA cases
		if (length(y)-sum(is.na(y))>3 & length(x)-sum(is.na(x))>3) {
			lm1 <- lm(y~x, na.action=na.omit)
			scovTs.lm[i,cols[j]] <- lm1$coef[2]
			scovTs.lm[i,cols[j]+1] <- summary(lm1)$coef[8]
		# Otherwise set NAs
		} else {
			scovTs.lm[i,cols[j]] <- NA
			scovTs.lm[i,cols[j]+1] <- NA
		}
	}
}
# What are the mean Betas for the x variables (in variable order 5, 10, 50)?
sapply(scovTs.lm, mean, na.rm=TRUE)
# 4.680507e-02  5.219132e-03 9.283453e-03  9.859399e-02 6.710935e-03  
# 4.413940e-01 2.103900e-01 1.202558e-01 6.448997e-02
# 3.718991e-02 4.815743e-03 1.179093e-02 1.414930e-01 5.251080e-03 
# 3.876121e-01 1.808581e-01 1.024569e-01 5.442215e-02
# 3.545417e-02  3.726690e-03 1.170069e-02 8.701335e-02 3.425994e-03 
# 4.362551e-01 1.865547e-01 1.009632e-01 5.136469e-02

# How many of these regressions are significant (5 10 50 indvar)?
sigp <- as.data.frame(scovTs.lm < 0.05)
sum(sigp$x1Pval, na.rm=TRUE) # 12 10 10 maxswe
sum(sigp$x2Pval, na.rm=TRUE) # 16 17 18 onsetdoy
sum(sigp$x3Pval, na.rm=TRUE) #  7  8  9 preonsetTair
sum(sigp$x4Pval, na.rm=TRUE) #  9 14 14 scovmat
sum(sigp$x5Pval, na.rm=TRUE) #  6  6  5 meltdoy
sum(sigp$x6Pval, na.rm=TRUE) # 23 28 30 novSWEmean
sum(sigp$x7Pval, na.rm=TRUE) # 37 39 29 decSWEmean
sum(sigp$x8Pval, na.rm=TRUE) # 40 31 29 janSWEmean
sum(sigp$x9Pval, na.rm=TRUE) # 14  8 13 febSWEmean

# BELOW-SNOW VWC - Create a data frame to store regression values
scovVWC.lm <- data.frame(site=sites,
			       nyrs=as.vector(table(reg.dat$siteClim)),
			       x1Beta=0,x1Pval=0,x2Beta=0,x2Pval=0,
			       x3Beta=0,x3Pval=0,x4Beta=0,x4Pval=0,
			       x5Beta=0,x5Pval=0,x6Beta=0,x6Pval=0,
			       x7Beta=0,x7Pval=0,x8Beta=0,x8Pval=0,
                   x9Beta=0,x9Pval=0)
		    
for (i in 1:length(sites)) {
	# Subset for the site
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	# Assign the independent variable
	y <- tmp$jfmVWC5mean
	# Create a list of independent variables
	# peak SWE, onsetdoy, preonset Tair, nov & dec SWE, etc
	xlist <- data.frame(tmp$maxswe,tmp$onsetdoy,tmp$preonsetTair,
			    tmp$scovmat,tmp$meltdoy,
			    tmp$novSWEmean,tmp$decSWEmean,tmp$janSWEmean,tmp$febSWEmean)
	cols <- seq(3, 19, 2)
	for (j in 1:9) {
		x <- xlist[,j] # Get the x variable from the list
		# Only do the linear model if there are 3 or more non-NA cases
		if (length(y)-sum(is.na(y))>3 & length(x)-sum(is.na(x))>3) {
			lm1 <- lm(y~x, na.action=na.omit)
			scovVWC.lm[i,cols[j]] <- lm1$coef[2]
			scovVWC.lm[i,cols[j]+1] <- summary(lm1)$coef[8]
		# Otherwise set NAs
		} else {
			scovVWC.lm[i,cols[j]] <- NA
			scovVWC.lm[i,cols[j]+1] <- NA
		}
	}
}
# What are the mean Betas for the x variables?
sapply(scovVWC.lm, mean, na.rm=TRUE)
# 1.071987e-02 2.519752e-03 3.541852e-03 7.933916e-02 2.214987e-03 
# 1.091579e-01 5.305423e-02 3.201914e-02 1.658974e-02
# 1.186218e-02 1.956167e-03 4.476801e-03 5.881795e-02 2.705888e-03 
# 1.243848e-01 5.420011e-02 3.192925e-02 1.781739e-02
# 7.113975e-03 6.892734e-04 8.616975e-03 7.631356e-02 2.362200e-03 
# 1.283124e-01 4.660537e-02 2.534936e-02 1.078104e-02

# How many of these regressions are significant?
sigp <- as.data.frame(scovVWC.lm < 0.05)
sum(sigp$x1Pval, na.rm=TRUE) # 10 12  6 maxswe
sum(sigp$x2Pval, na.rm=TRUE) # 18 10 12 onsetdoy
sum(sigp$x3Pval, na.rm=TRUE) #  9 10  5 preonsetTair
sum(sigp$x4Pval, na.rm=TRUE) # 14  9 10 scovmat
sum(sigp$x5Pval, na.rm=TRUE) #  5  9  9 meltdoy
sum(sigp$x6Pval, na.rm=TRUE) # 30 19 17 novSWEmean
sum(sigp$x7Pval, na.rm=TRUE) # 29 48 27 decSWEmean
sum(sigp$x8Pval, na.rm=TRUE) # 40 37 18 janSWEmean
sum(sigp$x9Pval, na.rm=TRUE) # 19 16  9 febSWEmean


# SUMMER (JAS) VWC - Create a data frame to store regression values
jasVWC.lm <- data.frame(site=sites,
			       nyrs=as.vector(table(reg.dat$siteClim)),
			       x1Beta=0,x1Pval=0,x2Beta=0,x2Pval=0,
			       x3Beta=0,x3Pval=0,x4Beta=0,x4Pval=0,
			       x5Beta=0,x5Pval=0)		    
for (i in 1:length(sites)) {
	# Subset for the site
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	# Assign the independent variable
	y <- tmp$jasVWC50mean
	# Create a list of independent variables
	# peak SWE, meltdoy, JAS Tair, JAS precip, winter Tsoil
	xlist <- data.frame(tmp$maxswe,tmp$meltdoy,tmp$jasMAT,
			    tmp$JASprecip,tmp$jfmTs5mean)
	cols <- seq(3, 11, 2)
	for (j in 1:5) {
		x <- xlist[,j] # Get the x variable from the list
		# Only do the linear model if there are 3 or more non-NA cases
		if (length(y)-sum(is.na(y))>3 & length(x)-sum(is.na(x))>3) {
			lm1 <- lm(y~x, na.action=na.omit)
			jasVWC.lm[i,cols[j]] <- lm1$coef[2]
			jasVWC.lm[i,cols[j]+1] <- summary(lm1)$coef[8]
		# Otherwise set NAs
		} else {
			jasVWC.lm[i,cols[j]] <- NA
			jasVWC.lm[i,cols[j]+1] <- NA
		}
	}
}
# What are the mean Betas for the x variables?
sapply(jasVWC.lm, mean, na.rm=TRUE)
# 0.005481079 0.001591878 -0.012759328 0.023487758 0.044567477 
# 0.006191316 0.001969057 -0.011746426 0.016697552 0.048679349 
#0.008187641 0.002715287 -0.011602165 0.005269692 0.058009332 

# How many of these regressions are significant?
sigp <- as.data.frame(jasVWC.lm < 0.05)
sum(sigp$x1Pval, na.rm=TRUE) # 13 16 22 maxswe
sum(sigp$x2Pval, na.rm=TRUE) #  9 12 16 meltdoy
sum(sigp$x3Pval, na.rm=TRUE) #  8 13 12 jasTair 
sum(sigp$x4Pval, na.rm=TRUE) # 26 18  6 jasPrecip
sum(sigp$x5Pval, na.rm=TRUE) # 10  9  7 jfmTs5mean

