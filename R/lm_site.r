source('getdata.r')
library('nlme')

options(contrasts=c("contr.sum","contr.poly"))


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
wgtZ.lm <- data.frame(site=sites,
  		       nyrs=as.vector(table(reg.dat$siteClim)),
			       x1Z=0,x1df=0,x2Z=0,x2df=0,
			       x3Z=0,x3df=0,x4Z=0,x4df=0,
			       x5Z=0,x5df=0,x6Z=0,x6df=0,
			       x7Z=0,x7df=0,x8Z=0,x8df=0,
                   x9Z=0,x9df=0)
for (i in 1:length(sites)) {
	# Subset for the site
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	# Assign the independent variable
	y <- tmp$snowcovTs50mean
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
      # Save the degrees of freedom to weight the combined test
      wgtZ.lm[i,cols[j]+1] <- lm1$df.residual
		# Otherwise set NAs
		} else {
			scovTs.lm[i,cols[j]] <- NA
			scovTs.lm[i,cols[j]+1] <- NA
		}
	}
}
# Get the mean Betas for the x variables 
meandat <- sapply(scovTs.lm, mean, na.rm=TRUE)

# How many of these regressions are significant (5 10 50 indvar)?
sigp <- as.data.frame(scovTs.lm < 0.05)
pdat <- sapply(sigp, sum, na.rm=TRUE)

#Show a matrix of these values
cbind(meandat[3:20], pdat[3:20])
# maxswe
# onsetdoy
# preonsetTair
# scovmat
# meltdoy
# novSWEmean
# decSWEmean
# janSWEmean
# febSWEmean

# Now do the weighted z transform (Whitlock 2005) to test if these
# are significant in a combine sense

for (i in 1:length(sites)) {
  # Subset pvalue data for the site
  tmp <- scovTs.lm[scovTs.lm$site==sites[i],]
  # Create a list of independent variables
  # peak SWE, onsetdoy, preonset Tair, nov & dec SWE
  #xlist <- data.frame(tmp$maxswe,tmp$onsetdoy,tmp$preonsetTair,
  #                    tmp$scovmat,tmp$meltdoy,tmp$novSWEmean,tmp$decSWEmean,
  #                    tmp$janSWEmean,tmp$febSWEmean)
  cols <- seq(3, 19, 2)
  for (j in 1:9) {
    p_mu <- meandat[cols[j] + 1] # Get the mean p value for that row
    wgtZ.lm[i,cols[j]] <- (p_mu - tmp[cols[j] + 1])/1
  }
}
# Now we need to compute the numerator and denominator for the test
#browser()

# Now the nlme analysis
y <- reg.dat$jasVWC50mean #reg.dat$snowcovTs5mean #jfm
# df <- data.frame(reg.dat$maxswe,reg.dat$onsetdoy,
#                  reg.dat$preonsetTair,reg.dat$scovmat,reg.dat$meltdoy,
#                  reg.dat$novSWEmean,reg.dat$decSWEmean,reg.dat$janSWEmean,
#                  reg.dat$febSWEmean)
df <- data.frame(reg.dat$maxswe,reg.dat$meltdoy,reg.dat$jasMAT,
                    reg.dat$JASprecip,reg.dat$jfmTs5mean)
sitevec <- reg.dat$siteClim
for (i in 1:9) {
  x <- df[,i] # Get the x variable from the list
  mod <- lme(y ~ x, random=~-1 + x| sitevec,
          na.action = na.omit)
  print(i)
  print(anova(mod))
}

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
	y <- tmp$jfmVWC50mean
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
# Get the mean Betas for the x variables 
meandat <- sapply(scovVWC.lm, mean, na.rm=TRUE)

# How many of these regressions are significant (5 10 50 indvar)?
sigp <- as.data.frame(scovVWC.lm < 0.05)
pdat <- sapply(sigp, sum, na.rm=TRUE)

#Show a matrix of these values
cbind(meandat[3:20], pdat[3:20])
# maxswe
# onsetdoy
# preonsetTair
# scovmat
# meltdoy
# novSWEmean
# decSWEmean
# janSWEmean
# febSWEmean

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
# Get the mean Betas for the x variables 
meandat <- sapply(jasVWC.lm, mean, na.rm=TRUE)

# How many of these regressions are significant (5 10 50 indvar)?
sigp <- as.data.frame(jasVWC.lm < 0.05)
pdat <- sapply(sigp, sum, na.rm=TRUE)

#Show a matrix of these values
cbind(meandat[3:12], pdat[3:12])
# maxswe
# meltdoy
# jasTair
# jasPrecip
# jfmTs5mean

