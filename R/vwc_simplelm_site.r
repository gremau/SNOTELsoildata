source('getdata.r')

# Make some dataframes to analyze with linear models
reg.dat <- cbind(climData.sub[,c('siteClim','yearClim','maxswe','meltdoy',
			     'onsetdoy','JASprecip','junPrecip','julPrecip',
			     'augPrecip','sepPrecip','jasMAT','julTairMean',
			     'augTairMean','sepTairMean','decSWEmean',
			     'jfmMAT')],
	     soilVWCData[,c('jasVWC5mean','julVWC5mean','augVWC5mean',
			    'sepVWC5mean','jasVWC20mean','julVWC20mean',
			    'augVWC20mean','sepVWC20mean','jasVWC50mean',
			    'julVWC50mean','augVWC50mean','sepVWC50mean',
			    'preonsetVWC20','jfmVWC20mean')])
# Create a list of lmLists for the simple linear regressions,
# but first change the precip columns
reg.dat$julPrecip <- (reg.dat$junPrecip + reg.dat$julPrecip)
reg.dat$augPrecip <- (reg.dat$julPrecip + reg.dat$augPrecip)
reg.dat$sepPrecip <- (reg.dat$augPrecip + reg.dat$sepPrecip)


# Screw it, nlme doesn't give the same results as lm
sites <- unique(reg.dat$siteClim)

# First do jas 5cm depth
simpledat.jas5 <- data.frame(site=sites,nyrs=as.vector(table(reg.dat$siteClim)),
		     precInt=0,precIntPval=0,precBeta=0,precPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     matInt=0,matIntPval=0,matBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$jasVWC5mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$JASprecip))>3) {
			lm1 <- lm(jasVWC5mean~JASprecip,data=tmp,
				  na.action=na.omit)
			simpledat.jas5[i,3] <- lm1$coef[1]
			simpledat.jas5[i,4] <- summary(lm1)$coef[7]
			simpledat.jas5[i,5] <- lm1$coef[2]
		       	simpledat.jas5[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.jas5[i,3] <- NA
			simpledat.jas5[i,4] <- NA
			simpledat.jas5[i,5] <- NA
			simpledat.jas5[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(jasVWC5mean~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.jas5[i,7] <- lm2$coef[1]
			simpledat.jas5[i,8] <- summary(lm2)$coef[7]
			simpledat.jas5[i,9] <- lm2$coef[2]
			simpledat.jas5[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.jas5[i,7] <- NA
			simpledat.jas5[i,8] <- NA
			simpledat.jas5[i,9] <- NA
			simpledat.jas5[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(jasVWC5mean~meltdoy,data=tmp,
				  na.action=na.omit)
			simpledat.jas5[i,11] <- lm3$coef[1]
			simpledat.jas5[i,12] <- summary(lm3)$coef[7]
			simpledat.jas5[i,13] <- lm3$coef[2]
			simpledat.jas5[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.jas5[i,11] <- NA
			simpledat.jas5[i,12] <- NA
			simpledat.jas5[i,13] <- NA
			simpledat.jas5[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(jasVWC5mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.jas5[i,15] <- lm4$coef[1]
			simpledat.jas5[i,16] <- summary(lm4)$coef[7]
			simpledat.jas5[i,17] <- lm4$coef[2]
			simpledat.jas5[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.jas5[i,15] <- NA
			simpledat.jas5[i,16] <- NA
			simpledat.jas5[i,17] <- NA
			simpledat.jas5[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$jasMAT))>3) {
			lm5 <- lm(jasVWC5mean~jasMAT,data=tmp,
				  na.action=na.omit)
			simpledat.jas5[i,19] <- lm5$coef[1]
			simpledat.jas5[i,20] <- summary(lm5)$coef[7]
			simpledat.jas5[i,21] <- lm5$coef[2]
			simpledat.jas5[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.jas5[i,19] <- NA
			simpledat.jas5[i,20] <- NA
			simpledat.jas5[i,21] <- NA
			simpledat.jas5[i,22] <- NA
		}
	} else {
		simpledat.jas5[i,-(1:2)] <- NA
	}
}


# Now do 20cm JAS
simpledat.jas20 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     precInt=0,precIntPval=0,precBeta=0,precPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     matInt=0,matIntPval=0,matBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$jasVWC20mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$JASprecip))>3) {
			lm1 <- lm(jasVWC20mean~JASprecip,data=tmp,
				  na.action=na.omit)
			simpledat.jas20[i,3] <- lm1$coef[1]
			simpledat.jas20[i,4] <- summary(lm1)$coef[7]
			simpledat.jas20[i,5] <- lm1$coef[2]
		       	simpledat.jas20[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.jas20[i,3] <- NA
			simpledat.jas20[i,4] <- NA
			simpledat.jas20[i,5] <- NA
			simpledat.jas20[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(jasVWC20mean~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.jas20[i,7] <- lm2$coef[1]
			simpledat.jas20[i,8] <- summary(lm2)$coef[7]
			simpledat.jas20[i,9] <- lm2$coef[2]
			simpledat.jas20[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.jas20[i,7] <- NA
			simpledat.jas20[i,8] <- NA
			simpledat.jas20[i,9] <- NA
			simpledat.jas20[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(jasVWC20mean~meltdoy,data=tmp,
				  na.action=na.omit)
			simpledat.jas20[i,11] <- lm3$coef[1]
			simpledat.jas20[i,12] <- summary(lm3)$coef[7]
			simpledat.jas20[i,13] <- lm3$coef[2]
			simpledat.jas20[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.jas20[i,11] <- NA
			simpledat.jas20[i,12] <- NA
			simpledat.jas20[i,13] <- NA
			simpledat.jas20[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(jasVWC20mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.jas20[i,15] <- lm4$coef[1]
			simpledat.jas20[i,16] <- summary(lm4)$coef[7]
			simpledat.jas20[i,17] <- lm4$coef[2]
			simpledat.jas20[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.jas20[i,15] <- NA
			simpledat.jas20[i,16] <- NA
			simpledat.jas20[i,17] <- NA
			simpledat.jas20[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$jasMAT))>3) {
			lm5 <- lm(jasVWC20mean~jasMAT,data=tmp,
				  na.action=na.omit)
			simpledat.jas20[i,19] <- lm5$coef[1]
			simpledat.jas20[i,20] <- summary(lm5)$coef[7]
			simpledat.jas20[i,21] <- lm5$coef[2]
			simpledat.jas20[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.jas20[i,19] <- NA
			simpledat.jas20[i,20] <- NA
			simpledat.jas20[i,21] <- NA
			simpledat.jas20[i,22] <- NA
		}
	} else {
		simpledat.jas20[i,-(1:2)] <- NA
	}
}


# Now do JAS 50cm 
simpledat.jas50 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     precInt=0,precIntPval=0,precBeta=0,precPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     matInt=0,matIntPval=0,matBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$jasVWC50mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$JASprecip))>3) {
			lm1 <- lm(jasVWC50mean~JASprecip,data=tmp,
				  na.action=na.omit)
			simpledat.jas50[i,3] <- lm1$coef[1]
			simpledat.jas50[i,4] <- summary(lm1)$coef[7]
			simpledat.jas50[i,5] <- lm1$coef[2]
		       	simpledat.jas50[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.jas50[i,3] <- NA
			simpledat.jas50[i,4] <- NA
			simpledat.jas50[i,5] <- NA
			simpledat.jas50[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(jasVWC50mean~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.jas50[i,7] <- lm2$coef[1]
			simpledat.jas50[i,8] <- summary(lm2)$coef[7]
			simpledat.jas50[i,9] <- lm2$coef[2]
			simpledat.jas50[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.jas50[i,7] <- NA
			simpledat.jas50[i,8] <- NA
			simpledat.jas50[i,9] <- NA
			simpledat.jas50[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(jasVWC50mean~meltdoy,data=tmp,
				  na.action=na.omit)
			#plot(tmp$meltdoy,tmp$jasVWC50mean,main=sites[i])
			#abline(lm3)
			#yo <- readline('Press return for next plot')
			simpledat.jas50[i,11] <- lm3$coef[1]
			simpledat.jas50[i,12] <- summary(lm3)$coef[7]
			simpledat.jas50[i,13] <- lm3$coef[2]
			simpledat.jas50[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.jas50[i,11] <- NA
			simpledat.jas50[i,12] <- NA
			simpledat.jas50[i,13] <- NA
			simpledat.jas50[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(jasVWC50mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.jas50[i,15] <- lm4$coef[1]
			simpledat.jas50[i,16] <- summary(lm4)$coef[7]
			simpledat.jas50[i,17] <- lm4$coef[2]
			simpledat.jas50[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.jas50[i,15] <- NA
			simpledat.jas50[i,16] <- NA
			simpledat.jas50[i,17] <- NA
			simpledat.jas50[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$jasMAT))>3) {
			lm5 <- lm(jasVWC50mean~jasMAT,data=tmp,
				  na.action=na.omit)
			simpledat.jas50[i,19] <- lm5$coef[1]
			simpledat.jas50[i,20] <- summary(lm5)$coef[7]
			simpledat.jas50[i,21] <- lm5$coef[2]
			simpledat.jas50[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.jas50[i,19] <- NA
			simpledat.jas50[i,20] <- NA
			simpledat.jas50[i,21] <- NA
			simpledat.jas50[i,22] <- NA
		}
	} else {
		simpledat.jas50[i,-(1:2)] <- NA
	}
}


# Now do JFM 20cm 
simpledat.jfm20 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     precInt=0,precIntPval=0,precBeta=0,precPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     matInt=0,matIntPval=0,matBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$jfmVWC20mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm1 <- lm(jfmVWC20mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.jfm20[i,3] <- lm1$coef[1]
			simpledat.jfm20[i,4] <- summary(lm1)$coef[7]
			simpledat.jfm20[i,5] <- lm1$coef[2]
		       	simpledat.jfm20[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.jfm20[i,3] <- NA
			simpledat.jfm20[i,4] <- NA
			simpledat.jfm20[i,5] <- NA
			simpledat.jfm20[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$preonsetVWC20))>3) {
			lm2 <- lm(jfmVWC20mean~preonsetVWC20,data=tmp,
			  na.action=na.omit)
			simpledat.jfm20[i,7] <- lm2$coef[1]
			simpledat.jfm20[i,8] <- summary(lm2)$coef[7]
			simpledat.jfm20[i,9] <- lm2$coef[2]
			simpledat.jfm20[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.jfm20[i,7] <- NA
			simpledat.jfm20[i,8] <- NA
			simpledat.jfm20[i,9] <- NA
			simpledat.jfm20[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(jfmVWC20mean~meltdoy,data=tmp,
				  na.action=na.omit)
			#plot(tmp$meltdoy,tmp$jasVWC50mean,main=sites[i])
			#abline(lm3)
			#yo <- readline('Press return for next plot')
			simpledat.jfm20[i,11] <- lm3$coef[1]
			simpledat.jfm20[i,12] <- summary(lm3)$coef[7]
			simpledat.jfm20[i,13] <- lm3$coef[2]
			simpledat.jfm20[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.jfm20[i,11] <- NA
			simpledat.jfm20[i,12] <- NA
			simpledat.jfm20[i,13] <- NA
			simpledat.jfm20[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$decSWEmean))>3) {
			lm4 <- lm(jfmVWC20mean~decSWEmean,data=tmp,
				  na.action=na.omit)
			simpledat.jfm20[i,15] <- lm4$coef[1]
			simpledat.jfm20[i,16] <- summary(lm4)$coef[7]
			simpledat.jfm20[i,17] <- lm4$coef[2]
			simpledat.jfm20[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.jfm20[i,15] <- NA
			simpledat.jfm20[i,16] <- NA
			simpledat.jfm20[i,17] <- NA
			simpledat.jfm20[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$jfmMAT))>3) {
			lm5 <- lm(jfmVWC20mean~jfmMAT,data=tmp,
				  na.action=na.omit)
			simpledat.jfm20[i,19] <- lm5$coef[1]
			simpledat.jfm20[i,20] <- summary(lm5)$coef[7]
			simpledat.jfm20[i,21] <- lm5$coef[2]
			simpledat.jfm20[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.jfm20[i,19] <- NA
			simpledat.jfm20[i,20] <- NA
			simpledat.jfm20[i,21] <- NA
			simpledat.jfm20[i,22] <- NA
		}
	} else {
		simpledat.jfm20[i,-(1:2)] <- NA
	}
}

# Now do 5cm August
simpledat.aug5 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     precInt=0,precIntPval=0,precBeta=0,precPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     matInt=0,matIntPval=0,matBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$augVWC5mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$augPrecip))>3) {
			lm1 <- lm(augVWC5mean~augPrecip,data=tmp,
				  na.action=na.omit)
			simpledat.aug5[i,3] <- lm1$coef[1]
			simpledat.aug5[i,4] <- summary(lm1)$coef[7]
			simpledat.aug5[i,5] <- lm1$coef[2]
		       	simpledat.aug5[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.aug5[i,3] <- NA
			simpledat.aug5[i,4] <- NA
			simpledat.aug5[i,5] <- NA
			simpledat.aug5[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(augVWC5mean~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.aug5[i,7] <- lm2$coef[1]
			simpledat.aug5[i,8] <- summary(lm2)$coef[7]
			simpledat.aug5[i,9] <- lm2$coef[2]
			simpledat.aug5[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.aug5[i,7] <- NA
			simpledat.aug5[i,8] <- NA
			simpledat.aug5[i,9] <- NA
			simpledat.aug5[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(augVWC5mean~meltdoy,data=tmp,
				  na.action=na.omit)
			simpledat.aug5[i,11] <- lm3$coef[1]
			simpledat.aug5[i,12] <- summary(lm3)$coef[7]
			simpledat.aug5[i,13] <- lm3$coef[2]
			simpledat.aug5[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.aug5[i,11] <- NA
			simpledat.aug5[i,12] <- NA
			simpledat.aug5[i,13] <- NA
			simpledat.aug5[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(augVWC5mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.aug5[i,15] <- lm4$coef[1]
			simpledat.aug5[i,16] <- summary(lm4)$coef[7]
			simpledat.aug5[i,17] <- lm4$coef[2]
			simpledat.aug5[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.aug5[i,15] <- NA
			simpledat.aug5[i,16] <- NA
			simpledat.aug5[i,17] <- NA
			simpledat.aug5[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$jasMAT))>3) {
			lm5 <- lm(augVWC5mean~jasMAT,data=tmp,
				  na.action=na.omit)
			simpledat.aug5[i,19] <- lm5$coef[1]
			simpledat.aug5[i,20] <- summary(lm5)$coef[7]
			simpledat.aug5[i,21] <- lm5$coef[2]
			simpledat.aug5[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.aug5[i,19] <- NA
			simpledat.aug5[i,20] <- NA
			simpledat.aug5[i,21] <- NA
			simpledat.aug5[i,22] <- NA
		}
	} else {
		simpledat.aug5[i,-(1:2)] <- NA
	}
}

# Now do 20cm August
simpledat.aug20 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     precInt=0,precIntPval=0,precBeta=0,precPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     matInt=0,matIntPval=0,matBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$augVWC20mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$augPrecip))>3) {
			lm1 <- lm(augVWC20mean~augPrecip,data=tmp,
				  na.action=na.omit)
			simpledat.aug20[i,3] <- lm1$coef[1]
			simpledat.aug20[i,4] <- summary(lm1)$coef[7]
			simpledat.aug20[i,5] <- lm1$coef[2]
		       	simpledat.aug20[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.aug20[i,3] <- NA
			simpledat.aug20[i,4] <- NA
			simpledat.aug20[i,5] <- NA
			simpledat.aug20[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(augVWC20mean~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.aug20[i,7] <- lm2$coef[1]
			simpledat.aug20[i,8] <- summary(lm2)$coef[7]
			simpledat.aug20[i,9] <- lm2$coef[2]
			simpledat.aug20[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.aug20[i,7] <- NA
			simpledat.aug20[i,8] <- NA
			simpledat.aug20[i,9] <- NA
			simpledat.aug20[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(augVWC20mean~meltdoy,data=tmp,
				  na.action=na.omit)
			simpledat.aug20[i,11] <- lm3$coef[1]
			simpledat.aug20[i,12] <- summary(lm3)$coef[7]
			simpledat.aug20[i,13] <- lm3$coef[2]
			simpledat.aug20[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.aug20[i,11] <- NA
			simpledat.aug20[i,12] <- NA
			simpledat.aug20[i,13] <- NA
			simpledat.aug20[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(augVWC20mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.aug20[i,15] <- lm4$coef[1]
			simpledat.aug20[i,16] <- summary(lm4)$coef[7]
			simpledat.aug20[i,17] <- lm4$coef[2]
			simpledat.aug20[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.aug20[i,15] <- NA
			simpledat.aug20[i,16] <- NA
			simpledat.aug20[i,17] <- NA
			simpledat.aug20[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$jasMAT))>3) {
			lm5 <- lm(augVWC20mean~jasMAT,data=tmp,
				  na.action=na.omit)
			simpledat.aug20[i,19] <- lm5$coef[1]
			simpledat.aug20[i,20] <- summary(lm5)$coef[7]
			simpledat.aug20[i,21] <- lm5$coef[2]
			simpledat.aug20[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.aug20[i,19] <- NA
			simpledat.aug20[i,20] <- NA
			simpledat.aug20[i,21] <- NA
			simpledat.aug20[i,22] <- NA
		}
	} else {
		simpledat.aug20[i,-(1:2)] <- NA
	}
}

# Now do 50cm August
simpledat.aug50 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     precInt=0,precIntPval=0,precBeta=0,precPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     matInt=0,matIntPval=0,matBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$augVWC50mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$augPrecip))>3) {
			lm1 <- lm(augVWC50mean~augPrecip,data=tmp,
				  na.action=na.omit)
			simpledat.aug50[i,3] <- lm1$coef[1]
			simpledat.aug50[i,4] <- summary(lm1)$coef[7]
			simpledat.aug50[i,5] <- lm1$coef[2]
		       	simpledat.aug50[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.aug50[i,3] <- NA
			simpledat.aug50[i,4] <- NA
			simpledat.aug50[i,5] <- NA
			simpledat.aug50[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(augVWC50mean~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.aug50[i,7] <- lm2$coef[1]
			simpledat.aug50[i,8] <- summary(lm2)$coef[7]
			simpledat.aug50[i,9] <- lm2$coef[2]
			simpledat.aug50[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.aug50[i,7] <- NA
			simpledat.aug50[i,8] <- NA
			simpledat.aug50[i,9] <- NA
			simpledat.aug50[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(augVWC50mean~meltdoy,data=tmp,
				  na.action=na.omit)
			simpledat.aug50[i,11] <- lm3$coef[1]
			simpledat.aug50[i,12] <- summary(lm3)$coef[7]
			simpledat.aug50[i,13] <- lm3$coef[2]
			simpledat.aug50[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.aug50[i,11] <- NA
			simpledat.aug50[i,12] <- NA
			simpledat.aug50[i,13] <- NA
			simpledat.aug50[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(augVWC50mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.aug50[i,15] <- lm4$coef[1]
			simpledat.aug50[i,16] <- summary(lm4)$coef[7]
			simpledat.aug50[i,17] <- lm4$coef[2]
			simpledat.aug50[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.aug50[i,15] <- NA
			simpledat.aug50[i,16] <- NA
			simpledat.aug50[i,17] <- NA
			simpledat.aug50[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$jasMAT))>3) {
			lm5 <- lm(augVWC50mean~jasMAT,data=tmp,
				  na.action=na.omit)
			simpledat.aug50[i,19] <- lm5$coef[1]
			simpledat.aug50[i,20] <- summary(lm5)$coef[7]
			simpledat.aug50[i,21] <- lm5$coef[2]
			simpledat.aug50[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.aug50[i,19] <- NA
			simpledat.aug50[i,20] <- NA
			simpledat.aug50[i,21] <- NA
			simpledat.aug50[i,22] <- NA
		}
	} else {
		simpledat.aug50[i,-(1:2)] <- NA
	}
}

rm(lm1, lm2, lm3, lm4, lm5, tmp, sites, reg.dat, i)


# Summarize the coefficients and significance of all these
rm(simpledat.all, simpledat.allnum)
datlist <- ls(pattern='simpledat')

simpledat.all <- data.frame(ydat=datlist,nyrs=0,
		     precInt=0,precIntPval=0,precBeta=0,precPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,
		     onsetPval=0,meltInt=0,meltIntPval=0,meltBeta=0,
		     meltPval=0,maxsweInt=0,maxsweIntPval=0,
		     maxsweBeta=0,maxswePval=0,matInt=0,matIntPval=0,
		     matBeta=0,matPval=0)
simpledat.allnum <- data.frame(ydat=datlist,precIntN=0,precBetaN=0,
	  		       onsetIntN=0,onsetBetaN=0,
			       meltIntN=0,meltBetaN=0,
			       maxsweIntN=0,maxsweBetaN=0,
			       matIntN=0,matBetaN=0)

for (i in 1:length(datlist)) {
	dat <- get(datlist[i])
	simpledat.all$nyrs[i] <- mean(dat[,2], na.rm=TRUE)
	numidx <- 2
	for(j in seq(4,22,by=2)) {
		sigtest <- dat[,j]<0.05
		totalreg <- sum(!is.na(sigtest))
		numsig <- sum(sigtest, na.rm=TRUE)
		sigcoeff <- mean(dat[sigtest,j-1], na.rm=TRUE)
		simpledat.all[i,j] <- numsig
		simpledat.all[i,j-1] <- sigcoeff
		simpledat.allnum[i,numidx] <- totalreg
		numidx <- numidx+1
	}
}






