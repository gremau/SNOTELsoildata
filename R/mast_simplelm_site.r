source('getdata.r')

climData.sub <- climData.sub[order(climData.sub$siteClim, 
				   climData.sub$yearClim),]
soilTData <- soilTData[order(soilTData$siteTsoil,
			     soilTData$yearTsoil),]

# Make some dataframes to analyze with linear models
reg.dat <- cbind(climData.sub[,c('siteClim','yearClim','maxswe','meltdoy',
			     'onsetdoy','maat','totaldaysSC','decTairMean',
			     'janTairMean','febTairMean','mayTairMean',
			     'augTairMean','novTairMean','freemat',
			     'scovmat','novSWEmed','decSWEmed','janSWEmed',
			     'febSWEmed','preonsetTair','decSWEmean')],
	     soilTData[,c('mast5cm','mast20cm','mast50cm','decTs20mean',
			  'janTs20mean','febTs20mean', 'mayTs20mean',
			  'augTs20mean','novTs20mean','ondTs20mean',
			  'jfmTs20mean','amjTs20mean','jasTs20mean',
			  'snowfreeTs20mean','snowcovTs20mean')])
# Create a list of lmLists for the simple linear regressions,
# but first change the precip columns
reg.dat$diff5cm <- (reg.dat$maat - reg.dat$mast5cm)
reg.dat$diff20cm <- (reg.dat$maat - reg.dat$mast20cm)
reg.dat$diff50cm <- (reg.dat$maat -  reg.dat$mast50cm)
reg.dat$ondMAT <- (climData.sub$octTairMean+climData.sub$novTairMean+
		   climData.sub$decTairMean)
reg.dat$jfmMAT <- (climData.sub$janTairMean+climData.sub$febTairMean+
		   climData.sub$marTairMean)
reg.dat$amjMAT <- (climData.sub$aprTairMean+climData.sub$mayTairMean+
		   climData.sub$junTairMean)

# Screw it, nlme doesn't give the same results as lm
sites <- unique(reg.dat$siteClim)

# First do mast @ 20cm depth
simpledat.mast20 <- data.frame(site=sites,
			       nyrs=as.vector(table(reg.dat$siteClim)),
		     sdaysInt=0,sdaysIntPval=0,sdaysBeta=0,sdaysPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     matInt=0,matIntPval=0,matBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$mast20cm))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$totaldaysSC))>3) {
			lm1 <- lm(mast20cm~totaldaysSC,data=tmp,
				  na.action=na.omit)
			simpledat.mast20[i,3] <- lm1$coef[1]
			simpledat.mast20[i,4] <- summary(lm1)$coef[7]
			simpledat.mast20[i,5] <- lm1$coef[2]
		       	simpledat.mast20[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.mast20[i,3] <- NA
			simpledat.mast20[i,4] <- NA
			simpledat.mast20[i,5] <- NA
			simpledat.mast20[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(mast20cm~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.mast20[i,7] <- lm2$coef[1]
			simpledat.mast20[i,8] <- summary(lm2)$coef[7]
			simpledat.mast20[i,9] <- lm2$coef[2]
			simpledat.mast20[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.mast20[i,7] <- NA
			simpledat.mast20[i,8] <- NA
			simpledat.mast20[i,9] <- NA
			simpledat.mast20[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(mast20cm~meltdoy,data=tmp,
				  na.action=na.omit)
			simpledat.mast20[i,11] <- lm3$coef[1]
			simpledat.mast20[i,12] <- summary(lm3)$coef[7]
			simpledat.mast20[i,13] <- lm3$coef[2]
			simpledat.mast20[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.mast20[i,11] <- NA
			simpledat.mast20[i,12] <- NA
			simpledat.mast20[i,13] <- NA
			simpledat.mast20[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(mast20cm~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.mast20[i,15] <- lm4$coef[1]
			simpledat.mast20[i,16] <- summary(lm4)$coef[7]
			simpledat.mast20[i,17] <- lm4$coef[2]
			simpledat.mast20[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.mast20[i,15] <- NA
			simpledat.mast20[i,16] <- NA
			simpledat.mast20[i,17] <- NA
			simpledat.mast20[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maat))>3) {
			lm5 <- lm(mast20cm~maat,data=tmp,
				  na.action=na.omit)
			simpledat.mast20[i,19] <- lm5$coef[1]
			simpledat.mast20[i,20] <- summary(lm5)$coef[7]
			simpledat.mast20[i,21] <- lm5$coef[2]
			simpledat.mast20[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.mast20[i,19] <- NA
			simpledat.mast20[i,20] <- NA
			simpledat.mast20[i,21] <- NA
			simpledat.mast20[i,22] <- NA
		}
	} else {
		simpledat.mast20[i,-(1:2)] <- NA
	}
}


# Now do diff at 20cm
simpledat.diff20 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     sdaysInt=0,sdaysIntPval=0,sdaysBeta=0,sdaysPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     matInt=0,matIntPval=0,matBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$diff20cm))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$totaldaysSC))>3) {
			lm1 <- lm(diff20cm~totaldaysSC,data=tmp,
				  na.action=na.omit)
			simpledat.diff20[i,3] <- lm1$coef[1]
			simpledat.diff20[i,4] <- summary(lm1)$coef[7]
			simpledat.diff20[i,5] <- lm1$coef[2]
		       	simpledat.diff20[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.diff20[i,3] <- NA
			simpledat.diff20[i,4] <- NA
			simpledat.diff20[i,5] <- NA
			simpledat.diff20[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(diff20cm~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.diff20[i,7] <- lm2$coef[1]
			simpledat.diff20[i,8] <- summary(lm2)$coef[7]
			simpledat.diff20[i,9] <- lm2$coef[2]
			simpledat.diff20[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.diff20[i,7] <- NA
			simpledat.diff20[i,8] <- NA
			simpledat.diff20[i,9] <- NA
			simpledat.diff20[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(diff20cm~meltdoy,data=tmp,
				  na.action=na.omit)
			simpledat.diff20[i,11] <- lm3$coef[1]
			simpledat.diff20[i,12] <- summary(lm3)$coef[7]
			simpledat.diff20[i,13] <- lm3$coef[2]
			simpledat.diff20[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.diff20[i,11] <- NA
			simpledat.diff20[i,12] <- NA
			simpledat.diff20[i,13] <- NA
			simpledat.diff20[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(diff20cm~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.diff20[i,15] <- lm4$coef[1]
			simpledat.diff20[i,16] <- summary(lm4)$coef[7]
			simpledat.diff20[i,17] <- lm4$coef[2]
			simpledat.diff20[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.diff20[i,15] <- NA
			simpledat.diff20[i,16] <- NA
			simpledat.diff20[i,17] <- NA
			simpledat.diff20[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maat))>3) {
			lm5 <- lm(diff20cm~maat,data=tmp,
				  na.action=na.omit)
			simpledat.diff20[i,19] <- lm5$coef[1]
			simpledat.diff20[i,20] <- summary(lm5)$coef[7]
			simpledat.diff20[i,21] <- lm5$coef[2]
			simpledat.diff20[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.diff20[i,19] <- NA
			simpledat.diff20[i,20] <- NA
			simpledat.diff20[i,21] <- NA
			simpledat.diff20[i,22] <- NA
		}
	} else {
		simpledat.diff20[i,-(1:2)] <- NA
	}
}

# Now do jan feb march 20cm soil 
simpledat.jfm20 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     sdaysInt=0,sdaysIntPval=0,sdaysBeta=0,sdaysPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     jfmmatInt=0,jfmmatIntPval=0,jfmmatBeta=0,jfmmatPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$jfmTs20mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$totaldaysSC))>3) {
			lm1 <- lm(jfmTs20mean~totaldaysSC,data=tmp,
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
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(jfmTs20mean~onsetdoy,data=tmp,
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
			lm3 <- lm(jfmTs20mean~meltdoy,data=tmp,
				  na.action=na.omit)
			#plot(tmp$meltdoy,tmp$jfmTs20mean,main=sites[i])
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
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(jfmTs20mean~maxswe,data=tmp,
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
			lm4 <- lm(jfmTs20mean~jfmMAT,data=tmp,
				  na.action=na.omit)
			simpledat.jfm20[i,19] <- lm4$coef[1]
			simpledat.jfm20[i,20] <- summary(lm4)$coef[7]
			simpledat.jfm20[i,21] <- lm4$coef[2]
			simpledat.jfm20[i,22] <- summary(lm4)$coef[8]
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

# Now do april may june 20cm soil 
simpledat.amj20 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     sdaysInt=0,sdaysIntPval=0,sdaysBeta=0,sdaysPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     amjmatInt=0,amjmatIntPval=0,amjmatBeta=0,matPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$amjTs20mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$totaldaysSC))>3) {
			lm1 <- lm(amjTs20mean~totaldaysSC,data=tmp,
				  na.action=na.omit)
			simpledat.amj20[i,3] <- lm1$coef[1]
			simpledat.amj20[i,4] <- summary(lm1)$coef[7]
			simpledat.amj20[i,5] <- lm1$coef[2]
		       	simpledat.amj20[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.amj20[i,3] <- NA
			simpledat.amj20[i,4] <- NA
			simpledat.amj20[i,5] <- NA
			simpledat.amj20[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(amjTs20mean~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.amj20[i,7] <- lm2$coef[1]
			simpledat.amj20[i,8] <- summary(lm2)$coef[7]
			simpledat.amj20[i,9] <- lm2$coef[2]
			simpledat.amj20[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.amj20[i,7] <- NA
			simpledat.amj20[i,8] <- NA
			simpledat.amj20[i,9] <- NA
			simpledat.amj20[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(amjTs20mean~meltdoy,data=tmp,
				  na.action=na.omit)
			simpledat.amj20[i,11] <- lm3$coef[1]
			simpledat.amj20[i,12] <- summary(lm3)$coef[7]
			simpledat.amj20[i,13] <- lm3$coef[2]
			simpledat.amj20[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.amj20[i,11] <- NA
			simpledat.amj20[i,12] <- NA
			simpledat.amj20[i,13] <- NA
			simpledat.amj20[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(amjTs20mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.amj20[i,15] <- lm4$coef[1]
			simpledat.amj20[i,16] <- summary(lm4)$coef[7]
			simpledat.amj20[i,17] <- lm4$coef[2]
			simpledat.amj20[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.amj20[i,15] <- NA
			simpledat.amj20[i,16] <- NA
			simpledat.amj20[i,17] <- NA
			simpledat.amj20[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$amjMAT))>3) {
			lm5 <- lm(amjTs20mean~amjMAT,data=tmp,
				  na.action=na.omit)
			simpledat.amj20[i,19] <- lm5$coef[1]
			simpledat.amj20[i,20] <- summary(lm5)$coef[7]
			simpledat.amj20[i,21] <- lm5$coef[2]
			simpledat.amj20[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.amj20[i,19] <- NA
			simpledat.amj20[i,20] <- NA
			simpledat.amj20[i,21] <- NA
			simpledat.amj20[i,22] <- NA
		}
	} else {
		simpledat.amj20[i,-(1:2)] <- NA
	}
}

# Now do snowcovered 20cm soil 
simpledat.snowcov20 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     preonsInt=0,preonsIntPval=0,preonsBeta=0,preonsPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     decsweInt=0,decsweIntPval=0,decsweBeta=0,decswePval=0,
		     scovmatInt=0,scovmatIntPval=0,scovmatBeta=0,scovPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$snowcovTs20mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$preonsetTair))>3) {
			lm1 <- lm(snowcovTs20mean~preonsetTair,data=tmp,
				  na.action=na.omit)
			simpledat.snowcov20[i,3] <- lm1$coef[1]
			simpledat.snowcov20[i,4] <- summary(lm1)$coef[7]
			simpledat.snowcov20[i,5] <- lm1$coef[2]
		       	simpledat.snowcov20[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.snowcov20[i,3] <- NA
			simpledat.snowcov20[i,4] <- NA
			simpledat.snowcov20[i,5] <- NA
			simpledat.snowcov20[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(snowcovTs20mean~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.snowcov20[i,7] <- lm2$coef[1]
			simpledat.snowcov20[i,8] <- summary(lm2)$coef[7]
			simpledat.snowcov20[i,9] <- lm2$coef[2]
			simpledat.snowcov20[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.snowcov20[i,7] <- NA
			simpledat.snowcov20[i,8] <- NA
			simpledat.snowcov20[i,9] <- NA
			simpledat.snowcov20[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm3 <- lm(snowcovTs20mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.snowcov20[i,11] <- lm3$coef[1]
			simpledat.snowcov20[i,12] <- summary(lm3)$coef[7]
			simpledat.snowcov20[i,13] <- lm3$coef[2]
			simpledat.snowcov20[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.snowcov20[i,11] <- NA
			simpledat.snowcov20[i,12] <- NA
			simpledat.snowcov20[i,13] <- NA
			simpledat.snowcov20[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$decSWEmean))>3) {
			lm4 <- lm(snowcovTs20mean~decSWEmean,data=tmp,
				  na.action=na.omit)
			simpledat.snowcov20[i,15] <- lm4$coef[1]
			simpledat.snowcov20[i,16] <- summary(lm4)$coef[7]
			simpledat.snowcov20[i,17] <- lm4$coef[2]
			simpledat.snowcov20[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.snowcov20[i,15] <- NA
			simpledat.snowcov20[i,16] <- NA
			simpledat.snowcov20[i,17] <- NA
			simpledat.snowcov20[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$scovmat))>3) {
			lm5 <- lm(snowcovTs20mean~scovmat,data=tmp,
				  na.action=na.omit)
			simpledat.snowcov20[i,19] <- lm5$coef[1]
			simpledat.snowcov20[i,20] <- summary(lm5)$coef[7]
			simpledat.snowcov20[i,21] <- lm5$coef[2]
			simpledat.snowcov20[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.snowcov20[i,19] <- NA
			simpledat.snowcov20[i,20] <- NA
			simpledat.snowcov20[i,21] <- NA
			simpledat.snowcov20[i,22] <- NA
		}
	} else {
		simpledat.snowcov20[i,-(1:2)] <- NA
	}
}

# Now do snowfree 20cm soil 
simpledat.snowfree20 <- data.frame(site=sites,
			      nyrs=as.vector(table(reg.dat$siteClim)),
		     sdaysInt=0,sdaysIntPval=0,sdaysBeta=0,sdaysPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,onsetPval=0,
		     meltInt=0,meltIntPval=0,meltBeta=0,meltPval=0,
		     maxsweInt=0,maxsweIntPval=0,maxsweBeta=0,maxswePval=0,
		     freematInt=0,freematIntPval=0,freematBeta=0,freematPval=0)
		    
for (i in 1:length(sites)) {
	tmp <- reg.dat[reg.dat$siteClim==sites[i],]
	if (dim(tmp)[1]-sum(is.na(tmp$snowfreeTs20mean))>3) {
		if (dim(tmp)[1]-sum(is.na(tmp$snowcovTs20mean))>3) {
			lm1 <- lm(snowfreeTs20mean~snowcovTs20mean,data=tmp,
				  na.action=na.omit)
			simpledat.snowfree20[i,3] <- lm1$coef[1]
			simpledat.snowfree20[i,4] <- summary(lm1)$coef[7]
			simpledat.snowfree20[i,5] <- lm1$coef[2]
		       	simpledat.snowfree20[i,6] <- summary(lm1)$coef[8]
		} else {
			simpledat.snowfree20[i,3] <- NA
			simpledat.snowfree20[i,4] <- NA
			simpledat.snowfree20[i,5] <- NA
			simpledat.snowfree20[i,6] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$onsetdoy))>3) {
			lm2 <- lm(snowfreeTs20mean~onsetdoy,data=tmp,
			  na.action=na.omit)
			simpledat.snowfree20[i,7] <- lm2$coef[1]
			simpledat.snowfree20[i,8] <- summary(lm2)$coef[7]
			simpledat.snowfree20[i,9] <- lm2$coef[2]
			simpledat.snowfree20[i,10] <- summary(lm2)$coef[8]
		} else {
			simpledat.snowfree20[i,7] <- NA
			simpledat.snowfree20[i,8] <- NA
			simpledat.snowfree20[i,9] <- NA
			simpledat.snowfree20[i,10] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$meltdoy))>3) {
			lm3 <- lm(snowfreeTs20mean~meltdoy,data=tmp,
				  na.action=na.omit)
			simpledat.snowfree20[i,11] <- lm3$coef[1]
			simpledat.snowfree20[i,12] <- summary(lm3)$coef[7]
			simpledat.snowfree20[i,13] <- lm3$coef[2]
			simpledat.snowfree20[i,14] <- summary(lm3)$coef[8]
		} else {
			simpledat.snowfree20[i,11] <- NA
			simpledat.snowfree20[i,12] <- NA
			simpledat.snowfree20[i,13] <- NA
			simpledat.snowfree20[i,14] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$maxswe))>3) {
			lm4 <- lm(snowfreeTs20mean~maxswe,data=tmp,
				  na.action=na.omit)
			simpledat.snowfree20[i,15] <- lm4$coef[1]
			simpledat.snowfree20[i,16] <- summary(lm4)$coef[7]
			simpledat.snowfree20[i,17] <- lm4$coef[2]
			simpledat.snowfree20[i,18] <- summary(lm4)$coef[8]
		} else {
			simpledat.snowfree20[i,15] <- NA
			simpledat.snowfree20[i,16] <- NA
			simpledat.snowfree20[i,17] <- NA
			simpledat.snowfree20[i,18] <- NA
		}
		if (dim(tmp)[1]-sum(is.na(tmp$freemat))>3) {
			lm5 <- lm(snowfreeTs20mean~freemat,data=tmp,
				  na.action=na.omit)
			simpledat.snowfree20[i,19] <- lm5$coef[1]
			simpledat.snowfree20[i,20] <- summary(lm5)$coef[7]
			simpledat.snowfree20[i,21] <- lm5$coef[2]
			simpledat.snowfree20[i,22] <- summary(lm5)$coef[8]
		} else {
			simpledat.snowfree20[i,19] <- NA
			simpledat.snowfree20[i,20] <- NA
			simpledat.snowfree20[i,21] <- NA
			simpledat.snowfree20[i,22] <- NA
		}
	} else {
		simpledat.snowfree20[i,-(1:2)] <- NA
	}
}

rm(lm1, lm2, lm3, lm4, lm5, tmp, sites, reg.dat, i)


# Summarize the coefficients and significance of all these
rm(simpledat.all, simpledat.allnum)
datlist <- ls(pattern='simpledat')

simpledat.all <- data.frame(ydat=datlist,nyrs=0,
		     sdaysInt=0,sdaysIntPval=0,sdaysBeta=0,sdaysPval=0,
		     onsetInt=0,onsetIntPval=0,onsetBeta=0,
		     onsetPval=0,meltInt=0,meltIntPval=0,meltBeta=0,
		     meltPval=0,maxsweInt=0,maxsweIntPval=0,
		     maxsweBeta=0,maxswePval=0,matInt=0,matIntPval=0,
		     matBeta=0,matPval=0)
simpledat.allnum <- data.frame(ydat=datlist,sdaysIntN=0,sdaysBetaN=0,
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






