source('getdata.r')

# We're interested in whether growing season VWC is related to mean
# air temperature, melt timing, snowpack size, or summer rain.
#
# First lets check whether any predictors are correlated. We already know that
# the snow/temperature variables are correlated, what about precip?
reg1 <- lm(climData$onsetdoy ~ climData$JASprecip)
plot(climData$JASprecip, climData$onsetdoy)
abline(reg1)
summary(reg1)
# Yes, most are highly correlated
#                      Estimate Std. Error t value Pr(>|t|)
# maxswe & JASprecip   -0.21337    0.05287  -4.036 5.51e-05 ***
# maat & JASprecip     -0.08815    0.01032  -8.541   <2e-16 ***
# meltdoy & JASprecip   -0.9491     0.1176  -8.067 8.77e-16 ***
# onsetdoy & JASprecip   0.5569     0.0810   6.875 6.88e-12 ***

# Is JAS vwc predicted by maxswe?
reg2 <- lm(soilVWCData$jasVWC20mean ~ climData.sub$maxswe)
plot(climData.sub$maxswe, soilVWCData$jasVWC20mean)
abline(reg2)
summary(reg2)
# Not really - only a weak correlation
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.2670955  0.0107016  24.958   <2e-16 ***
# climData.sub$maxswe 0.0005806  0.0004971   1.168    0.243
AIC(reg2) # -729.8922

# Is JAS vwc predicted by snowmelt date?
reg3 <- lm(soilVWCData$jasVWC20mean ~ climData.sub$meltdoy)
plot(climData.sub$meltdoy, soilVWCData$jasVWC20mean)
abline(reg3)
summary(reg3)
# Yep, a pretty strong positive correlation, but the intercept fails
#                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -0.0255317  0.0556508  -0.459    0.646    
# climData.sub$meltdoy  0.0013035  0.0002379   5.478 5.28e-08 ***
# Multiple R-squared: 0.0255,	Adjusted R-squared: 0.02465 
# F-statistic: 30.01 on 1 and 1147 DF,  p-value: 5.283e-08 
AIC(reg3) # -758.2007

# Is JAS vwc predicted by JAS precip?
reg4 <- lm(soilVWCData$jasVWC20mean ~ climData.sub$JASprecip)
plot(climData.sub$JASprecip, soilVWCData$jasVWC20mean)
abline(reg4)
summary(reg4)
# Yeah, totally, not a great fit though
#                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.174535   0.009256   18.86   <2e-16 ***
# climData.sub$JASprecip 0.024204   0.001869   12.95   <2e-16 ***
# Multiple R-squared: 0.1268,	Adjusted R-squared: 0.126 
# F-statistic: 167.7 on 1 and 1155 DF,  p-value: < 2.2e-16
AIC(reg4) # -892.509

# Is JAS vwc predicted by jasMAT?
reg5 <- lm(soilVWCData$jasVWC20mean ~ climData.sub$jasMAT)
plot(climData.sub$jasMAT, soilVWCData$jasVWC20mean)
abline(reg5)
summary(reg5)
# Yes, but still a pretty crappy fit
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.6417030  0.0305979   20.97   <2e-16 ***
# jasMAT      -0.0089487  0.0007153  -12.51   <2e-16 ***
# Multiple R-squared: 0.1421,	Adjusted R-squared: 0.1412 
# F-statistic: 156.5 on 1 and 945 DF,  p-value: < 2.2e-16 
AIC(reg5) # -832.8904

# Is JAS vwc predicted by onsetdoy?
reg6 <- lm(soilVWCData$jasVWC20mean ~ climData.sub$onsetdoy)
plot(climData.sub$onsetdoy, soilVWCData$jasVWC20mean)
abline(reg6)
summary(reg6)
# Yes, but barely
#                         Estimate Std. Error t value Pr(>|t|)
# (Intercept)            0.2906356  0.0089778   32.37   <2e-16 ***
# climData.sub$onsetdoy -0.0005395  0.0003137   -1.72   0.0858 .
AIC(reg6) # -731.4851

# One other thing - is JAS vwc predicted by maat?
reg7 <- lm(soilVWCData$jasVWC20mean ~ climData.sub$maat)
plot(climData.sub$maat, soilVWCData$jasVWC20mean)
abline(reg7)
summary(reg7)
# Yes, very, but maybe this is because it correlates with so many other things?
#                         Estimate Std. Error t value Pr(>|t|)
# (Intercept)        0.406174   0.011500   35.32   <2e-16 ***
# climData.sub$maat -0.029746   0.002382  -12.49   <2e-16 ***
# Multiple R-squared: 0.1203,	Adjusted R-squared: 0.1195 
# F-statistic: 155.9 on 1 and 1140 DF,  p-value: < 2.2e-16 
AIC(reg7) # -872.2423

# Try some multiple regression, starting with the 2 best from simple reg:
regX <- lm(soilVWCData$jasVWC20mean ~ climData.sub$jasMAT+
	   climData.sub$JASprecip)
summary(regX)
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.4822059  0.0396124  12.173  < 2e-16 ***
# climData.sub$jasMAT    -0.0065009  0.0008069  -8.057 2.36e-15 ***
# climData.sub$JASprecip  0.0135980  0.0022019   6.176 9.79e-10 ***
AIC(regX) # -863.314
regX <- lm(soilVWCData$jasVWC20mean ~ climData.sub$jasMAT+
	   climData.sub$meltdoy)
AIC(regX) # -823.3053, meltdoy not significant
regX <- lm(soilVWCData$jasVWC20mean ~ climData.sub$jasMAT+
	   climData.sub$JASprecip+climData.sub$meltdoy)
AIC(regX) # -849.5018, meltdoy not significant
regX <- lm(soilVWCData$jasVWC20mean ~ climData.sub$jasMAT+
	   climData.sub$JASprecip+climData.sub$maxswe)
summary(regX)
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.5156288  0.0439008  11.745  < 2e-16 ***
# climData.sub$jasMAT    -0.0069945  0.0008408  -8.319 3.13e-16 ***
# climData.sub$JASprecip  0.0129720  0.0022355   5.803 8.95e-09 ***
# climData.sub$maxswe    -0.0005637  0.0005305  -1.063    0.288 
AIC(regX) # -850.6299
regX <- lm(soilVWCData$jasVWC20mean ~ climData.sub$jasMAT+
	   climData.sub$JASprecip+climData.sub$maxswe+climData.sub$meltdoy)
AIC(regX) # -849.3455, maxswe and meltdoy not significant
regX <- lm(soilVWCData$jasVWC20mean ~ climData.sub$jasMAT+
	   climData.sub$JASprecip+climData.sub$maxswe+climData.sub$meltdoy
	   +climData.sub$onsetdoy)
summary(regX)
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.3054816  0.1124245   2.717  0.00671 ** 
# climData.sub$jasMAT    -0.0068948  0.0009415  -7.323 5.26e-13 ***
# climData.sub$JASprecip  0.0133079  0.0022807   5.835 7.42e-09 ***
# climData.sub$maxswe    -0.0011189  0.0007070  -1.583  0.11384    
# climData.sub$meltdoy    0.0008246  0.0003982   2.071  0.03864 *  
# climData.sub$onsetdoy   0.0010524  0.0003512   2.996  0.00281 ** 
# Multiple R-squared: 0.1853,	Adjusted R-squared: 0.1809 
# F-statistic: 42.21 on 5 and 928 DF,  p-value: < 2.2e-16
AIC(regX) # -856.3377

# Getting to be a pretty complex model, but at least the coefficients for
# precip and jasMAT are staying the same. Can we drop anything?
drop1(regX, test='F')
#                        Df Sum of Sq    RSS     AIC F value    Pr(>F)    
# <none>                              21.537 -3508.9                      
# climData.sub$jasMAT     1   1.24452 22.781 -3458.4 53.6253 5.261e-13 ***
# climData.sub$JASprecip  1   0.79018 22.327 -3477.3 34.0481 7.424e-09 ***
# climData.sub$maxswe     1   0.05813 21.595 -3508.4  2.5048  0.113842    
# climData.sub$meltdoy    1   0.09954 21.636 -3506.6  4.2890  0.038637 *  
# climData.sub$onsetdoy   1   0.20835 21.745 -3501.9  8.9775  0.002806 **
# So... we could probably drop maxswe without much problem.
# I'm baffled as to why onsetday matters

# Now lets see if site effects matter
regSite <- lm(soilVWCData$jasVWC20mean ~ climData.sub$jasMAT+
	   climData.sub$JASprecip+climData.sub$maxswe+climData.sub$meltdoy
	   +climData.sub$onsetdoy+as.factor(climData.sub$siteClim))
summary(regSite)
#                                       Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                           0.3165641  0.0983303   3.219 0.001341 ** 
# climData.sub$jasMAT                 -0.0080112  0.0013135  -6.099 1.72e-09 ***
# climData.sub$JASprecip               0.0167360  0.0018068   9.263  < 2e-16 ***
# climData.sub$maxswe                  0.0019999  0.0005424   3.687 0.000244 ***
# climData.sub$meltdoy                 0.0003540  0.0002628   1.347 0.178383    
# climData.sub$onsetdoy               -0.0006406  0.0002163  -2.962 0.003158 ** 
# as.factor(climData.sub$siteClim)306  0.0196982  0.0620471   0.317 0.750976    
# as.factor(climData.sub$siteClim)310  0.1095817  0.0673734   1.626 0.104274
# ...
# Multiple R-squared: 0.8295,	Adjusted R-squared: 0.7841 
# F-statistic: 18.29 on 196 and 737 DF,  p-value: < 2.2e-16
AIC(regSite) # -1957.535
# Damn, so site is also a big player. Interestingly, adding site makes the
# effect of meltdoy disappear. Not sure how to interpret this yet.

drop1(regSite, test='F')
# climData.sub$jasMAT              1.722e-09 ***
# climData.sub$JASprecip           < 2.2e-16 ***
# climData.sub$maxswe              0.0002438 ***
# climData.sub$meltdoy             0.1783829    
# climData.sub$onsetdoy            0.0031578 ** 
# as.factor(climData.sub$siteClim) < 2.2e-16 ***
# This says we can drop meltdoy

# Lets check these results just using monthly data
# July
precip <- climData.sub$junPrecip+climData.sub$julPrecip
regJul <- lm(soilVWCData$julVWC20mean ~ climData.sub$julTairMean+
	   precip+climData.sub$maxswe+climData.sub$meltdoy
	   +climData.sub$onsetdoy)
summary(regJul)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -0.0054848  0.1284304  -0.043    0.966    
# climData.sub$julTairMean -0.0191162  0.0028315  -6.751 2.35e-11 ***
# precip                    0.0266734  0.0036092   7.390 2.86e-13 ***
# climData.sub$maxswe       0.0003451  0.0007949   0.434    0.664    
# climData.sub$meltdoy      0.0022837  0.0004565   5.002 6.58e-07 ***
# climData.sub$onsetdoy     0.0017993  0.0004194   4.290 1.94e-05 ***
# Multiple R-squared:   0.2,	Adjusted R-squared: 0.1965 
# F-statistic: 55.86 on 5 and 1117 DF,  p-value: < 2.2e-16
AIC(regJul) # -492.5251
# August
precip <- climData.sub$julPrecip+climData.sub$augPrecip
regAug <- lm(soilVWCData$augVWC20mean ~ climData.sub$augTairMean+
	   precip+climData.sub$maxswe+climData.sub$meltdoy
	   +climData.sub$onsetdoy)
summary(regAug)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               0.0226036  0.0971146   0.233 0.815997    
# climData.sub$augTairMean -0.0130461  0.0024422  -5.342 1.11e-07 ***
# precip                    0.0392830  0.0027524  14.272  < 2e-16 ***
# climData.sub$maxswe      -0.0021629  0.0006554  -3.300 0.000996 ***
# climData.sub$meltdoy      0.0013986  0.0003643   3.840 0.000130 ***
# climData.sub$onsetdoy     0.0013341  0.0003532   3.777 0.000167 ***
# Multiple R-squared: 0.2861,	Adjusted R-squared: 0.2829 
# F-statistic: 90.48 on 5 and 1129 DF,  p-value: < 2.2e-16 
AIC(regAug) # -890.8842

#September
precip <- climData.sub$augPrecip+climData.sub$sepPrecip
regSep <- lm(soilVWCData$sepVWC20mean ~ climData.sub$sepTairMean+
	   precip+climData.sub$maxswe+climData.sub$meltdoy
	   +climData.sub$onsetdoy)
summary(regSep)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               0.2242165  0.1022368   2.193   0.0285 *  
# climData.sub$sepTairMean -0.0159694  0.0024634  -6.483 1.44e-10 ***
# precip                    0.0429576  0.0031795  13.511  < 2e-16 ***
# climData.sub$maxswe      -0.0017799  0.0007022  -2.535   0.0114 *  
# climData.sub$meltdoy      0.0002806  0.0003917   0.716   0.4740    
# climData.sub$onsetdoy     0.0006133  0.0003484   1.761   0.0786 .  
# Multiple R-squared: 0.3115,	Adjusted R-squared: 0.3079 
# F-statistic: 86.43 on 5 and 955 DF,  p-value: < 2.2e-16
AIC(regSep) # -869.9365

# Seems a little inconclusive. What about doing this on individual sites?
climData.sub <- climData.sub[order(climData.sub$siteClim, 
				   climData.sub$yearClim),]
soilVWCData <- soilVWCData[order(soilVWCData$siteVWC, 
				   soilVWCData$yearVWC),]

library(nlme)

# Make some dataframes to analyze with linear models
reg.dat <- cbind(climData.sub[,c('siteClim','yearClim','maxswe','meltdoy',
			     'onsetdoy','JASprecip','junPrecip','julPrecip',
			     'augPrecip','sepPrecip','jasMAT','julTairMean',
			     'augTairMean','sepTairMean')],
	     soilVWCData[,c('jasVWC5mean','julVWC5mean','augVWC5mean',
			    'sepVWC5mean','jasVWC20mean','julVWC20mean',
			    'augVWC20mean','sepVWC20mean','jasVWC50mean',
			    'julVWC50mean','augVWC50mean','sepVWC50mean')])
# Create a list of lmLists for the simple linear regressions,
# but first change the precip columns
reg.dat$julPrecip <- (reg.dat$junPrecip + reg.dat$julPrecip)
reg.dat$augPrecip <- (reg.dat$julPrecip + reg.dat$augPrecip)
reg.dat$sepPrecip <- (reg.dat$augPrecip + reg.dat$sepPrecip)

simple.list <- list(lmList(jasVWC5mean~JASprecip | siteClim,data=reg.dat,
			   na.action=na.omit), 
		    lmList(jasVWC5mean~onsetdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(jasVWC5mean~meltdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(jasVWC5mean~maxswe | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(jasVWC5mean~jasMAT | as.factor(siteClim),
			   data=reg.dat,na.action=na.omit),
                    lmList(jasVWC20mean~JASprecip | siteClim,data=reg.dat,
			   na.action=na.omit), 
		    lmList(jasVWC20mean~onsetdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(jasVWC20mean~meltdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(jasVWC20mean~maxswe | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(jasVWC20mean~jasMAT | siteClim,data=reg.dat,
			   na.action=na.omit),
                    lmList(jasVWC50mean~JASprecip | siteClim,data=reg.dat,
			   na.action=na.omit), 
		    lmList(jasVWC50mean~onsetdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(jasVWC50mean~meltdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(jasVWC50mean~maxswe | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(jasVWC50mean~jasMAT | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(julVWC20mean~julPrecip | siteClim,data=reg.dat,
			   na.action=na.omit), 
		    lmList(julVWC20mean~onsetdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(julVWC20mean~meltdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(julVWC20mean~maxswe | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(julVWC20mean~julTairMean | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(augVWC20mean~augPrecip | siteClim,data=reg.dat,
			   na.action=na.omit), 
		    lmList(augVWC20mean~onsetdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(augVWC20mean~meltdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(augVWC20mean~maxswe | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(augVWC20mean~augTairMean | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(sepVWC20mean~sepPrecip | siteClim,data=reg.dat,
			   na.action=na.omit), 
		    lmList(sepVWC20mean~onsetdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(sepVWC20mean~meltdoy | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(sepVWC20mean~maxswe | siteClim,data=reg.dat,
			   na.action=na.omit),
		    lmList(sepVWC20mean~sepTairMean | siteClim,data=reg.dat,
			   na.action=na.omit))
# Each lmlist contains coeffs and p-values for each individual site regression
# Get them and put them in a table
simple.dat <- data.frame(num_lm=1:length(simple.list),noRegI=0,intAv=0,
			 noSig1=0,noRegB=0,beta1Av=0,noSig2=0)
for (i in 1:length(simple.list)) {
	simple <- simple.list[[i]]
	int1 <- summary(simple)$coefficients[,1,1]
	intStdE <- summary(simple)$coefficients[,2,1]
	intPvals <- summary(simple)$coefficients[,4,1]
        sigtest1 <- intPvals < 0.05
	beta1 <- summary(simple)$coefficients[,1,2]
	beta1StdE <- summary(simple)$coefficients[,2,2]
	beta1Pvals <- summary(simple)$coefficients[,4,2]
	sigtest2 <- beta1Pvals < 0.05
        simple.dat[i,'lmcall'] <- paste(summary(simple)$call[2])
	simple.dat[i,'noRegI'] <- sum(!is.na(intStdE))
	simple.dat[i,'intAv'] <- mean(int1[sigtest1],na.rm=TRUE)
	simple.dat[i,'noSig1'] <- sum(sigtest1,na.rm=TRUE)
	simple.dat[i,'noRegB'] <- sum(!is.na(beta1StdE))
	simple.dat[i,'beta1Av'] <- mean(beta1[sigtest2],na.rm=TRUE)
	simple.dat[i,'noSig2'] <- sum(sigtest2,na.rm=TRUE)
}

# Now lets do some multiple regressions
multi.list <- list(lmList(jasVWC5mean~JASprecip+jasMAT+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(jasVWC20mean~JASprecip+jasMAT+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(jasVWC50mean~JASprecip+jasMAT+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(julVWC5mean~julPrecip+julTairMean+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(augVWC5mean~augPrecip+augTairMean+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(sepVWC5mean~sepPrecip+sepTairMean+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(julVWC20mean~julPrecip+julTairMean+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(augVWC20mean~augPrecip+augTairMean+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(sepVWC20mean~sepPrecip+sepTairMean+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(julVWC50mean~julPrecip+julTairMean+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(augVWC50mean~augPrecip+augTairMean+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit),
		   lmList(sepVWC50mean~sepPrecip+sepTairMean+onsetdoy
			  | siteClim,data=reg.dat,na.action=na.omit))

# Each lmlist contains coeffs and p-values for each individual site regression
# Get them and put them in a table
multi.dat <- data.frame(num_lm=1:length(multi.list),noRegI=0,intAv=0,noSigI=0,
			noReg1=0,beta1Av=0,noSig1=0,
			noReg2=0,beta2Av=0,noSig2=0,
			noReg3=0,beta3Av=0,noSig3=0)
			
for (i in 1:length(multi.list)) {
	multi <- multi.list[[i]]
	multi.dat[i,'lmcall'] <- paste(summary(multi)$call[2])
	colindex <- seq(2, 11, by=3)
	for(j in 1:4) {
		coeffs <- summary(multi)$coefficients[,1,j]
		StdErrs <- summary(multi)$coefficients[,2,j]
		Pvals <- summary(multi)$coefficients[,4,j]
		sigtest <- Pvals < 0.05
		multi.dat[i,colindex[j]] <- sum(!is.na(StdErrs))
		multi.dat[i,colindex[j]+1] <- mean(coeffs[sigtest],na.rm=TRUE)
		multi.dat[i,colindex[j]+2] <- sum(sigtest,na.rm=TRUE)
	}
}
# onsetdoy and meltdoy are similarly significant, however there are generally
# only 60-70 sites with valid regressions and 1-6 sites with significant
# coefficients for any of the predictor variables.
#
# A major problem here is that there are only a few datapoints for each site,
# so there is little power to do multiple regression. There are only a few
# sites with significant regressions
#
# Another major problem is that lmList doesn't seem to give the right linear
# model results, or at least they don't compare to what I get with plain old
# lm() in a loop. For this reason I am moving this job to a new script and will
# use those results rather than simple.dat




plot(density(regdat$scdSlope, na.rm=TRUE))
plot(density(regdat$scdPval, na.rm=TRUE))
plot(density(regdat$lmSlope1, na.rm=TRUE))
plot(density(regdat$modelPval, na.rm=TRUE))
