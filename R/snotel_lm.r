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

# We're interested in whether MAST and MAST-MAT are related to mean
# annual air temperature or duration of snowcover (todaldaysSC). First
# lets check whether these 2 predictors are correlated.
reg1 <- lm(climData$totaldaysSC ~ climData$maat)
plot(climData$maat, climData$totaldaysSC)
abline(reg1)
summary(reg1)
# Yes, totaldaysSC and maat are highly correlated
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   248.7383     0.7340   338.9   <2e-16 ***
# climData$maat -12.8645     0.1835   -70.1   <2e-16 ***

# Is MAST predicted by total snowcovered days?
reg2 <- lm(soilTData$mast20cm ~ climData.sub$totaldaysSC)
plot(climData.sub$totaldaysSC, soilTData$mast20cm)
abline(reg2)
summary(reg2)
# Yep, a pretty strong negative correlation
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              11.30053    0.18971   59.57   <2e-16 ***
# climData.sub$totaldaysSC -0.03163    0.00095  -33.30   <2e-16 ***

# Is MAST predicted by MAT?
reg3 <- lm(soilTData$mast20cm ~ climData.sub$maat)
plot(climData.sub$maat, soilTData$mast20cm)
abline(reg3)
summary(reg3)
# Yep, a pretty strong positive correlation
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.48089    0.08520   29.12   <2e-16 ***
# climData.sub$maat  0.60769    0.01767   34.38   <2e-16 ***                  

# What about the difference (MAST-MAT)?
diff <- (soilTData$mast20cm - climData.sub$maat)
reg4 <-lm (diff ~ climData.sub$totaldaysSC)
plot(climData.sub$totaldaysSC, diff)
abline(reg4)
summary(reg4)
# Interesting, totaldaysSC is highly significant, intercept isn't. Slope is
# weak and has high uncertainty.
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -0.072132   0.234126  -0.308 0.758068    
# climData.sub$totaldaysSC  0.004271   0.001172   3.643 0.000281 ***

reg5 <-lm (diff ~ climData.sub$maat)
plot(climData.sub$maat, diff)
abline(reg5)
summary(reg5)
# This looks good, maat is highly significant
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.48089    0.08520   29.12   <2e-16 ***
# climData.sub$maat -0.39231    0.01767  -22.20   <2e-16 ***

# How about a 2 term model with MAT and totaldays SC
reg6 <- lm(diff ~ climData.sub$totaldaysSC + climData.sub$maat)
summary(reg6)
# This looks good, both terms are highly significant
# (Intercept)               7.475158   0.301050   24.83   <2e-16 ***
# climData.sub$totaldaysSC -0.019561   0.001166  -16.77   <2e-16 ***
# climData.sub$maat        -0.663246   0.021595  -30.71   <2e-16 ***
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(reg6) #(I don't know what these diagnostics mean!)
# However, maat and totalSCdays are highly correlated (see reg1). now what?

# Look at the interaction
reg7 <-lm (diff ~ climData.sub$totaldaysSC + climData.sub$maat + 
	   climData.sub$totaldaysSC*climData.sub$maat)
summary(reg7)
# Interaction term doesn't look very significant (p=0.075)

# What about pcr - principle components regression
library(pls)
x <- cbind(climData.sub$totaldaysSC, climData.sub$maat)
natest <- is.na(x[,1]) | is.na(x[,2])
x <- x[!natest,]
y <- as.matrix(diff)[!natest,]
xxy.pcr <- pcr(y ~ x, ncomp = 2, method = "svdpc", validation='CV')
summary(xxy.pcr)
# VALIDATION: RMSEP
# Cross-validated using 10 random segments.
#        (Intercept)  1 comps  2 comps
# CV           1.484    1.477    1.098
# adjCV        1.484    1.476    1.097
# 
# TRAINING: % variance explained
#    1 comps  2 comps
# X   99.838   100.00
# y    1.152    45.49
# Not quite sure what the validation means but it looks like one predictor,
# and it could be either one (I switched order in X and got same result),
# explains 99.38% of the variance?

# Well that was a bit inconclusive, so what about trying this on individual
#sites?

climData.sub <- climData.sub[order(climData.sub$siteClim, 
				   climData.sub$yearClim),]
soilTData <- soilTData[order(soilTData$siteTsoil, 
				   soilTData$yearTsoil),]
soilVWCData <- soilVWCData[order(soilVWCData$siteVWC, 
				   soilVWCData$yearVWC),]

sites = unique(soilTData$siteTsoil)
regdat <- data.frame(site=sites,nyrs=as.vector(table(soilTData$siteTsoil)),
		     lmInt=0,lmSlope1=0,lmSlope2=0,PvalInt=0,Pval1=0,Pval2=0,
		     lmR2=0,modelPval=0,comp1=0,scdInt=0,scdSlope=0,scdPval=0,
		     scdR2=0)
for (i in 1:length(sites)) {
    tmpC <- climData.sub[climData.sub$siteClim==sites[i],]
    tmpT <- soilTData[soilTData$siteTsoil==sites[i],]
    tmpDiff <- tmpT$mast20cm - tmpC$maat
    tmp2X <- cbind(tmpC$totaldaysSC, tmpC$maat)
    natest <- is.na(tmp2X[,1]) | is.na(tmp2X[,2]) | is.na(tmpDiff)
    tmp2X <- tmp2X[!natest,]
    tmpy <- as.matrix(tmpDiff)[!natest,]
    if (length(tmpy) > 4){ 
        #tmpC <- climData.sub[climData.sub$siteClim==sites[i],]
        #tmpT <- soilTData[soilTData$siteTsoil==sites[i],]
        #tmpDiff <- tmpT$mast20cm - tmpC$maat
        #tmp2X <- cbind(tmpC$totaldaysSC, tmpC$maat)
        #natest <- is.na(tmp2X[,1]) | is.na(tmp2X[,2])
        #tmp2X <- tmp2X[!natest,]
        #tmpy <- as.matrix(tmpDiff)[!natest,]
        # Do the pcr for the site
        diff.pcr <- pcr(tmpy ~ tmp2X, ncomp=2, method='svdpc')
        regdat[i,'comp1'] <- (diff.pcr$Xvar/diff.pcr$Xtotvar)[1]
        # Do the multiple regression and get F, pvals, and r2 
        diff.lm <- lm(tmpDiff ~ tmpC$totaldaysSC + tmpC$maat)
        fs <- summary(diff.lm)$fstatistic #Get F, df1, df2
        regdat[i,'modelPval'] <- pf(fs[1], fs[2], fs[3], lower.tail=FALSE)
	# Do the simple regression (on totaldaysSC and get stats
	simp.lm <- lm(tmpDiff ~ tmpC$totaldaysSC)
        fs <- summary(simp.lm)$fstatistic #Get F, df1, df2
        regdat[i,'scdPval'] <- pf(fs[1], fs[2], fs[3], lower.tail=FALSE)
        # Add stats to regdat
        regdat[i,'lmInt'] <- diff.lm$coef[1]
        regdat[i,'lmSlope1'] <- diff.lm$coef[2]
        regdat[i,'lmSlope2'] <- diff.lm$coef[3]
        regdat[i,'PvalInt'] <- summary(diff.lm)$coef[10]
        regdat[i,'Pval1'] <- summary(diff.lm)$coef[11]
        regdat[i,'Pval2'] <- summary(diff.lm)$coef[12]
        regdat[i,'lmR2'] <- summary(diff.lm)$r.squared
	regdat[i,'scdInt'] <- simp.lm$coef[1]
	regdat[i,'scdSlope'] <- simp.lm$coef[2]
        regdat[i,'scdR2'] <- summary(simp.lm)$r.squared

    } else {
	regdat[i,-(1:2)] <- NA
    }
}


plot(density(regdat$scdSlope, na.rm=TRUE))
plot(density(regdat$scdPval, na.rm=TRUE))
plot(density(regdat$lmSlope1, na.rm=TRUE))
plot(density(regdat$modelPval, na.rm=TRUE))
