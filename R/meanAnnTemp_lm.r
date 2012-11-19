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

# Is MAST better predicted by a 2 term model of MAT and totaldaysSC?
reg4 <- lm(soilTData$mast20cm ~ climData.sub$maat + climData.sub$totaldaysSC)
summary(reg4)
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               7.475158   0.301050   24.83   <2e-16 ***
# climData.sub$maat         0.336754   0.021595   15.59   <2e-16 ***
# climData.sub$totaldaysSC -0.019561   0.001166  -16.77   <2e-16 ***
# Multiple R-squared: 0.5745,	Adjusted R-squared: 0.5738 
# F-statistic: 782.6 on 2 and 1159 DF,  p-value: < 2.2e-16 
#
# So, both are highly significant here and have opposing slopes. The model is
# significant too, and I checked all depths and got the same result. r2 is a
# bit higher for 50cm depth



# What about the difference (MAST-MAT)?
diff <- (soilTData$mast20cm - climData.sub$maat)
reg5 <-lm (diff ~ climData.sub$totaldaysSC)
plot(climData.sub$totaldaysSC, diff)
abline(reg5)
summary(reg5)
# Interesting, totaldaysSC is highly significant, intercept isn't. Slope is
# weak and has high uncertainty. Maybe because it goes to 0?
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -0.072132   0.234126  -0.308 0.758068    
# climData.sub$totaldaysSC  0.004271   0.001172   3.643 0.000281 ***

reg6 <-lm (diff ~ climData.sub$maat)
plot(climData.sub$maat, diff)
abline(reg6)
summary(reg6)
# This looks good, maat is highly significant
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.48089    0.08520   29.12   <2e-16 ***
# climData.sub$maat -0.39231    0.01767  -22.20   <2e-16 ***

# How about a 2 term model with MAT and totaldays SC
reg7 <- lm(diff ~ climData.sub$totaldaysSC + climData.sub$maat)
summary(reg7)
# (Intercept)               7.475158   0.301050   24.83   <2e-16 ***
# climData.sub$totaldaysSC -0.019561   0.001166  -16.77   <2e-16 ***
# climData.sub$maat        -0.663246   0.021595  -30.71   <2e-16 ***
# Multiple R-squared: 0.4549,	Adjusted R-squared: 0.454 
# F-statistic: 483.7 on 2 and 1159 DF,  p-value: < 2.2e-16
# This looks good, both terms are highly significant and the model is
# significant too. Checked other depths and they look good also, 50cm model
# may fit a bit better again
# However, maat and totalSCdays are highly correlated (see reg1). now what?

layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(reg7) #(I don't know what these diagnostics mean!)

# Look at the interaction
reg8 <-lm (diff ~ climData.sub$totaldaysSC + climData.sub$maat + 
	   climData.sub$totaldaysSC*climData.sub$maat)
summary(reg8)
# Interaction term doesn't look very significant (p=0.075)


# What about a complicated model with meltdoy and onsetdoy?
reg9 <- lm(diff ~ climData.sub$totaldaysSC + climData.sub$maat + 
	  climData.sub$meltdoy + climData.sub$onsetdoy)
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               8.573531   0.475752  18.021  < 2e-16 ***
# climData.sub$totaldaysSC -0.013619   0.002105  -6.469 1.46e-10 ***
# climData.sub$maat        -0.697586   0.019862 -35.121  < 2e-16 ***
# climData.sub$meltdoy     -0.009992   0.002896  -3.451  0.00058 ***
# climData.sub$onsetdoy     0.005221   0.002369   2.204  0.02774 *  
# Residual standard error: 1.017 on 1146 degrees of freedom
#   (439 observations deleted due to missingness)
# Multiple R-squared: 0.5321,	Adjusted R-squared: 0.5305 
# F-statistic: 325.8 on 4 and 1146 DF,  p-value: < 2.2e-16 
#
# Crap - all are significant, but I did correlations between them all and
# they are all sifnificantly correlated (even onsetdoy and meltdoy!)
# Maybe we can lump all the snow parameters into one?



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
#
# Actually I think this means that the predictor is a combination of MAT and
# total snowcovered days.

# Try a more complicated PCR
x <- cbind(climData.sub$totaldaysSC, climData.sub$maat, climData.sub$meltdoy,
	   climData.sub$onsetdoy, climData.sub$maxswe)
natest <- (is.na(x[,1]) | is.na(x[,2]) | is.na(x[,3]) | is.na(x[,4]) | 
			       is.na(x[,1]))
x <- x[!natest,]
y <- as.matrix(diff)[!natest,]
xxxxxy.pcr <- pcr(y ~ x, ncomp = 5, method = "svdpc", validation='CV')
summary(xxxxxy.pcr)
# VALIDATION: RMSEP
# Cross-validated using 10 random segments.
#        (Intercept)  1 comps  2 comps  3 comps  4 comps  5 comps
# CV           1.485    1.477    1.476    1.455    1.436    1.022
# adjCV        1.485    1.476    1.476    1.455    1.435    1.022
# 
# TRAINING: % variance explained
#    1 comps  2 comps  3 comps  4 comps  5 comps
# X   86.220   94.334   97.722   99.911   100.00
# y    1.216    1.524    4.683    7.252    53.32

# Seems a little inconclusive. What about doing this on all sites?
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
