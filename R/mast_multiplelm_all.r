source('getdata.r')

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
# Multiple R-squared: 0.4838,	Adjusted R-squared: 0.4834 
# F-statistic:  1109 on 1 and 1183 DF,  p-value: < 2.2e-16
AIC(reg2) # 3803.425

# Is MAST predicted by MAT?
reg3 <- lm(soilTData$mast20cm ~ climData.sub$maat)
plot(climData.sub$maat, soilTData$mast20cm)
abline(reg3)
summary(reg3)
# Yep, a pretty strong positive correlation
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.48089    0.08520   29.12   <2e-16 ***
# climData.sub$maat  0.60769    0.01767   34.38   <2e-16 ***                  
# Multiple R-squared: 0.5015,	Adjusted R-squared: 0.5011 
# F-statistic:  1182 on 1 and 1175 DF,  p-value: < 2.2e-1
AIC(reg3) # 3847.411

# Is MAST better predicted by a 2 term model of MAT and totaldaysSC?
reg4 <- lm(soilTData$mast20cm ~ climData.sub$maat+climData.sub$totaldaysSC)
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
AIC(reg4) # 3515.17

# What if we toss in some other parameters
regX <- lm(soilTData$mast20cm ~ climData.sub$maat+climData.sub$totaldaysSC
	   +climData.sub$meltdoy)
AIC(regX) # 3510.032
regX <- lm(soilTData$mast20cm ~ climData.sub$maat+climData.sub$totaldaysSC
	   +climData.sub$onsetdoy)
AIC(regX) # 3513.686
regX <- lm(soilTData$mast20cm ~ climData.sub$maat+climData.sub$totaldaysSC
	   +climData.sub$maxswe)
AIC(regX) # 3515.88
regX <- lm(soilTData$mast20cm ~ climData.sub$maat+climData.sub$totaldaysSC
	   +climData.sub$meltdoy+climData.sub$onsetdoy)
AIC(regX) # 3506.315
regX <- lm(soilTData$mast20cm ~ climData.sub$maat+climData.sub$totaldaysSC
	   +climData.sub$meltdoy+climData.sub$maxswe)
AIC(regX) # 3506.785
regX <- lm(soilTData$mast20cm ~ climData.sub$maat+climData.sub$totaldaysSC
	   +climData.sub$meltdoy+climData.sub$onsetdoy+climData.sub$maxswe)
AIC(regX) # 3504.815
regX <- lm(soilTData$mast20cm ~ climData.sub$maat+
	   +climData.sub$meltdoy+climData.sub$onsetdoy)
AIC(regX) # 3535.027


# This is getting pretty complex, but all these seem to influence the fit of the model. Try drop1:
drop1(regX, test='F')
#                          Df Sum of Sq    RSS    AIC  F value    Pr(>F)    
# <none>                                1372.2 205.20                       
# climData.sub$maat         1   244.838 1617.0 393.98 206.2636 < 2.2e-16 ***
# climData.sub$totaldaysSC  1    40.907 1413.1 237.34  34.4620 5.673e-09 ***
# climData.sub$meltdoy      1    14.687 1386.9 215.57  12.3728 0.0004525 ***
# climData.sub$onsetdoy     1     4.697 1376.9 207.17   3.9566 0.0469237 *  
# climData.sub$maxswe       1     4.140 1376.3 206.70   3.4876 0.0620841 .
# So... we could probably drop maxswe without much problem

# Now lets see if site effects matter
regSite <- lm(soilTData$mast20cm ~ climData.sub$maat+climData.sub$totaldaysSC
	   +climData.sub$meltdoy+climData.sub$onsetdoy
	   +as.factor(climData.sub$siteClim))
# Multiple R-squared: 0.9248,	Adjusted R-squared: 0.9072 
# F-statistic: 52.58 on 220 and 941 DF,  p-value: < 2.2e-16
AIC(regSite) # 1937.815
# Damn, so site is also a big player. Interestingly, adding site makes the
# effect of totaldaysSC disappear. Shows up in drop1 too. Not sure how to
# interpret this yet.


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
AIC(reg5) # 4205.112

reg6 <-lm (diff ~ climData.sub$maat)
plot(climData.sub$maat, diff)
abline(reg6)
summary(reg6)
# This looks good, maat is highly significant
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.48089    0.08520   29.12   <2e-16 ***
# climData.sub$maat -0.39231    0.01767  -22.20   <2e-16 ***
AIC(reg6) # 3847.411

# How about a 2 term model with MAT and totaldays SC
reg7 <- lm(diff ~ climData.sub$totaldaysSC + climData.sub$maat)
summary(reg7)
# (Intercept)               7.475158   0.301050   24.83   <2e-16 ***
# climData.sub$totaldaysSC -0.019561   0.001166  -16.77   <2e-16 ***
# climData.sub$maat        -0.663246   0.021595  -30.71   <2e-16 ***
# Multiple R-squared: 0.4549,	Adjusted R-squared: 0.454 
# F-statistic: 483.7 on 2 and 1159 DF,  p-value: < 2.2e-16
AIC(reg7) # 3515.17
# This looks good, both terms are highly significant and the model is
# significant too. Checked other depths and they look good also, 50cm model
# may fit a bit better again
# Try a few others
regX <- lm(diff ~ climData.sub$totaldaysSC+climData.sub$maat
	   +climData.sub$meltdoy)
AIC(regX) #3510.32
regX <- lm(diff ~ climData.sub$totaldaysSC+climData.sub$maat
	   +climData.sub$meltdoy+climData.sub$onsetdoy)
AIC(regX) #3506.315
regX <- lm(diff ~ climData.sub$totaldaysSC+climData.sub$maat
	   +climData.sub$meltdoy+climData.sub$onsetdoy+climData.sub$maxswe)
AIC(regX) #3504.815

# So the story is about the same for diff - best fit is 4 or 5 parameters
# However, maat and totalSCdays are highly correlated (see reg1). now what?


# Look at the interaction
reg8 <-lm (diff ~ climData.sub$totaldaysSC + climData.sub$maat + 
	   climData.sub$totaldaysSC*climData.sub$maat)
summary(reg8)
# Interaction term doesn't look very significant (p=0.075)


# What about a complicated model with meltdoy and onsetdoy (like with MAST)?
reg9 <- lm(diff ~ climData.sub$totaldaysSC + climData.sub$maat + 
	  climData.sub$meltdoy + climData.sub$onsetdoy)
summary(reg9)
# NOTE: this is for diff at 50cm - 20cm looks similar
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
AIC(reg9) # 3312.676
# Crap - all are significant, but I did correlations between them all and
# they are all sifnificantly correlated (even onsetdoy and meltdoy!)
# Maybe we can lump all the snow parameters into one?

# Now lets see if site effects matter
regSiteD <- lm(diff ~ climData.sub$maat+climData.sub$totaldaysSC
	   +climData.sub$meltdoy+climData.sub$onsetdoy
	   +as.factor(climData.sub$siteClim))
# Multiple R-squared: 0.9036,	Adjusted R-squared: 0.8811 
# F-statistic:  40.1 on 220 and 941 DF,  p-value: < 2.2e-16
AIC(regSiteD) # 1937.815
# OK - site is also a big player with diff too (same AIC value, wow!).
# totaldaysSC dropped out here too. 

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
	   climData.sub$onsetdoy, climData.sub$siteClim)
natest <- (is.na(x[,1]) | is.na(x[,2]) | is.na(x[,3]) | is.na(x[,4]) | 
			       is.na(x[,5]))
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

# Seems a little inconclusive. What if we do this for snowcovered and snowfree
# soil temperature
# Start with a complicated model with meltdoy and onsetdoy (like with MAST)?
reg10 <- lm(soilTData$snowcovTs20mean ~ climData.sub$scovmat+
	    climData.sub$onsetdoy+climData.sub$novSWEmean)
summary(reg10)
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.821839   0.056801  14.469  < 2e-16 ***
# climData.sub$scovmat    0.155109   0.014835  10.456  < 2e-16 ***
# climData.sub$onsetdoy   0.002160   0.001407   1.535    0.125    
# climData.sub$novSWEmean 0.100237   0.014237   7.041 3.26e-12 ***
AIC(reg10) # 2765.598

reg11 <- lm(soilTData$snowcovTs20mean ~ climData.sub$preonsetTair+
	    climData.sub$onsetdoy+climData.sub$novSWEmean)
summary(reg11)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               -0.289589   0.092167  -3.142  0.00172 ** 
# climData.sub$preonsetTair  0.110551   0.010619  10.411  < 2e-16 ***
# climData.sub$onsetdoy      0.014630   0.001756   8.333  < 2e-16 ***
# climData.sub$novSWEmean    0.105528   0.014301   7.379 3.03e-13 ***
AIC(reg11) # 2754.598

reg12 <- lm(soilTData$snowcovTs20mean ~ climData.sub$preonsetTair+
	    climData.sub$onsetdoy+climData.sub$novSWEmean+climData.sub$scovmat)
summary(reg12)
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               0.194799   0.133551   1.459    0.145    
# climData.sub$preonsetTair 0.069420   0.013391   5.184 2.56e-07 ***
# climData.sub$onsetdoy     0.009503   0.002018   4.710 2.78e-06 ***
# climData.sub$novSWEmean   0.102999   0.014171   7.268 6.70e-13 ***
# climData.sub$scovmat      0.094754   0.018832   5.031 5.64e-07 ***
AIC(reg12) # 2703.864

reg13 <- lm(soilTData$snowcovTs20mean ~ climData.sub$preonsetTair+
	    climData.sub$onsetdoy+climData.sub$novSWEmean+climData.sub$scovmat
	    +climData.sub$maxswe)
summary(reg13)
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               0.063985   0.142650   0.449   0.6538    
# climData.sub$preonsetTair 0.068491   0.013363   5.125 3.48e-07 ***
# climData.sub$onsetdoy     0.010169   0.002029   5.011 6.27e-07 ***
# climData.sub$novSWEmean   0.080347   0.016667   4.821 1.62e-06 ***
# climData.sub$scovmat      0.083240   0.019315   4.310 1.78e-05 ***
# climData.sub$maxswe       0.007048   0.002747   2.566   0.0104 *  

AIC(reg13) # 2699.265

reg14 <- lm(soilTData$snowcovTs20mean ~ climData.sub$preonsetTair+
	    climData.sub$onsetdoy+climData.sub$decSWEmean+climData.sub$scovmat
	    +climData.sub$maxswe)
summary(reg14)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.044769   0.144313   0.310    0.756    
# climData.sub$preonsetTair  0.069637   0.013468   5.171 2.75e-07 ***
# climData.sub$onsetdoy      0.009370   0.002025   4.628 4.11e-06 ***
# climData.sub$decSWEmean    0.079630   0.013495   5.901 4.76e-09 ***
# climData.sub$scovmat       0.078999   0.019673   4.016 6.32e-05 ***
# climData.sub$maxswe       -0.004424   0.003918  -1.129    0.259
AIC(reg14) # 2699.669

reg15 <- lm(soilTData$snowcovTs20mean ~ climData.sub$preonsetTair+
	    climData.sub$onsetdoy+climData.sub$decSWEmean+climData.sub$scovmat)
summary(reg15)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               -0.008190   0.136495  -0.060 0.952167    
# climData.sub$preonsetTair  0.070529   0.013446   5.245 1.86e-07 ***
# climData.sub$onsetdoy      0.009816   0.001986   4.943 8.82e-07 ***
# climData.sub$decSWEmean    0.067384   0.008032   8.390  < 2e-16 ***
# climData.sub$scovmat       0.073321   0.019022   3.855 0.000122 ***

AIC(reg15) # 2698.95

# This is a good model including site effects
regSite <- lm(soilTData$snowcovTs20mean ~ climData.sub$preonsetTair+
	    climData.sub$decSWEmean+climData.sub$scovmat+climData.sub$onsetdoy
	    +as.factor(climData.sub$siteClim))
summary(regSite)
#                                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                     -0.3554677  0.3016446  -1.178 0.238925    
# climData.sub$preonsetTair        0.0484763  0.0098238   4.935 9.51e-07 ***
# climData.sub$decSWEmean          0.0867157  0.0075848  11.433  < 2e-16 ***
# climData.sub$scovmat            -0.0581096  0.0234971  -2.473 0.013574 *  
# climData.sub$onsetdoy            0.0029860  0.0016407   1.820 0.069092 .
AIC(regSite) #1825.254

# This is the best model including site effects - dropped another parameter
regSite <- lm(soilTData$snowcovTs20mean ~ climData.sub$preonsetTair+
	    climData.sub$decSWEmean+climData.sub$scovmat
	    +as.factor(climData.sub$siteClim))
summary(regSite)
#                                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                     -0.2381520  0.2950412  -0.807 0.419767    
# climData.sub$preonsetTair        0.0361972  0.0071491   5.063 4.97e-07 ***
# climData.sub$decSWEmean          0.0901358  0.0073574  12.251  < 2e-16 ***
# climData.sub$scovmat            -0.0540554  0.0234202  -2.308 0.021213 * 
AIC(regSite) # 1827.348


# However, there are some interactions to pay attention to also
# First - between winter temp and snowcover
regSite <- lm(soilTData$snowcovTs20mean ~ as.factor(climData.sub$siteClim)
	      +climData.sub$preonsetTair+
	    climData.sub$decSWEmean*climData.sub$scovmat)
summary(regSite)
AIC(regSite) # 1819.243

# Add early season interaction
regSite <- lm(soilTData$snowcovTs20mean ~ as.factor(climData.sub$siteClim)
	      +climData.sub$preonsetTair*climData.sub$decSWEmean + 
	      climData.sub$maxswe*climData.sub$scovmat)
summary(regSite)
AIC(regSite) # 1786.243

# Remove later season interaction
regSite <- lm(soilTData$snowcovTs20mean ~ as.factor(climData.sub$siteClim)
	      +climData.sub$preonsetTair*climData.sub$decSWEmean + 
	      climData.sub$maxswe + climData.sub$scovmat)
summary(regSite)
#                               Estimate Std. Error t value
# (Intercept)                 -0.3732476  0.2940521  -1.26   0.2046
# climData.sub$preonsetTair    0.0765750  0.0105201   7.279  7.15e-13 ***
# climData.sub$decSWEmean      0.1701198  0.0153027  11.117  < 2e-16  ***
# climData.sub$maxswe         -0.0098241  0.0031394  -3.129  0.001806 **
# climData.sub$scovmat        -0.0496339  0.0230384  -2.154  0.031465 *
# preonsetTair:decSWEmean     -0.0106136  0.0021401  -4.959  8.40e-07 ***
AIC(regSite) # 1788.063

regSite <- lm(soilTData$snowcovTs20mean ~
	      climData.sub$preonsetTair*climData.sub$decSWEmean + 
	      climData.sub$maxswe + climData.sub$scovmat)
summary(regSite)



reg16 <- lm(soilTData$snowcovTs20mean ~ climData.sub$preonsetTair+
	    climData.sub$onsetdoy+climData.sub$decSWEmean)
summary(reg16)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               -0.389774   0.092774  -4.201 2.86e-05 ***
# climData.sub$preonsetTair  0.102205   0.010608   9.635  < 2e-16 ***
# climData.sub$onsetdoy      0.013776   0.001705   8.079 1.63e-15 ***
# climData.sub$decSWEmean    0.069918   0.008057   8.678  < 2e-16 ***
AIC(reg16) # 2739.421

# Now lets do snowfree
reg17 <- lm(soilTData$snowfreeTs20mean ~ climData.sub$freemat+
	    climData.sub$meltdoy)
summary(reg17)
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           9.944001   0.628199   15.83   <2e-16 ***
# climData.sub$freemat  0.560325   0.031435   17.82   <2e-16 ***
# climData.sub$meltdoy -0.027482   0.002388  -11.51   <2e-16 ***
AIC(reg17) # 4532.077

reg18 <- lm(soilTData$snowfreeTs20mean ~ climData.sub$freemat+
	    climData.sub$meltdoy+ climData.sub$maxswe)
summary(reg18)
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          10.640595   0.789727  13.474   <2e-16 ***
# climData.sub$freemat  0.548845   0.032396  16.942   <2e-16 ***
# climData.sub$meltdoy -0.030708   0.003258  -9.425   <2e-16 ***
# climData.sub$maxswe   0.009949   0.006841   1.454    0.146
AIC(reg18) # 4531.956

reg19 <- lm(soilTData$snowfreeTs20mean ~ climData.sub$freemat+
	    climData.sub$meltdoy+soilTData$snowcovTs20mean)
summary(reg19)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               10.894673   0.582347   18.71   <2e-16 ***
# climData.sub$freemat       0.431111   0.030320   14.22   <2e-16 ***
# climData.sub$meltdoy      -0.027951   0.002199  -12.71   <2e-16 ***
# soilTData$snowcovTs20mean  0.846378   0.058716   14.41   <2e-16 ***
AIC(reg19) # 4331.002

# Add site effects to this
reg20 <- lm(soilTData$snowfreeTs20mean ~ climData.sub$freemat+
	    climData.sub$meltdoy+soilTData$snowcovTs20mean+
	    as.factor(climData.sub$siteClim))
summary(reg20)
#                                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                       5.080677   0.463234  10.968  < 2e-16 ***
# climData.sub$freemat              0.612046   0.021581  28.361  < 2e-16 ***
# climData.sub$meltdoy             -0.007741   0.001376  -5.627 2.43e-08 ***
# soilTData$snowcovTs20mean         0.719859   0.033096  21.751  < 2e-16 ***
AIC(reg20) # 2024.006


reg21 <- lm(soilTData$snowfreeTs20mean ~ climData.sub$freemat+
	    climData.sub$meltdoy+soilTData$snowcovTs20mean+
	    climData.sub$totaldaysSC)
summary(reg21)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                9.883580   0.707140  13.977  < 2e-16 ***
# climData.sub$freemat       0.431290   0.030250  14.257  < 2e-16 ***
# climData.sub$meltdoy      -0.018094   0.004501  -4.020  6.2e-05 ***
# soilTData$snowcovTs20mean  0.843083   0.058595  14.388  < 2e-16 ***
# climData.sub$totaldaysSC  -0.006567   0.002618  -2.508   0.0123 *
AIC(reg21) # 4326.7

