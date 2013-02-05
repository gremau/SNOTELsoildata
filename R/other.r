# Other interesting regression results

source('getdata.r')
reg1 <- lm(soilTData$snowcovTs20mean~climData.sub$maxswe)
plot(climData.sub$maxswe, soilTData$snowcovTs20mean)
abline(reg1)
summary(reg1) 
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.455721   0.049266   9.250  < 2e-16 ***
# climData.sub$maxswe 0.016407   0.002287   7.173 1.29e-12 ***


reg2 <- lm(soilTData$snowcovTs20mean~climData.sub$decSWEmean)
plot(climData.sub$decSWEmean, soilTData$snowcovTs20mean)
abline(reg2)
summary(reg2)
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.436521   0.043144  10.118   <2e-16 ***
# climData.sub$decSWEmean 0.072310   0.008058   8.973   <2e-16 ***

reg3 <- lm(soilTData$snowfreeTs20mean~climData.sub$meltdoy)
plot(climData.sub$meltdoy, soilTData$snowfreeTs20mean)
abline(reg3)
summary(reg3)
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          15.207909   0.622216  24.442   <2e-16 ***
# climData.sub$meltdoy -0.022625   0.002656  -8.519   <2e-16 ***


reg4<- lm(soilTData$snowfreeTs20mean~climData.sub$totaldaysSC)
plot(climData.sub$totaldaysSC, soilTData$snowfreeTs20mean)
abline(reg4)
summary(reg4)
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              12.506778   0.308885   40.49   <2e-16 ***
# climData.sub$totaldaysSC -0.013133   0.001547   -8.49   <2e-16 ***
