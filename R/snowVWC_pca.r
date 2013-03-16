# snowVWC_pca.r
# This is the script for performing a pca and further inferential tests (using
# pca results) on below-snow soilVWC across the SNOTEL network. This pca is 
# designed for testing the predictors of mean below-snow soil VWC at 20cm depth.

source('getdata.r')

# climData.sub should have the same sites/years in the same order as soilTData
# and soilVWCData
# sum(soilVWCData[,1]-climData.sub[,1] )

# Create a dataframe of variables for pca - include monthly means
varframe <- data.frame(climData.sub$totaldaysSC, climData.sub$elev, 
		climData.sub$meltdoy,climData.sub$onsetdoy,
		climData.sub$maxswe, climData.sub$scovmat,
		climData.sub$octTairMean,climData.sub$novTairMean,
	       	climData.sub$decTairMean,climData.sub$janTairMean,
		climData.sub$febTairMean,climData.sub$marTairMean,
		climData.sub$aprTairMean,climData.sub$mayTairMean,
		climData.sub$octSWEmean,climData.sub$novSWEmean,
	       	climData.sub$decSWEmean,climData.sub$janSWEmean,
		climData.sub$febSWEmean,climData.sub$marSWEmean,
		climData.sub$aprSWEmean,climData.sub$maySWEmean,
		soilVWCData$preonsetVWC20,soilTData$preonsetTs20,
		climData.sub$preonsetTair)


# Shorten the variable names in varframe
names(varframe) = sub("climData.sub.","",names(varframe))
names(varframe) = sub("soilTData.","",names(varframe))
names(varframe) = sub("soilVWCData.","",names(varframe))

# What about the difference (MAST-MAT)?
#varframe$Diff20cm <- (soilTData$mast20cm - climData.sub$maat)

# Write varframe to file
# write.csv(varframe, file='vwcMV.dat', col.names=TRUE)

# Run the pca on centered, scaled data - make it a na-free matrix first
varframe.cmplt <- complete.cases(varframe) # Boolean indicating NA-free rows
vwc20.pca <- prcomp(data.matrix(varframe[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(vwc20.pca)
# Make a screeplot of this
screeplot(vwc20.pca)

# Now make some useful variables
sd <- vwc20.pca$sdev
loadings <- vwc20.pca$rotation
rownames(loadings) <- colnames(varframe)
scores <- vwc20.pca$x
rownames(scores) <- 1:dim(scores)[1]
# Look at loadings
loadings[,1:4]

# Make a biplot of PCA 1 & 2 - note that the scaling factor of 10 may need to
# be changed depending on the data set
plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2",
     type="n", xlim=c(min(scores[,1:2]), max(scores[,1:2])),
     ylim=c(min(scores[,1:2]), max(scores[,1:2])))
arrows(0,0,loadings[,1]*29,loadings[,2]*10, length=0.1,angle=20, col="red")
# Add the names of the variables to the plot - 1.2 scaling insures that labels
# are plotted just beyond the arrows
text(loadings[,1]*29*1.2,loadings[,2]*10*1.2, rownames(loadings), col="red", 
     cex=0.7)
# Now add the scores with rownames
text(scores[,1],scores[,2], rownames(scores), col="blue", cex=0.7)

# Do the same plot (sorta) using built-in biplot function
biplot(scores[,1:2], loadings[,1:2], cex=0.7)

#setEPS()
#postscript("biplot13.eps", onefile="FALSE", horizontal="FALSE", 
#	   height = 5, width = 5)
# Built-in biplot function - axes 1 & 3
biplot(cbind(scores[,1], scores[,3]), cbind(loadings[,1], loadings[,3]),
       cex=0.7, xlab='PC1', ylab='PC3')
#dev.off()

# Built-in biplot function - axes 2 & 3
biplot(scores[,2:3], loadings[,2:3], xlab=rownames(scores), cex=0.7)

# Are axes 1-4 normal?
qqnorm(scores[,1])
qqline(scores[,1])
plot(density(scores[,1]))# mostly, a little negative skew
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # yes
qqnorm(scores[,4])
qqline(scores[,4]) # Some positive skew 

# Assign scores for pc1-4 to thier observation (rows)
soilVWCData$pc1score <- NA
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc4score <- NA
soilVWCData$pc1score[varframe.cmplt] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt] <- scores[,3]
soilVWCData$pc4score[varframe.cmplt] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 4
# principal components
lm1 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilVWCData)
AIC(lm1) # -264.7826

# remove PC4 
lm2 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score, data=soilVWCData)
AIC(lm2) # 228.7072 - so it seems like pc4 is important

# Add site
lm3 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score+pc4score
	 + as.factor(siteVWC), data=soilVWCData)
AIC(lm3) # -1034.22

# Add year
lm4 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score+pc4score
	  + as.factor(yearVWC) + as.factor(siteVWC), data=soilVWCData)
AIC(lm4) # -1222.722

summary(lm1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.465341   0.006680  69.666  < 2e-16 ***
# pc1score     0.004880   0.002111   2.311    0.021 *  
# pc2score    -0.043088   0.002817 -15.294  < 2e-16 ***
# pc3score     0.021986   0.004316   5.094 4.18e-07 ***
# pc4score     0.035502   0.005713   6.214 7.58e-10 ***
# Multiple R-squared: 0.2341,	Adjusted R-squared: 0.231 
# F-statistic: 76.03 on 4 and 995 DF,  p-value: < 2.2e-16 

summary(lm2)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.465357   0.006805  68.389  < 2e-16 ***
# pc1score     0.004884   0.002151   2.271   0.0234 *  
# pc2score    -0.043073   0.002870 -15.008  < 2e-16 ***
# pc3score     0.022000   0.004397   5.004 6.63e-07 ***
# Multiple R-squared: 0.2044,	Adjusted R-squared: 0.202 
# F-statistic: 85.28 on 3 and 996 DF,  p-value: < 2.2e-16

summary(lm3) # not showing site
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.311795   0.076481   4.077 5.03e-05 ***
# pc1score               -0.012738   0.003202  -3.978 7.57e-05 ***
# pc2score               -0.038488   0.004785  -8.043 3.20e-15 ***
# pc3score                0.013536   0.003221   4.203 2.94e-05 ***
# pc4score                0.034209   0.004094   8.356 2.91e-16 ***
# Multiple R-squared: 0.7645,	Adjusted R-squared: 0.7022 
# F-statistic: 12.27 on 209 and 790 DF,  p-value: < 2.2e-16 

summary(lm4) # not showing site or year
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.316082   0.123782   2.554 0.010854 *  
# pc1score               -0.010507   0.005412  -1.941 0.052584 .  
# pc2score               -0.008836   0.007320  -1.207 0.227760    
# pc3score               -0.000189   0.003510  -0.054 0.957066    
# pc4score                0.037071   0.006217   5.963 3.75e-09 ***
# Multiple R-squared:  0.81,	Adjusted R-squared: 0.7557 
# F-statistic: 14.92 on 222 and 777 DF,  p-value: < 2.2e-16


# One thing that is unclear here is that year seems to have a big effect.
# It might be best to try this with some selected years, then exclude year
# from the multiple regression model.
# 2007
varframe.07 <- varframe[climData.sub$yearClim==2007,] # Boolean 2007 rows
varframe.cmplt <- complete.cases(varframe.07)
vwc20.07.pca <- prcomp(data.matrix(varframe.07[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(vwc20.07.pca)
# Make a screeplot of this
screeplot(vwc20.07.pca)

# Now make some useful variables
sd <- vwc20.07.pca$sdev
loadings <- vwc20.07.pca$rotation
rownames(loadings) <- colnames(varframe.07)
scores <- vwc20.07.pca$x
rownames(scores) <- 1:dim(scores)[1]
# Look at loadings
loadings[,1:4]

# Built-in biplot function - axes 1 & 2
biplot(scores[,1:2], loadings[,1:2], xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 1 & 3
biplot(cbind(scores[,1], scores[,3]), cbind(loadings[,1], loadings[,3]),
       xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 2 & 3
biplot(scores[,2:3], loadings[,2:3], xlab=rownames(scores), cex=0.7)

# Are axes 1-4 normal?
qqnorm(scores[,1])
qqline(scores[,1])
plot(density(scores[,1]))# mostly, a little negative skew
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # yes except for at the extremes
qqnorm(scores[,4])
qqline(scores[,4]) # mostly, but with some negative skew 

# Assign scores for pc1-4 to thier observation (rows)
varframe.cmplt.07 <- (soilVWCData$yearVWC==2007 & complete.cases(varframe))
soilVWCData$pc1score <- NA
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc4score <- NA
soilVWCData$pc1score[varframe.cmplt.07] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt.07] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt.07] <- scores[,3]
soilVWCData$pc4score[varframe.cmplt.07] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 4
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilVWCData)
AIC(lm1) # -77.78469
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.488524   0.015308  31.913  < 2e-16 ***
# pc1score    0.004897   0.004425   1.107  0.27061    
# pc2score    0.020251   0.006325   3.202  0.00173 ** 
# pc3score    0.059976   0.009899   6.059 1.47e-08 ***
# pc4score    0.086471   0.014463   5.979 2.16e-08 ***
# Multiple R-squared: 0.4031,	Adjusted R-squared: 0.3842 
# F-statistic: 21.27 on 4 and 126 DF,  p-value: 2.01e-13 

lm2 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score, data=soilVWCData)
AIC(lm2) # -47.06656
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.491005   0.017270  28.432  < 2e-16 ***
# pc1score    0.005315   0.004993   1.064  0.28919    
# pc2score    0.021192   0.007136   2.970  0.00357 ** 
# pc3score    0.059972   0.011172   5.368 3.65e-07 ***
# Multiple R-squared: 0.2338,	Adjusted R-squared: 0.2157 
# F-statistic: 12.91 on 3 and 127 DF,  p-value: 2.036e-07 

# 2009
varframe.09 <- varframe[climData.sub$yearClim==2009,] # Boolean 2009 rows
varframe.cmplt <- complete.cases(varframe.09)
vwc20.09.pca <- prcomp(data.matrix(varframe.09[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(vwc20.09.pca)
# Make a screeplot of this
screeplot(vwc20.09.pca)

# Now make some useful variables
sd <- vwc20.09.pca$sdev
loadings <- vwc20.09.pca$rotation
rownames(loadings) <- colnames(varframe.09)
scores <- vwc20.09.pca$x
rownames(scores) <- 1:dim(scores)[1]
# Look at loadings
loadings[,1:4]

# Built-in biplot function - axes 1 & 2
biplot(scores[,1:2], loadings[,1:2], xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 1 & 3
biplot(cbind(scores[,1], scores[,3]), cbind(loadings[,1], loadings[,3]),
       xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 2 & 3
biplot(scores[,2:3], loadings[,2:3], xlab=rownames(scores), cex=0.7)

# Are axes 1-4 normal?
qqnorm(scores[,1])
qqline(scores[,1])
plot(density(scores[,1]))# mostly, a little negative skew
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # bit of positive skew
qqnorm(scores[,4])
qqline(scores[,4]) # yes, except at extremes 

# Assign scores for pc1-4 to thier observation (rows)
varframe.cmplt.09 <- (soilVWCData$yearVWC==2009 & complete.cases(varframe))
soilVWCData$pc1score <- NA
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc4score <- NA
soilVWCData$pc1score[varframe.cmplt.09] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt.09] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt.09] <- scores[,3]
soilVWCData$pc4score[varframe.cmplt.09] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 4
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilVWCData)
AIC(lm1) # -82.17018
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.444328   0.014780  30.063  < 2e-16 ***
# pc1score    0.010609   0.004303   2.466   0.0149 *  
# pc2score    0.045937   0.005912   7.770 1.52e-12 ***
# pc3score    0.021464   0.009738   2.204   0.0291 *  
# pc4score    0.082046   0.015082   5.440 2.31e-07 ***
# Multiple R-squared: 0.4189,	Adjusted R-squared: 0.4023 
# F-statistic: 25.23 on 4 and 140 DF,  p-value: 9.561e-16 

lm2 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score, data=soilVWCData)
AIC(lm2) # -56.36391
summary(lm2)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.444328   0.016209  27.412  < 2e-16 ***
# pc1score    0.010609   0.004719   2.248   0.0261 *  
# pc2score    0.045937   0.006483   7.085 6.08e-11 ***
# pc3score    0.021464   0.010680   2.010   0.0464 *  
# Multiple R-squared: 0.296,	Adjusted R-squared: 0.2811 
# F-statistic: 19.76 on 3 and 141 DF,  p-value: 9.417e-11 


# 2011
varframe.11 <- varframe[climData.sub$yearClim==2011,] # Boolean 2011 rows
varframe.cmplt <- complete.cases(varframe.11)
vwc20.11.pca <- prcomp(data.matrix(varframe.11[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(vwc20.11.pca)
# Make a screeplot of this
screeplot(vwc20.11.pca)

# Now make some useful variables
sd <- vwc20.11.pca$sdev
loadings <- vwc20.11.pca$rotation
rownames(loadings) <- colnames(varframe.11)
scores <- vwc20.11.pca$x
rownames(scores) <- 1:dim(scores)[1]
# Look at loadings
loadings[,1:4]

# Built-in biplot function - axes 1 & 2
biplot(scores[,1:2], loadings[,1:2], xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 1 & 3
biplot(cbind(scores[,1], scores[,3]), cbind(loadings[,1], loadings[,3]),
       xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 2 & 3
biplot(scores[,2:3], loadings[,2:3], xlab=rownames(scores), cex=0.7)

# Are axes 1-4 normal?
qqnorm(scores[,1])
qqline(scores[,1])
plot(density(scores[,1]))# yes
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # yes except for at the extremes, some neg. skew
qqnorm(scores[,4])
qqline(scores[,4]) # yes

# Assign scores for pc1-4 to thier observation (rows)
varframe.cmplt.11 <- (soilVWCData$yearVWC==2011 & complete.cases(varframe))
soilVWCData$pc1score <- NA
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc4score <- NA
soilVWCData$pc1score[varframe.cmplt.11] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt.11] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt.11] <- scores[,3]
soilVWCData$pc4score[varframe.cmplt.11] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilVWCData)
AIC(lm1) # -68.46771
summary(lm1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.569571   0.014867  38.310  < 2e-16 ***
# pc1score     0.019922   0.004251   4.687 5.79e-06 ***
# pc2score     0.042148   0.006005   7.018 5.54e-11 ***
# pc3score     0.016914   0.009290   1.821   0.0705 .  
# pc4score    -0.018553   0.016086  -1.153   0.2504    
# Multiple R-squared: 0.315,	Adjusted R-squared: 0.2984 
# F-statistic: 18.97 on 4 and 165 DF,  p-value: 7.539e-13 

lm2 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score, data=soilVWCData)
AIC(lm2) # -69.10273
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.569571   0.014882  38.272  < 2e-16 ***
# pc1score    0.019922   0.004255   4.682 5.88e-06 ***
# pc2score    0.042148   0.006011   7.011 5.67e-11 ***
# pc3score    0.016914   0.009299   1.819   0.0707 .  
# Multiple R-squared: 0.3094,	Adjusted R-squared: 0.297 
# F-statistic:  24.8 on 3 and 166 DF,  p-value: 2.62e-13

