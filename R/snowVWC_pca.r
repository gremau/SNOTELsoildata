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
loadings[,1:3]

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


# Are axes 1, 2, and 3 normal?
qqnorm(scores[,1])
qqline(scores[,1])
plot(density(scores[,1]))# mostly, a little negative skew
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # yes

# Assign scores for pc1, 2, & 3 to thier observation (rows)
soilVWCData$pc1score <- NA
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc4score <- NA
soilVWCData$pc1score[varframe.cmplt] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt] <- scores[,3]
soilVWCData$pc4score[varframe.cmplt] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 3
# principal components
lm1 <- lm(jfmVWC20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilVWCData)

lm2 <- lm(jfmVWC20mean ~ pc1score + pc2score + pc3score
	 + as.factor(siteVWC), data=soilVWCData)

lm3 <- lm(jfmVWC20mean ~ pc1score + pc2score, data=soilVWCData)

lm4 <- lm(jfmVWC20mean ~ pc1score + pc2score
	 + as.factor(siteVWC), data=soilVWCData)

lm5 <- lm(jfmVWC20mean ~ pc1score + pc2score + pc3score
	  + as.factor(yearVWC) + as.factor(siteVWC), data=soilVWCData)

lm6 <- lm(jfmVWC20mean ~ pc1score + as.factor(yearVWC)
	 + as.factor(siteVWC) , data=soilVWCData)


AIC(lm1) # -228.7072
AIC(lm2) # -951.5271
AIC(lm3) # -205.8769
AIC(lm4) # -936.2304
AIC(lm5) # -1179.974
AIC(lm6) # -1179.317

summary(lm1) 
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.465357   0.006805  68.389  < 2e-16 ***
# pc1score     0.004884   0.002151   2.271   0.0234 *  
# pc2score    -0.043073   0.002870 -15.008  < 2e-16 ***
# pc3score     0.022000   0.004397   5.004 6.63e-07 ***
# Multiple R-squared: 0.2044,	Adjusted R-squared: 0.202 
# F-statistic: 85.28 on 3 and 996 DF,  p-value: < 2.2e-16 

summary(lm2) # site not shown
                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.3284683  0.0797118   4.121 4.18e-05 ***
# pc1score               -0.0085704  0.0032974  -2.599 0.009519 ** 
# pc2score               -0.0485671  0.0048278 -10.060  < 2e-16 ***
# pc3score                0.0124643  0.0033552   3.715 0.000218 ***
# Multiple R-squared: 0.7437,	Adjusted R-squared: 0.6763 
# F-statistic: 11.04 on 208 and 791 DF,  p-value: < 2.2e-16 

summary(lm5)) # site not shown
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.298800   0.126467   2.363 0.018389 *  
# pc1score               -0.015388   0.005467  -2.815 0.005008 ** 
# pc2score               -0.010161   0.007477  -1.359 0.174558    
# pc3score               -0.004259   0.003519  -1.210 0.226524
# Multiple R-squared: 0.8013,	Adjusted R-squared: 0.7448 
# F-statistic:  14.2 on 221 and 778 DF,  p-value: < 2.2e-16 

summary(lm6)) # site not shown
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.249046   0.118970   2.093 0.036640 *  
# pc1score               -0.020286   0.003759  -5.397 8.98e-08 ***
# Multiple R-squared: 0.8004,	Adjusted R-squared: 0.7443 
# F-statistic: 14.28 on 219 and 780 DF,  p-value: < 2.2e-16 

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

# Built-in biplot function - axes 1 & 2
biplot(scores[,1:2], loadings[,1:2], xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 1 & 3
biplot(cbind(scores[,1], scores[,3]), cbind(loadings[,1], loadings[,3]),
       xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 2 & 3
biplot(scores[,2:3], loadings[,2:3], xlab=rownames(scores), cex=0.7)

# Are axes 1, 2, and 3 normal?
qqnorm(scores[,1])
qqline(scores[,1])
plot(density(scores[,1]))# yes
qqnorm(scores[,2])
qqline(scores[,2]) # yes, with a little positive skew
qqnorm(scores[,3])
qqline(scores[,3]) # sorta, not at extremes though

# Assign scores for pc1, 2, & 3 to thier observation (rows)
varframe.cmplt.07 <- (soilVWCData$yearVWC==2007 & complete.cases(varframe))
soilVWCData$pc1score <- NA
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc4score <- NA
soilVWCData$pc1score[varframe.cmplt.07] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt.07] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt.07] <- scores[,3]
soilVWCData$pc4score[varframe.cmplt.07] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(jfmVWC20mean ~ pc1score + pc2score + pc3score, data=soilVWCData)
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.491005   0.017270  28.432  < 2e-16 ***
# pc1score    0.005315   0.004993   1.064  0.28919    
# pc2score    0.021192   0.007136   2.970  0.00357 ** 
# pc3score    0.059972   0.011172   5.368 3.65e-07 ***
# Multiple R-squared: 0.2338,	Adjusted R-squared: 0.2157 
# F-statistic: 12.91 on 3 and 127 DF,  p-value: 2.036e-07 

lm2 <- lm(jfmVWC20mean ~ pc2score + pc3score, data=soilVWCData)
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.49092    0.01728  28.413  < 2e-16 ***
# pc2score     0.02116    0.00714   2.964  0.00363 ** 
# pc3score     0.05997    0.01118   5.365 3.66e-07 ***
# Multiple R-squared: 0.2269,	Adjusted R-squared: 0.2148 
# F-statistic: 18.79 on 2 and 128 DF,  p-value: 7.021e-08 

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

loadings[,1:3]

# Built-in biplot function - axes 1 & 2
biplot(scores[,1:2], loadings[,1:2], xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 1 & 3
biplot(cbind(scores[,1], scores[,3]), cbind(loadings[,1], loadings[,3]),
       xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 2 & 3
biplot(scores[,2:3], loadings[,2:3], xlab=rownames(scores), cex=0.7)

# Are axes 1, 2, and 3 normal?
qqnorm(scores[,1])
qqline(scores[,1])
plot(density(scores[,1]))# mostly
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # sorta, definately positive skew

# Assign scores for pc1, 2, & 3 to thier observation (rows)
varframe.cmplt.09 <- (soilVWCData$yearVWC==2009 & complete.cases(varframe))
soilVWCData$pc1score <- NA
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc4score <- NA
soilVWCData$pc1score[varframe.cmplt.09] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt.09] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt.09] <- scores[,3]
soilVWCData$pc4score[varframe.cmplt.09] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(jfmVWC20mean ~ pc1score + pc2score + pc3score, data=soilVWCData)
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
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

loadings[,1:3]

# Built-in biplot function - axes 1 & 2
biplot(scores[,1:2], loadings[,1:2], xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 1 & 3
biplot(cbind(scores[,1], scores[,3]), cbind(loadings[,1], loadings[,3]),
       xlab=rownames(scores), cex=0.7)
# Built-in biplot function - axes 2 & 3
biplot(scores[,2:3], loadings[,2:3], xlab=rownames(scores), cex=0.7)

# Are axes 1, 2, and 3 normal?
qqnorm(scores[,1])
qqline(scores[,1])
plot(density(scores[,1]))# yes
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # yes except for at the extremes

# Assign scores for pc1, 2, & 3 to thier observation (rows)
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
lm1 <- lm(jfmVWC20mean ~ pc1score + pc2score + pc3score, data=soilVWCData)
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.569571   0.014882  38.272  < 2e-16 ***
# pc1score    0.019922   0.004255   4.682 5.88e-06 ***
# pc2score    0.042148   0.006011   7.011 5.67e-11 ***
# pc3score    0.016914   0.009299   1.819   0.0707 .  
# Multiple R-squared: 0.3094,	Adjusted R-squared: 0.297 
# F-statistic:  24.8 on 3 and 166 DF,  p-value: 2.62e-13 

lm2 <- lm(jfmVWC20mean ~ pc1score + pc2score, data=soilVWCData)
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.569571   0.014985  38.010  < 2e-16 ***
# pc1score    0.019922   0.004285   4.650 6.72e-06 ***
# pc2score    0.042148   0.006053   6.963 7.27e-11 ***
# Multiple R-squared: 0.2957,	Adjusted R-squared: 0.2872 
# F-statistic: 35.05 on 2 and 167 DF,  p-value: 1.944e-13 
