# gsVWC_pca
# This is the script for performing a pca and further inferential tests (using
# pca results) on growing season soil VWC across the SNOTEL network. This pca 
# is designed for testing the predictors of mean growing season soil VWC at 
# 20cm depth.

source('getdata.r')

# climData.sub should have the same sites/years in the same order as soilVWCData
# and soilVWCData
# sum(soilVWCData[,1]-climData.sub[,1] )

# Create a dataframe of variables for pca - include monthly means
varframe <- data.frame(climData.sub$totaldaysSC,climData.sub$jasMAT,
		       climData.sub$meltdoy,climData.sub$maxswe, 
		       climData.sub$scovmat,soilTData$jfmTs5mean,
		       climData.sub$aprTairMean,climData.sub$mayTairMean,
		       climData.sub$junTairMean,climData.sub$julTairMean,
		       climData.sub$augTairMean,climData.sub$sepTairMean,
		       climData.sub$JASprecip,climData.sub$mayPrecip,
		       climData.sub$junPrecip,climData.sub$julPrecip,
		       climData.sub$augPrecip,climData.sub$sepPrecip,
		       soilVWCData$jfmVWC20mean,climData.sub$elev)


# Shorten the variable names in varframe
names(varframe) = sub("climData.sub.","",names(varframe))
names(varframe) = sub("soilTData.","",names(varframe))
names(varframe) = sub("soilVWCData.","",names(varframe))

# What about the difference (MAST-MAT)?
#varframe$Diff20cm <- (soilVWCData$mast20cm - climData.sub$maat)

# Write varframe to file
# write.csv(varframe, file='TsoilMV.dat', col.names=TRUE)

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
qqline(scores[,2]) # mostly, some positive skew
qqnorm(scores[,3])
qqline(scores[,3]) # yes, slight negative skew

# Assign scores for pc1, 2, & 3 to thier observation (rows)
soilVWCData$pc1score <- NA 
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc1score[varframe.cmplt] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components
lm1 <- lm(jasVWC20mean ~ pc1score + pc2score + pc3score, data=soilVWCData)

# Add site
lm2 <- lm(jasVWC20mean ~ pc1score + pc2score + pc3score
	 + as.factor(siteVWC), data=soilVWCData)

# Add year
lm3 <- lm(jasVWC20mean ~ pc1score + pc2score + pc3score + as.factor(yearVWC)
	 + as.factor(siteVWC), data=soilVWCData)

AIC(lm1) # -876.1643
AIC(lm2) # -2093.832
AIC(lm3) # -2173.575

summary(lm1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.275297   0.004917  55.989  < 2e-16 ***
# pc1score    -0.021617   0.001786 -12.105  < 2e-16 ***
# pc2score     0.007996   0.002564   3.118  0.00187 ** 
# pc3score    -0.019927   0.003331  -5.983 3.04e-09 ***
# Multiple R-squared: 0.1589,	Adjusted R-squared: 0.1564 
# F-statistic:  63.9 on 3 and 1015 DF,  p-value: < 2.2e-16 

summary(lm2)
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.1415014  0.0462187   3.062 0.002274 ** 
# pc1score               -0.0373956  0.0035810 -10.443  < 2e-16 ***
# pc2score                0.0010156  0.0024574   0.413 0.679518    
# pc3score               -0.0386070  0.0030351 -12.720  < 2e-16 ***
# Multiple R-squared: 0.8287,	Adjusted R-squared: 0.7855 
# F-statistic: 19.19 on 205 and 813 DF,  p-value: < 2.2e-16 

summary(lm3)
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.0489963  0.1032221   0.475 0.635153    
# pc1score               -0.0439191  0.0057158  -7.684 4.51e-14 ***
# pc2score                0.0136366  0.0036435   3.743 0.000195 ***
# pc3score               -0.0389700  0.0043436  -8.972  < 2e-16 ***
# Multiple R-squared: 0.8453,	Adjusted R-squared: 0.8034 
# F-statistic: 20.17 on 217 and 801 DF,  p-value: < 2.2e-16 

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
qqline(scores[,2]) # mostly, tiny skew
qqnorm(scores[,3])
qqline(scores[,3]) # yes except for at the extremes

# Assign scores for pc1, 2, & 3 to thier observation (rows)
varframe.cmplt.07 <- (soilVWCData$yearVWC==2007 & complete.cases(varframe))
soilVWCData$pc1score <- NA 
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc1score[varframe.cmplt.07] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt.07] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt.07] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(jasVWC20mean ~ pc1score + pc2score + pc3score, data=soilVWCData)
summary(lm1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.229824   0.011550  19.898  < 2e-16 ***
# pc1score    -0.023541   0.003709  -6.348 3.24e-09 ***
# pc2score     0.021379   0.006751   3.167  0.00191 ** 
# pc3score    -0.002858   0.007715  -0.370  0.71165    
# Multiple R-squared: 0.2765,	Adjusted R-squared: 0.2601 
# F-statistic: 16.82 on 3 and 132 DF,  p-value: 2.599e-09

lm2 <- lm(jasVWC20mean ~ pc1score + pc2score, data=soilVWCData)
summary(lm2)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.229824   0.011513  19.963  < 2e-16 ***
# pc1score    -0.023541   0.003697  -6.368 2.87e-09 ***
# pc2score     0.021379   0.006729   3.177  0.00185 ** 
# Multiple R-squared: 0.2758,	Adjusted R-squared: 0.2649 
# F-statistic: 25.32 on 2 and 133 DF,  p-value: 4.795e-10

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
plot(density(scores[,1]))# yes, except for a little at the extremes
qqnorm(scores[,2])
qqline(scores[,2]) # yes, tiny positive skew
qqnorm(scores[,3])
qqline(scores[,3]) # yes

# Assign scores for pc1, 2, & 3 to thier observation (rows)
varframe.cmplt.09 <- (soilVWCData$yearVWC==2009 & complete.cases(varframe))
soilVWCData$pc1score <- NA 
soilVWCData$pc2score <- NA
soilVWCData$pc3score <- NA
soilVWCData$pc1score[varframe.cmplt.09] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt.09] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt.09] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)

lm1 <- lm(jasVWC20mean ~ pc1score + pc2score + pc3score, data=soilVWCData)
summary(lm1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.271016   0.011623  23.317  < 2e-16 ***
# pc1score    -0.013402   0.003795  -3.531 0.000556 ***
# pc2score    -0.001580   0.006509  -0.243 0.808547    
# pc3score     0.009857   0.008500   1.160 0.248128    
# Multiple R-squared: 0.08787,	Adjusted R-squared: 0.06887 
# F-statistic: 4.624 on 3 and 144 DF,  p-value: 0.004049

lm2 <- lm(jasVWC20mean ~ pc1score, data=soilVWCData)
summary(lm2)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.271016   0.011599  23.365  < 2e-16 ***
# pc1score    -0.013402   0.003787  -3.538  0.00054 ***
# Multiple R-squared: 0.07898,	Adjusted R-squared: 0.07268 
# F-statistic: 12.52 on 1 and 146 DF,  p-value: 0.0005401 

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
soilVWCData$pc1score[varframe.cmplt.11] <- scores[,1]
soilVWCData$pc2score[varframe.cmplt.11] <- scores[,2]
soilVWCData$pc3score[varframe.cmplt.11] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(jasVWC20mean ~ pc1score + pc2score + pc3score, data=soilVWCData)
summary(lm1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.339582   0.012924  26.276  < 2e-16 ***
# pc1score    -0.025304   0.004442  -5.696 5.24e-08 ***
# pc2score    -0.008621   0.006249  -1.380  0.16949    
# pc3score    -0.023206   0.008335  -2.784  0.00597 ** 
# Multiple R-squared: 0.1976,	Adjusted R-squared: 0.1835 
# F-statistic: 14.04 on 3 and 171 DF,  p-value: 3.195e-08

lm2 <- lm(jasVWC20mean ~ pc1score + pc3score, data=soilVWCData)
summary(lm2)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.339582   0.012957  26.207  < 2e-16 ***
# pc1score    -0.025304   0.004454  -5.682 5.59e-08 ***
# pc3score    -0.023206   0.008357  -2.777   0.0061 ** 
# Multiple R-squared: 0.1886,	Adjusted R-squared: 0.1792 
# F-statistic:    20 on 2 and 172 DF,  p-value: 1.556e-08 
