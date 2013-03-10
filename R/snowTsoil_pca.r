# snowTsoil_pca
# This is the script for performing a pca and further inferential tests (using
# pca results) on below-snow Tsoil across the SNOTEL network. This pca is 
# designed for testing the predictors of mean below-snow Tsoil at 20cm depth.

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
# write.csv(varframe, file='TsoilMV.dat', col.names=TRUE)

# Run the pca on centered, scaled data - make it a na-free matrix first
varframe.cmplt <- complete.cases(varframe) # Boolean indicating NA-free rows
ts20.pca <- prcomp(data.matrix(varframe[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(ts20.pca)
# Make a screeplot of this
screeplot(ts20.pca)

# Now make some useful variables
sd <- ts20.pca$sdev
loadings <- ts20.pca$rotation
rownames(loadings) <- colnames(varframe)
scores <- ts20.pca$x
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
qqline(scores[,3]) # yes
qqnorm(scores[,4])
qqline(scores[,4]) # Some positive skew 

# Assign scores for pc1-4 to thier observation (rows)
soilTData$pc1score <- NA
soilTData$pc2score <- NA
soilTData$pc3score <- NA
soilTData$pc4score <- NA
soilTData$pc1score[varframe.cmplt] <- scores[,1]
soilTData$pc2score[varframe.cmplt] <- scores[,2]
soilTData$pc3score[varframe.cmplt] <- scores[,3]
soilTData$pc4score[varframe.cmplt] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 4
# principal components
lm1 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilTData)
AIC(lm1) # 2170.249

# remove PC4 
lm2 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score, data=soilTData)
AIC(lm2) # 2208.249 - so it seems like pc4 is important

# Add site
lm3 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score+pc4score
	 + as.factor(siteTsoil), data=soilTData)
AIC(lm3) # 1594.67

# Add year
lm4 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score+pc4score
	  + as.factor(yearTsoil) + as.factor(siteTsoil), data=soilTData)
AIC(lm4) # 1349.05


summary(lm1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.734750   0.022570  32.555  < 2e-16 ***
# pc1score    -0.020627   0.007134  -2.892  0.00392 ** 
# pc2score    -0.140904   0.009520 -14.801  < 2e-16 ***
# pc3score    -0.032726   0.014582  -2.244  0.02504 *  
# pc4score     0.122926   0.019304   6.368 2.92e-10 ***
# Multiple R-squared: 0.215,	Adjusted R-squared: 0.2119 
# F-statistic: 68.14 on 4 and 995 DF,  p-value: < 2.2e-16

summary(lm2)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.734803   0.023013  31.930  < 2e-16 ***
# pc1score    -0.020614   0.007274  -2.834  0.00469 ** 
# pc2score    -0.140855   0.009707 -14.511  < 2e-16 ***
# pc3score    -0.032677   0.014869  -2.198  0.02820 *  
# Multiple R-squared: 0.183,	Adjusted R-squared: 0.1806 
# F-statistic: 74.39 on 3 and 996 DF,  p-value: < 2.2e-16

summary(lm3) # not showing site
                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               0.536762   0.284715   1.885 0.059762 .  
# pc1score                 -0.049195   0.011919  -4.127 4.06e-05 ***
# pc2score                 -0.080106   0.017813  -4.497 7.93e-06 ***
# pc3score                  0.003148   0.011990   0.263 0.792977    
# pc4score                  0.119407   0.015240   7.835 1.51e-14 ***
# Multiple R-squared: 0.707,	Adjusted R-squared: 0.6295 
# F-statistic: 9.123 on 209 and 790 DF,  p-value: < 2.2e-16 

summary(lm4) # not showing site or year
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               1.353485   0.447829   3.022 0.002591 ** 
# pc1score                 -0.118783   0.019581  -6.066 2.04e-09 ***
# pc2score                  0.045074   0.026482   1.702 0.089145 .  
# pc3score                 -0.066593   0.012699  -5.244 2.03e-07 ***
# pc4score                  0.050914   0.022491   2.264 0.023868 *  
# Multiple R-squared: 0.7767,	Adjusted R-squared: 0.7129 
# F-statistic: 12.18 on 222 and 777 DF,  p-value: < 2.2e-16 



# One thing that is unclear here is that year seems to have a big effect.
# It might be best to try this with some selected years, then exclude year
# from the multiple regression model.
# 2007
varframe.07 <- varframe[climData.sub$yearClim==2007,] # Boolean 2007 rows
varframe.cmplt <- complete.cases(varframe.07)
ts20.07.pca <- prcomp(data.matrix(varframe.07[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(ts20.07.pca)
# Make a screeplot of this
screeplot(ts20.07.pca)

# Now make some useful variables
sd <- ts20.07.pca$sdev
loadings <- ts20.07.pca$rotation
rownames(loadings) <- colnames(varframe.07)
scores <- ts20.07.pca$x
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
varframe.cmplt.07 <- (soilTData$yearTsoil==2007 & complete.cases(varframe))
soilTData$pc1score <- NA
soilTData$pc2score <- NA
soilTData$pc3score <- NA
soilTData$pc4score <- NA
soilTData$pc1score[varframe.cmplt.07] <- scores[,1]
soilTData$pc2score[varframe.cmplt.07] <- scores[,2]
soilTData$pc3score[varframe.cmplt.07] <- scores[,3]
soilTData$pc4score[varframe.cmplt.07] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 4
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilTData)
AIC(lm1) # 235.0844
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.65950    0.05053  13.052  < 2e-16 ***
# pc1score    -0.01521    0.01461  -1.041    0.300    
# pc2score     0.08776    0.02088   4.203 4.95e-05 ***
# pc3score    -0.04411    0.03268  -1.350    0.179    
# pc4score     0.25760    0.04774   5.396 3.25e-07 ***
# Multiple R-squared: 0.287,	Adjusted R-squared: 0.2643 
# F-statistic: 12.68 on 4 and 126 DF,  p-value: 1.061e-08 


lm2 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score, data=soilTData)
AIC(lm2) # 260.3205 
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.66690    0.05582  11.947  < 2e-16 ***
# pc1score    -0.01396    0.01614  -0.865 0.388714    
# pc2score     0.09056    0.02307   3.926 0.000141 ***
# pc3score    -0.04413    0.03611  -1.222 0.223984    
# Multiple R-squared: 0.1222,	Adjusted R-squared: 0.1015 
# F-statistic: 5.894 on 3 and 127 DF,  p-value: 0.0008469 

# 2009
varframe.09 <- varframe[climData.sub$yearClim==2009,] # Boolean 2009 rows
varframe.cmplt <- complete.cases(varframe.09)
ts20.09.pca <- prcomp(data.matrix(varframe.09[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(ts20.09.pca)
# Make a screeplot of this
screeplot(ts20.09.pca)

# Now make some useful variables
sd <- ts20.09.pca$sdev
loadings <- ts20.09.pca$rotation
rownames(loadings) <- colnames(varframe.09)
scores <- ts20.09.pca$x
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
varframe.cmplt.09 <- (soilTData$yearTsoil==2009 & complete.cases(varframe))
soilTData$pc1score <- NA
soilTData$pc2score <- NA
soilTData$pc3score <- NA
soilTData$pc4score <- NA
soilTData$pc1score[varframe.cmplt.09] <- scores[,1]
soilTData$pc2score[varframe.cmplt.09] <- scores[,2]
soilTData$pc3score[varframe.cmplt.09] <- scores[,3]
soilTData$pc4score[varframe.cmplt.09] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 4
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilTData)
AIC(lm1) # 307.7004
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.54842    0.05669   9.674  < 2e-16 ***
# pc1score    -0.01330    0.01650  -0.806 0.421834    
# pc2score     0.15700    0.02268   6.924 1.46e-10 ***
# pc3score    -0.14195    0.03735  -3.800 0.000215 ***
# pc4score     0.17070    0.05785   2.951 0.003719 ** 
# Multiple R-squared: 0.3388,	Adjusted R-squared: 0.3199 
# F-statistic: 17.93 on 4 and 140 DF,  p-value: 6.556e-12 

lm2 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score, data=soilTData)
AIC(lm2) # 314.4483
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.54842    0.05822   9.420  < 2e-16 ***
# pc1score    -0.01330    0.01695  -0.784 0.434087    
# pc2score     0.15700    0.02329   6.742 3.71e-10 ***
# pc3score    -0.14195    0.03836  -3.700 0.000308 ***
# Multiple R-squared: 0.2977,	Adjusted R-squared: 0.2827 
# F-statistic: 19.92 on 3 and 141 DF,  p-value: 8.012e-11 

# 2011
varframe.11 <- varframe[climData.sub$yearClim==2011,] # Boolean 2011 rows
varframe.cmplt <- complete.cases(varframe.11)
ts20.11.pca <- prcomp(data.matrix(varframe.11[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(ts20.11.pca)
# Make a screeplot of this
screeplot(ts20.11.pca)

# Now make some useful variables
sd <- ts20.11.pca$sdev
loadings <- ts20.11.pca$rotation
rownames(loadings) <- colnames(varframe.11)
scores <- ts20.11.pca$x
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
varframe.cmplt.11 <- (soilTData$yearTsoil==2011 & complete.cases(varframe))
soilTData$pc1score <- NA
soilTData$pc2score <- NA
soilTData$pc3score <- NA
soilTData$pc4score <- NA
soilTData$pc1score[varframe.cmplt.11] <- scores[,1]
soilTData$pc2score[varframe.cmplt.11] <- scores[,2]
soilTData$pc3score[varframe.cmplt.11] <- scores[,3]
soilTData$pc4score[varframe.cmplt.11] <- scores[,4]

# Now it should be possible to do a multiple regression model with the 4
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilTData)
AIC(lm1) # 240.6223
summary(lm1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.877093   0.036902  23.768  < 2e-16 ***
# pc1score     0.034941   0.010551   3.312  0.00114 ** 
# pc2score     0.109620   0.014906   7.354 8.53e-12 ***
# pc3score    -0.049654   0.023059  -2.153  0.03274 *  
# pc4score    -0.006519   0.039927  -0.163  0.87050    
# Multiple R-squared: 0.297,	Adjusted R-squared:  0.28 
# F-statistic: 17.43 on 4 and 165 DF,  p-value: 6.021e-12

lm2 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score, data=soilTData)
AIC(lm2) # 238.6498
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.87709    0.03679  23.838  < 2e-16 ***
# pc1score     0.03494    0.01052   3.321   0.0011 ** 
# pc2score     0.10962    0.01486   7.376  7.4e-12 ***
# pc3score    -0.04965    0.02299  -2.160   0.0322 *  
# Multiple R-squared: 0.2969,	Adjusted R-squared: 0.2842 
# F-statistic: 23.37 on 3 and 166 DF,  p-value: 1.145e-12 
