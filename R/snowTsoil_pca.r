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
plot(density(scores[,1]))# mostly, a little negative skew
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # yes

# Assign scores for pc1, 2, & 3 to thier observation (rows)
soilTData$pc1score <- NA
soilTData$pc2score <- NA
soilTData$pc3score <- NA
soilTData$pc1score[varframe.cmplt] <- scores[,1]
soilTData$pc2score[varframe.cmplt] <- scores[,2]
soilTData$pc3score[varframe.cmplt] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components
lm1 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score, data=soilTData)

# Add site
lm2 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score
	 + as.factor(siteTsoil), data=soilTData)

# Add year
lm3 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score
	  + as.factor(yearTsoil) + as.factor(siteTsoil), data=soilTData)

# Remove pc2 (keep site)
lm4 <- lm(snowcovTs20mean ~ pc1score + pc3score + as.factor(yearTsoil)
	 + as.factor(siteTsoil), data=soilTData)

# Add year
lm5 <- lm(snowcovTs20mean ~ pc1score + pc2score + as.factor(yearTsoil)
	 + as.factor(siteTsoil) , data=soilTData)
# Add pc3 back in
lm6 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score
	  + as.factor(yearTsoil) + as.factor(siteTsoil), data=soilTData)

AIC(lm1) # 2208.249
AIC(lm2) # 1667.501
AIC(lm3) # 1353.623
AIC(lm4) # 1355.031

summary(lm1)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.734803   0.023013  31.930  < 2e-16 ***
# pc1score    -0.020614   0.007274  -2.834  0.00469 ** 
# pc2score    -0.140855   0.009707 -14.511  < 2e-16 ***
# pc3score    -0.032677   0.014869  -2.198  0.02820 *  
# Multiple R-squared: 0.183,	Adjusted R-squared: 0.1806 
# F-statistic: 74.39 on 3 and 996 DF,  p-value: < 2.2e-16

summary(lm2) # not showing site
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               0.5949591  0.2952824   2.015  0.04425 *  
# pc1score                 -0.0346480  0.0122147  -2.837  0.00468 ** 
# pc2score                 -0.1152891  0.0178841  -6.446 1.99e-10 ***
# pc3score                 -0.0005928  0.0124289  -0.048  0.96197 
# Multiple R-squared: 0.6843,	Adjusted R-squared: 0.6013 
# F-statistic: 8.242 on 208 and 791 DF,  p-value: < 2.2e-16

summary(lm3) # not showing site or year
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               1.329751   0.448891   2.962 0.003147 ** 
# pc1score                 -0.125487   0.019407  -6.466 1.77e-10 ***
# pc2score                  0.043254   0.026540   1.630 0.103556    
# pc3score                 -0.072183   0.012490  -5.779 1.08e-08 ***
# Multiple R-squared: 0.7753,	Adjusted R-squared: 0.7114 
# F-statistic: 12.14 on 221 and 778 DF,  p-value: < 2.2e-16 

summary(lm4) # not showing site or year
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               1.5772353  0.4228752   3.730 0.000206 ***
# pc1score                 -0.1025905  0.0134023  -7.655 5.74e-14 ***
# pc3score                 -0.0703813  0.0124538  -5.651 2.23e-08 ***
# Multiple R-squared: 0.7745,	Adjusted R-squared: 0.7108 
# F-statistic: 12.16 on 220 and 779 DF,  p-value: < 2.2e-16 

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
plot(density(scores[,1]))# mostly, a little negative skew
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # yes except for at the extremes

# Assign scores for pc1, 2, & 3 to thier observation (rows)
varframe.cmplt.07 <- (soilTData$yearTsoil==2007 & complete.cases(varframe))
soilTData$pc1score <- NA
soilTData$pc2score <- NA
soilTData$pc3score <- NA
soilTData$pc1score[varframe.cmplt.07] <- scores[,1]
soilTData$pc2score[varframe.cmplt.07] <- scores[,2]
soilTData$pc3score[varframe.cmplt.07] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score, data=soilTData)
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.66690    0.05582  11.947  < 2e-16 ***
# pc1score    -0.01396    0.01614  -0.865 0.388714    
# pc2score     0.09056    0.02307   3.926 0.000141 ***
# pc3score    -0.04413    0.03611  -1.222 0.223984    
# Multiple R-squared: 0.1222,	Adjusted R-squared: 0.1015 
# F-statistic: 5.894 on 3 and 127 DF,  p-value: 0.0008469 

lm2 <- lm(snowcovTs20mean ~ pc2score, data=soilTData)
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.66711    0.05587  11.939  < 2e-16 ***
# pc2score     0.09064    0.02309   3.926  0.00014 ***
# Multiple R-squared: 0.1067,	Adjusted R-squared: 0.09979 
# F-statistic: 15.41 on 1 and 129 DF,  p-value: 0.00014 


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
plot(density(scores[,1]))# mostly, a little negative skew
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # bit of positive skew

# Assign scores for pc1, 2, & 3 to thier observation (rows)
varframe.cmplt.09 <- (soilTData$yearTsoil==2009 & complete.cases(varframe))
soilTData$pc1score <- NA
soilTData$pc2score <- NA
soilTData$pc3score <- NA
soilTData$pc1score[varframe.cmplt.09] <- scores[,1]
soilTData$pc2score[varframe.cmplt.09] <- scores[,2]
soilTData$pc3score[varframe.cmplt.09] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score, data=soilTData)
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.54842    0.05822   9.420  < 2e-16 ***
# pc1score    -0.01330    0.01695  -0.784 0.434087    
# pc2score     0.15700    0.02329   6.742 3.71e-10 ***
# pc3score    -0.14195    0.03836  -3.700 0.000308 ***
# Multiple R-squared: 0.2977,	Adjusted R-squared: 0.2827 
# F-statistic: 19.92 on 3 and 141 DF,  p-value: 8.012e-11 

lm2 <- lm(snowcovTs20mean ~ pc2score + pc3score, data=soilTData)
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.54842    0.05814   9.432  < 2e-16 ***
# pc2score     0.15700    0.02326   6.751 3.47e-10 ***
# pc3score    -0.14195    0.03831  -3.706 0.000301 ***
# Multiple R-squared: 0.2946,	Adjusted R-squared: 0.2847 
# F-statistic: 29.65 on 2 and 142 DF,  p-value: 1.732e-11 

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
qqline(scores[,3]) # yes except for at the extremes, some neg. skew

# Assign scores for pc1, 2, & 3 to thier observation (rows)
varframe.cmplt.11 <- (soilTData$yearTsoil==2011 & complete.cases(varframe))
soilTData$pc1score <- NA
soilTData$pc2score <- NA
soilTData$pc3score <- NA
soilTData$pc1score[varframe.cmplt.11] <- scores[,1]
soilTData$pc2score[varframe.cmplt.11] <- scores[,2]
soilTData$pc3score[varframe.cmplt.11] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score, data=soilTData)
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.87709    0.03679  23.838  < 2e-16 ***
# pc1score     0.03494    0.01052   3.321   0.0011 ** 
# pc2score     0.10962    0.01486   7.376  7.4e-12 ***
# pc3score    -0.04965    0.02299  -2.160   0.0322 *  
# Multiple R-squared: 0.2969,	Adjusted R-squared: 0.2842 
# F-statistic: 23.37 on 3 and 166 DF,  p-value: 1.145e-12

