# snowTsoil_pca
# This is the script for performing a pca and further inferential tests (using
# pca results) on below-snow Tsoil across the SNOTEL network. This pca is 
# designed for testing the predictors of mean below-snow Tsoil at 20cm depth.

source('getdata.r')

# climData.sub should have the same sites/years in the same order as soilTData
# and soilVWCData
# sum(soilVWCData[,1]-climData.sub[,1] )

# Create a dataframe of variables for pca
varframe <- data.frame(climData.sub$totaldaysSC, climData.sub$maat, 
		climData.sub$meltdoy,climData.sub$onsetdoy,
		climData.sub$maxswe, climData.sub$scovmat,
		soilVWCData$preonsetVWC20,soilTData$preonsetTs20,
		climData.sub$preonsetTair,climData.sub$elev)

# Create a dataframe of variables for pca - include monthly means
varframe <- data.frame(climData.sub$totaldaysSC, #climData.sub$maat, 
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
		climData.sub$preonsetTair,climData.sub$elev)


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
biplot(scores[,1:2], loadings[,1:2], xlab=rownames(scores),
       ylab=rownames(loadings), cex=0.7)

# Are axes 1, 2, and 3 normal?
qqnorm(scores[,1])
qqline(scores[,1])
density(scores[,1])# mostly, a little negative skew
qqnorm(scores[,2])
qqline(scores[,2]) # yes
qqnorm(scores[,3])
qqline(scores[,3]) # yes

# Assign scores for pc1, 2, & 3 to thier observation (rows)
soilTData$pc1score[varframe.cmplt] <- scores[,1]
soilTData$pc2score[varframe.cmplt] <- scores[,2]
soilTData$pc3score[varframe.cmplt] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components
lm1 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score, data=soilTData)

lm2 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score
	 + as.factor(siteTsoil), data=soilTData)

lm3 <- lm(snowcovTs20mean ~ pc1score + pc2score, data=soilTData)

lm4 <- lm(snowcovTs20mean ~ pc1score + pc2score
	 + as.factor(siteTsoil), data=soilTData)

lm5 <- lm(snowcovTs20mean ~ pc1score + pc2score + as.factor(yearTsoil)
	 + as.factor(siteTsoil) , data=soilTData)

lm6 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score
	  + as.factor(yearTsoil) + as.factor(siteTsoil), data=soilTData)

AIC(lm1) # 2149.277
AIC(lm2) # 1600.425
AIC(lm3) # 2150.76
AIC(lm4) # 1598.536
AIC(lm5) # 1269.112
AIC(lm6) # 1236.874

summary(lm4)
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               0.612258   0.290035   2.111 0.035093 *  
# pc1score                 -0.034856   0.011940  -2.919 0.003610 ** 
# pc2score                 -0.122687   0.017712  -6.927 9.04e-12 ***
# Multiple R-squared: 0.6923,	Adjusted R-squared: 0.612 
# F-statistic: 8.613 on 203 and 777 DF,  p-value: < 2.2e-16 

summary(lm5)
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               4.328343   0.554666   7.804 1.98e-14 ***
# pc1score                 -0.102177   0.018624  -5.486 5.58e-08 ***
# pc2score                  0.017746   0.025999   0.683 0.495095
# Multiple R-squared: 0.7858,	Adjusted R-squared: 0.7253 
# F-statistic: 12.98 on 216 and 764 DF,  p-value: < 2.2e-16

summary(lm6)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               4.0248562  0.5485341   7.337 5.58e-13 ***
# pc1score                 -0.1152548  0.0184850  -6.235 7.47e-10 ***
# pc2score                  0.0310184  0.0256931   1.207 0.227703    
# pc3score                 -0.0624008  0.0119869  -5.206 2.49e-07 ***
# Multiple R-squared: 0.7932,	Adjusted R-squared: 0.7344 
# F-statistic: 13.49 on 217 and 763 DF,  p-value: < 2.2e-16 


# One thing that is unclear here is that year seems to have a big effect.
# It might be best to try this with some selected years, then exclude year
# from the multiple regression model.
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
qqline(scores[,3]) # yes except for at the extremes

# Assign scores for pc1, 2, & 3 to thier observation (rows)
varframe.cmplt.09 <- (soilTData$yearTsoil==2009 & complete.cases(varframe))
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
qqline(scores[,3]) # yes except for at the extremes

# Assign scores for pc1, 2, & 3 to thier observation (rows)
varframe.cmplt.11 <- (soilTData$yearTsoil==2011 & complete.cases(varframe))
soilTData$pc1score[varframe.cmplt.11] <- scores[,1]
soilTData$pc2score[varframe.cmplt.11] <- scores[,2]
soilTData$pc3score[varframe.cmplt.11] <- scores[,3]

# Now it should be possible to do a multiple regression model with the 3
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(snowcovTs20mean ~ pc1score + pc2score + pc3score, data=soilTData)
summary(lm1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.73004    0.03526  20.703  < 2e-16 ***
# pc1score     0.01301    0.01017   1.280    0.202    
# pc2score     0.12919    0.01417   9.118  < 2e-16 ***
# pc3score    -0.09523    0.02276  -4.185 3.73e-05 ***
# Multiple R-squared: 0.2499,	Adjusted R-squared: 0.2426 
# F -statistic:  34.1 on 3 and 307 DF,  p-value: < 2.2e-16 

lm2 <- lm(snowcovTs20mean ~ pc2score + pc3score, data=soilTData)
summary(lm2)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.73004    0.03530  20.682  < 2e-16 ***
# pc2score     0.12919    0.01418   9.108  < 2e-16 ***
# pc3score    -0.09523    0.02278  -4.180  3.8e-05 ***
# Multiple R-squared: 0.2459,	Adjusted R-squared: 0.241 
# F-statistic: 50.22 on 2 and 308 DF,  p-value: < 2.2e-16 

