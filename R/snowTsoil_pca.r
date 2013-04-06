# snowTsoil_pca
# This is the script for performing a pca and further inferential tests (using
# pca results) on below-snow Tsoil across the SNOTEL network. This pca is 
# designed for testing the predictors of mean below-snow Tsoil at 20cm depth.

source('getdata.r')

# climData.sub should have the same sites/years in the same order as soilTData
# and soilVWCData
# sum(soilVWCData[,1]-climData.sub[,1] )

# Create a dataframe of variables for pca - include monthly means
varframe <- data.frame(climData.sub$elev,climData.sub$totaldaysSC,
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
# We'll use any pc with SD > 1
var_exp <- summary(ts20.pca)$importance[,1:4]

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

# Put them in a matrix
pc_loadings <- loadings[,1:4]

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

# Put them in a matrix
lm_coeffs <- summary(lm1)$coefficients

summary(lm2)

summary(lm3) # not showing site

summary(lm4) # not showing site or year

# One thing that is unclear here is that year seems to have a big effect.
# It might be best to try this with some selected years, then exclude year
# from the multiple regression model.
#
# 2007
#
varframe.07 <- varframe[climData.sub$yearClim==2007,] # Boolean 2007 rows
varframe.cmplt <- complete.cases(varframe.07)
ts20.07.pca <- prcomp(data.matrix(varframe.07[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(ts20.07.pca)

# We'll use any pc with SD > 1
var_exp <- cbind(var_exp, summary(ts20.07.pca)$importance[,1:4])

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

# Put them in a matrix
pc_loadings <- cbind(pc_loadings, loadings[,1:4])

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

# Put them in a matrix
lm_coeffs <- cbind(lm_coeffs, summary(lm1)$coefficients)

lm2 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score, data=soilTData)
AIC(lm2) # 260.3205 
summary(lm2)

# 2009
#
varframe.09 <- varframe[climData.sub$yearClim==2009,] # Boolean 2009 rows
varframe.cmplt <- complete.cases(varframe.09)
ts20.09.pca <- prcomp(data.matrix(varframe.09[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(ts20.09.pca)

# We'll use any pc with SD > 1
var_exp <- cbind(var_exp, summary(ts20.09.pca)$importance[,1:4])

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

# Put them in a matrix
pc_loadings <- cbind(pc_loadings, loadings[,1:4])

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

# Put them in a matrix
lm_coeffs <- cbind(lm_coeffs, summary(lm1)$coefficients)

lm2 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score, data=soilTData)
AIC(lm2) # 314.4483
summary(lm2)

# 2011
#
varframe.11 <- varframe[climData.sub$yearClim==2011,] # Boolean 2011 rows
varframe.cmplt <- complete.cases(varframe.11)
ts20.11.pca <- prcomp(data.matrix(varframe.11[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(ts20.11.pca)

# We'll use any pc with SD > 1
var_exp <- cbind(var_exp, summary(ts20.11.pca)$importance[,1:4])

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

# Put them in a matrix
pc_loadings <- cbind(pc_loadings, loadings[,1:4])

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

# Put them in a matrix
lm_coeffs <- cbind(lm_coeffs, summary(lm1)$coefficients)

lm2 <- lm(snowcovTs20mean ~ pc1score+pc2score+pc3score, data=soilTData)
AIC(lm2) # 238.6498
summary(lm2)


# DONE - now show the tables we created

# The variance explained:
var_exp

# The loadings
pc_loadings

# The regression coefficients
lm_coeffs
