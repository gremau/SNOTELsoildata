# gsVWC_pca
# This is the script for performing a pca and further inferential tests (using
# pca results) on growing season soil VWC across the SNOTEL network. This pca 
# is designed for testing the predictors of mean growing season soil VWC at 
# 20cm depth.
setwd('/home/greg/data/current/SNOTELsoil-climate/data_analysis/R/')
source('getdata.r')
library(xtable)

# climData.sub should have the same sites/years in the same order as soilVWCData
# and soilVWCData
# sum(soilVWCData[,1]-climData.sub[,1] )

# Create a dataframe of variables for pca - include monthly means
varframe <- data.frame(climData.sub$elev,climData.sub$totaldaysSC,
                       climData.sub$meltdoy,climData.sub$maxswe,
                       climData.sub$scovmat,climData.sub$jasMAT,
                       climData.sub$aprTairmean,climData.sub$mayTairmean,
                       climData.sub$junTairmean,climData.sub$julTairmean,
                       climData.sub$augTairmean,climData.sub$sepTairmean,
                       climData.sub$JASprecip,climData.sub$mayPrecip,
                       climData.sub$junPrecip,climData.sub$julPrecip,
                       climData.sub$augPrecip,climData.sub$sepPrecip,
                       soilTData$jfmTs5mean,soilVWCData$jfmVWC20mean)


# Shorten the variable names in varframe
names(varframe) = sub("climData.sub.","",names(varframe))
names(varframe) = sub("soilTData.","",names(varframe))
names(varframe) = sub("soilVWCData.","",names(varframe))

# We are going to make dataframes that will become a table later
varexp.table <- data.frame(PC1all=numeric(), PC107=numeric(), 
                           PC109=numeric(), PC111=numeric(), PC2all=numeric(),
                           P207=numeric(), PC209=numeric(), PC211=numeric(),
                           PC3all=numeric(), PC307=numeric(), PC309=numeric(),
                           PC311=numeric())

loadings.table <- data.frame(PC1all=numeric(), PC107=numeric(), 
                             PC109=numeric(), PC111=numeric(), PC2all=numeric(),
                             P207=numeric(), PC209=numeric(), PC211=numeric(),
                             PC3all=numeric(), PC307=numeric(), PC309=numeric(),
                             PC311=numeric())

yearcount <- 0 # a year counter to set the correct column
colref <- c(1,5,9)

addtotable <- function(tab, toadd) {
  colref <- colref + yearcount
  tab[1:nrow(toadd),colref] <- toadd[,1:ncol(toadd)]
  return(tab)
}

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

# We'll use any pc with SD > 1 (dropped 4 after further analysis)
var_exp <- summary(vwc20.pca)$importance[,1:3]
# Put them in the table (don't increment year yet)
varexp.table <- addtotable(varexp.table, var_exp[,1:3])

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

# Put them in the table and increment year
loadings.table <- addtotable(loadings.table, loadings[,1:3])
yearcount <- yearcount + 1

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
qqline(scores[,2]) # yes, slight positive skew
qqnorm(scores[,3])
qqline(scores[,3]) # yes, slight negative skew
qqnorm(scores[,4])
qqline(scores[,4]) # yes

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
lm1 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilVWCData)
AIC(lm1) # -881.0855

# remove PC4 
lm2 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score, data=soilVWCData)
AIC(lm2) # -876.2319 - so it seems like pc4 is important

# Add site
lm3 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score+pc4score
	 + as.factor(siteVWC), data=soilVWCData)
AIC(lm3) # -2098.235

# Add year
lm4 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score+pc4score
	  + as.factor(yearVWC) + as.factor(siteVWC), data=soilVWCData)
AIC(lm4) # -2179.495

summary(lm1)

# Put them in a matrix
lm_coeffs <- summary(lm2)$coefficients
lm_coeffs <- rbind(lm_coeffs[,c(1,2,4)],
                   c(NA, summary(lm2)$adj.r.squared, NA))

summary(lm2)

summary(lm3)

summary(lm4)


# One thing that is unclear here is that year seems to have a big effect.
# It might be best to try this with some selected years, then exclude year
# from the multiple regression model.
#
# 2007
#
varframe.07 <- varframe[climData.sub$yearClim==2007,] # Boolean 2007 rows
varframe.cmplt <- complete.cases(varframe.07)
vwc20.07.pca <- prcomp(data.matrix(varframe.07[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(vwc20.07.pca)

# We'll use any pc with SD > 1
var_exp <- summary(vwc20.07.pca)$importance[,1:3]
# Put them in the table (don't increment year yet)
varexp.table <- addtotable(varexp.table, var_exp[,1:3])

# Make a screeplot of this
screeplot(vwc20.07.pca)

# Now make some useful variables
sd <- vwc20.07.pca$sdev
loadings <- vwc20.07.pca$rotation
rownames(loadings) <- colnames(varframe.07)
scores <- vwc20.07.pca$x
rownames(scores) <- 1:dim(scores)[1]

# Look at loadings
loadings[,1:3]

# Put them in the table and increment year
loadings.table <- addtotable(loadings.table, loadings[,1:3])
yearcount <- yearcount + 1

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
qqline(scores[,2]) # mostly, tiny skew
qqnorm(scores[,3])
qqline(scores[,3]) # yes except for at the extremes
qqnorm(scores[,4])
qqline(scores[,4]) # mostly, bit of negative skew

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
lm1 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilVWCData)
AIC(lm1) # -153.8788
summary(lm1)

lm2 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score, data=soilVWCData)
AIC(lm2) # -153.1189
summary(lm2)

# Put them in a matrix
lm_coeffs <- cbind(lm_coeffs, rbind(summary(lm2)$coefficients[,c(1,2,4)],
                                    c(NA, summary(lm2)$adj.r.squared, NA)))

# 2009
#
varframe.09 <- varframe[climData.sub$yearClim==2009,] # Boolean 2009 rows
varframe.cmplt <- complete.cases(varframe.09)
vwc20.09.pca <- prcomp(data.matrix(varframe.09[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(vwc20.09.pca)

# We'll use any pc with SD > 1
var_exp <- summary(vwc20.09.pca)$importance[,1:3]
# Put them in the table (don't increment year yet)
varexp.table <- addtotable(varexp.table, var_exp[,1:3])

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

# Put them in the table and increment year
loadings.table <- addtotable(loadings.table, loadings[,1:3])
yearcount <- yearcount + 1

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
plot(density(scores[,1]))# yes, except for a little at the extremes
qqnorm(scores[,2])
qqline(scores[,2]) # yes, tiny positive skew
qqnorm(scores[,3])
qqline(scores[,3]) # yes
qqnorm(scores[,4])
qqline(scores[,4]) # some positive skew 

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
lm1 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilVWCData)
AIC(lm1) # -161.7993
summary(lm1)

lm2 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score, data=soilVWCData)
AIC(lm2) # -153.3085
summary(lm2)

# Put them in a matrix
lm_coeffs <- cbind(lm_coeffs, rbind(summary(lm2)$coefficients[,c(1,2,4)],
                                    c(NA, summary(lm2)$adj.r.squared, NA)))

# 2011
#
varframe.11 <- varframe[climData.sub$yearClim==2011,] # Boolean 2011 rows
varframe.cmplt <- complete.cases(varframe.11)
vwc20.11.pca <- prcomp(data.matrix(varframe.11[varframe.cmplt,]), 
		   center=TRUE, scale=TRUE, retx=TRUE)

# Get a summary of the variance explained by each principal component
summary(vwc20.11.pca)

# We'll use any pc with SD > 1
var_exp <- summary(vwc20.11.pca)$importance[,1:3]
# Put them in the table (don't increment year)
varexp.table <- addtotable(varexp.table, var_exp[,1:3])

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

# Put them in the table and increment year
loadings.table <- addtotable(loadings.table, loadings[,1:3])
yearcount <- yearcount + 1

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
qqline(scores[,3]) # yes except for at the extremes
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

# Now it should be possible to do a multiple regression model with the 4
# principal components. Note that using site as a predictor results in
# all singularities (invalid model)
lm1 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score+pc4score, data=soilVWCData)
AIC(lm1) # -128.7211
summary(lm1)

lm2 <- lm(jasVWC20mean ~ pc1score+pc2score+pc3score, data=soilVWCData)
AIC(lm2) # -115.6261 
summary(lm2)

# Put them in a matrix
lm_coeffs <- cbind(lm_coeffs, rbind(summary(lm2)$coefficients[,c(1,2,4)],
                                    c(NA, summary(lm2)$adj.r.squared, NA)))

# DONE - now show the tables we created

# The variance explained:
varexp.table
rownames(varexp.table) <- c('Std. Deviation', '% Var. Explained',
                            'Cum. Var. Explained')
print(xtable(varexp.table, floating=T), file='../tables/rawtableB5.tex')
print(xtable(varexp.table, floating=T),
      file='../../manuscript_1/tables/rawtableB5.tex')

# The loadings
loadings.table
rownames(loadings.table) <- c('Elevation', 'Total snow-covered days', 'Snow-free day',
                            'Peak SWE', 'Below-snow period T\textsubscript{air}',
                            'Summer quarter mean T\textsubscript{air}',
                            'Apr. Mean T\textsubscript{air}',
                            'May. Mean T\textsubscript{air}',
                            'Jun. Mean T\textsubscript{air}',
                            'Jul. Mean T\textsubscript{air}',
                            'Aug. Mean T\textsubscript{air}',
                            'Sep. Mean T\textsubscript{air}',
                            'Summer quarter precip.',
                            'May precip.', 'Jun. precip.', 'Jul. precip.',
                            'Aug. precip.', 'Sep. precip.',
                            'Winter quarter 5 cm T\textsubscript{soil}',
                            'Winter quarter 20 cm θ')

print(xtable(loadings.table, floating=T), file='../tables/rawtableB6.tex')
print(xtable(loadings.table, floating=T),
      file='../../manuscript_1/tables/rawtableB6.tex')

# The regression coefficients
lm_coeffs
rownames(lm_coeffs) <- c('(Intercept)', 'PC 1', 'PC 2', 'PC 3',
                         'Model adj. R\textsuperscript{2}')
print(xtable(lm_coeffs, floating=T, digits=3), file='../tables/rawtableB7.tex')
print(xtable(lm_coeffs, floating=T, digits=3),
      file='../../manuscript_1/tables/rawtableB7.tex')