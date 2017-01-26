library(devtools)
library(Biobase)
install.packages(c("devtools"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase"))

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file= con)
close(con)
mp = montpick.eset
pdata = pData(mp)
edata = as.data.frame(exprs(mp))
fdata = fData(mp)
ls()
edata = edata[rowMeans(edata) > 100,] # to make it easier to visualize, substract out all of the rowmeans value which is less than 100, so reduce the size of data set 
dim(edata)
edata = log2(edata +1) # apply log2 transform so it is easier to work with 
# next thing to do is center the data because when we are doing singular decomposition, if we dont centre the data, means we dton remove the row means of the data, then the first singular value of your vector will be the mean value 
edata_centered = edata - rowMeans(edata)
# svd function to calculate the singular value decomposition 
svd1= svd(edata_centered)
names(svd1)
# singular value decomposition has 3 matriced : d,u,v ( diagonal matrix, 
svd1$d
svd1$u
dim(edata)

dim(svd1$u)
dim(svd1$v)
# plot the singular values 
plot(svd1$d,ylab="Singular value",col=2)
par ( pch =19)
plot(svd1$d,ylab="Singular value",col=2)
plot(svd1$d^2/sum(svd1$d^2),ylab="Percent Variance Explained",col=2)
plot (svd1$d^2/sum(svd1$d^2), ylab = "Percent Variance Explained", col = 2)
par(mfrow=c(1,2))
#Plot top two principal components
plot(svd1$v[,1], col = 2, ylab = "1st PC")
plot(svd1$v[,2], col = 3, ylab = "2nd PC")
plot(svd1$v[,1],svd1$v[,2],col=2,ylab="2nd PC",xlab="1st PC")
par(mfrow=c(1,1))
plot(svd1$v[,1],svd1$v[,2],col=2,ylab="2nd PC",xlab="1st PC")
plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC", xlab="1st PC",col=as.numeric(pdata$study))
# plot the 2 singular values ( principal components) as the color of pdata study 
# now tow studies has different variance 
boxplot (svd1$v[,1] ~ pdata$study, border = c(1,2))
# overlay the data point on top of the boxplot by using point function, plotting the same sigular vector vesus a jittered verson of  the study variable. 
points 
points(svd1$v[,1] ~ jitter(as.numeric(pdata$study)), col = as.numeric(pdata$study))
# jitter (base) as noise to numbers, as a small amount of noise to a numeric vector, should be numeric vector, as.numeric (pdata$study)
# now you can see there is a big difference in the first principal component between these two studies. 
pc1= prcomp(edata)
plot (pc1$rotation[,1], svd1$v[,1]) # now if we plot the first PC1 vesus the first singular vector, they are not the same since we have not scaled them in the same way
# center the col means, substracts the col means  rather than subtracting the row means, then we have the data set centered by col  instead of centered by row
edata_centered2 = t(t(edata) - colMeans(edata))
svd2 = svd(edata_centered2)
plot(pc1$rotation[,1],svd2$v[,1],col=2)
# why subtracting by column, because the principal component are calculating the variability between the columns 
# What happens if we introduce a single outlying gene
edata_outlier = edata_centered
edata_outlier[6,] = edata_centered[6,] * 10000 # take six genes, and make it really outlined
svd3 = svd(edata_outlier)
par(mfrow=c(1,2))
plot(svd1$v[,1],col=1,main="Without outlier")
plot(svd3$v[,1],col=2,main="With outlier")

plot ( svd1$v[,1], svd3$v[,1], xlab = "without outlier", ylab = " with Outlier")
# two data set with and without outlier dont neccessarily match to each other in term of  their value decomposition
plot(svd3$v[,6],edata_outlier[6,],col=4)
edata_outlier[6,] = edata_centered[6,] * 10000
plot(svd3$v[,1],edata_outlier[6,],col=4)
# plot the first singular value decomposition of the new decomposition with the outlying values itself you can see the highly corelated each other 
# the decomposition is looking for the patterns of variation, one gene is way higher expressed, then it's going to drive the most variation in the data set in the same way, so it'll be corelated 
# need to be careful when using these decomposition, make sure that you pick the centering and scaling so all of the deffierent measurements for all of the features  are on the common scale 