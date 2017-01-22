library(devtools)
library(Biobase)
library(preprocessCore)
install.packages(c("devtools"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase","preprocessCore"))
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
head(edata)
head(pdata)
head(fdata)
fdata = fData(mp)
ls()
edata = log2(edata +1)
edata = edata[rowMeans(edata) >3,]
dim(edata) # edata has 5862 genes and 129 samples 
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:20){lines(density(edata[,i]),lwd=3,col=colramp[i])}
# make the loop from 2:20, make the lines of the other samples on top of the first sample
# now some of the samples have the similar distribution but some of the samples have big ditributional differences between samples , that likely due to the technology, not due to biology 
# how di distinguish that the distributional differences come from the technical issue , not from biology issue
norm_edata = normalize.quantiles(as.matrix(edata))
dim(edata)
dim(norm_edata)
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:20){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}
plot(norm_edata[1,],col=as.numeric(pdata$study))
par(pch = 19)
svd1 = svd(norm_edata - rowMeans(norm_edata))
plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",col=as.numeric(pdata$study))
