# Data Transforms 
library(devtools)
library(Biobase)
install.packages(c("devtools"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase"))
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()
hist(rnorm(1000), col=2)
hist(edata[,1], col=2)
# one of the transformation is putting things on log scale. easier for visualize the data 
hist(log(edata[,1]), col=2, breaks= 1000)
# value can be -inf because you have mostly zeros in the data set
# add small number add one, log transformation defined for the rest of them
hist(log(edata[,1] +1), col=2, breaks= 1000)
hist(log2(edata[,1] +1), col=2, breaks= 100)
#log 2 transform can be used when you want to compare the 2 data set
# zoom in using the xlim and ylim value c(a,b)
hist(log2(edata[,1] +1), col=2, breaks= 100, xlim = c(1, 15), ylim = c(0,400) )
hist(rowSums(edata==0), col=2 ) # there are a lot of genes that bacsically have all zeros value 
low_genes = rowMeans(edata)<5
table(low_genes)
filt_edata = filter(as.data.frame(edata), !low_genes) # filter edata set, not keeping the low_genes 
dim(filt_edata)
summary(edata)
low_genes2 = rowMedians(as.matrix(edata)) <5
table(low_genes2)

# transform and filtering out some low count observation

