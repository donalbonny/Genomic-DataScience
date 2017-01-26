library(Biobase)
library(GenomicRanges)
data(sample.ExpressionSet, package = "Biobase")

se = makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)


library(plotrix)
source("http://bioconductor.org/biocLite.R")
biocLite("plottrix")
library(plotrix)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)

library(plotrix)
pie3D(pdata_bm$num.tech.reps,labels=pdata_bm$tissue.type)
row_sums = rowSums(edata)
edata = edata[order(-row_sums),]
index = 1:500
index
heatmap(edata[index,],Rowv=NA,Colv=NA)


row_sums = rowSums(edata)
index = which(rank(row_sums) < 500 )
heatmap(edata[index,],Colv=NA)
row_sums = rowSums(edata)
row_sums
edata = edata[order(row_sums),]
index = which(rank(-row_sums) < 500 )
index
heatmap(edata[index,],Rowv=NA,Colv=NA)
heatmap(edata[index,],Rowv=NA,Colv=NA)

row_sums = rowSums(edata)
index = which(rank(-row_sums) < 500 )
heatmap(edata[index,],Rowv=NA,Colv=NA)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata = pData(bm)
edata = exprs(bm)
mm = log2(edata[,1] +1) - log2(edata[,2] +1)
aa = log2(edata[,1] +1) + log2(edata[,2] +1)
plot (aa, mm , col=2)
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
par(mfrow=c(1,2))
mm = log2(edata[,1] +1) - log2(edata[,2] +1)
aa = log2(edata[,1] +1) + log2(edata[,2] +1)
plot (aa, mm , col=2)

help("rlog")
??rlog

rlogTransformation(as.matrix(edata), blind = TRUE, intercept, betaPriorVar, fitType = "parametric")
(edata, blind = TRUE, intercept, betaPriorVar, fitType = "parametric")
mm1 = log2(edata[,1] +1) - log2(edata[,2] +1)
aa1 = log2(edata[,1] +1) + log2(edata[,2] +1)
plot (aa, mm , col=2)



con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
mp
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = edata[rowMeans(edata) <100,]
dim(edata)
dist1 = dist(t(edata))
hclust1 = hclust(dist1)
plot(hclust1)
edata = log2(edata + 1)
dist1 = dist(t(edata))
hclust1 = hclust(dist1)
edata = edata[rowMeans(edata) <100,]
dist1 = dist(t(edata))
dim(dist1)
dim(edata)
dist1 = dist(t(edata))
hclust1 = hclust(dist1)
edata = log2(edata + 1)
dist1 = dist(t(edata))
hclust1 = hclust(dist1)
?cutree
