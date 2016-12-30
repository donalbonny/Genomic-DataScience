tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
# tell R that I want to use the set of color for my plots, using pallete command
par(pch=19) # the circles are filled in, rather than open circles when making plots.

library(gplots)
library(devtools)
library(Biobase)
library(RSkittleBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
source("http://bioconductor.org/biocLite.R")
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
biocLite(c("Biobase","org.Hs.eg.db","AnnotationDbi"))
biocLite("alyssafrazee/RSkittleBrewer")
library(RSkittleBrewer)
trop = RSkittleBrewer("tropical")

con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")

load(file = con)


close(con)
ls()
bm = bodymap.eset
pdata= pData(bm)
edata= exprs(bm)
fdata= fData(bm)
# now we have new variables that we have defined 
ls()
dim(pdata)
head(pdata, 5)
head(edata)
head(fdata, 5)
table(pdata$gender)
table(pdata$gender, pdata$race)
summary(edata) # show the summary of every column in expression data
table(pdata$age, useNA = 'ifany') # to show NA value in the table 
sum(pdata$age == " ") # meaning there are some NA values in the age variable 
sum(pdata$age == " ", na.rm = TRUE) # remove any NA values when doing this sum 
is.na(edata) [1,] # check any missing value in edata, its FALSE if there is not missing, its TRUE if there is missing value 
sum (is.na(edata)) # sum up all of the missing value data in the data set 
gene_na = rowSums(is.na(edata))
table(gene_na)
sample_na = colSums(is.na(edata))
table (sample_na)
# the number of rows in the phenotype data should match the number of columns in genomic data
# the number of rows in feature data should match the number of rows in the expression data since feature data describes gene, and phenotype data describes samples. 
# check all the missing values and dimention match 
boxplot(edata[, 1]) # boxplot of edata and first colum 
# so now the most of values are down at zero so we can do the transform to the data , take a log transform of the data 
boxplot(log2(edata[, 1] + 1))
boxplot(log2(edata+1))
par(mfrow=c(1,2)) # set up onw row and 2 columns of plots, set up side by side plots
hist(log2(edata [,1] +1), col = 2) # histogram of all the value in edata, almost the values are equal to zero 

par(mfrow=c(1,1))
plot(density(log2(edata[,1]+1)), col=2)
lines(density(log2(edata[,2] +1)), col = 3)
lines(density(log2(edata[,3] +1)), col = 4)
qqplot(log2(edata[,1] +1), log2(edata[,2] +1), col= 4) # all of the pcercentile of samples are displayed by the dots 
abline(c(0,1)) # add the 45 degree line into the plot
# qq plot to see how the distribution compared to each other
# MA plot when you want to see especially for the technical replicates.
mm = log2(edata[,1] +1) - log2(edata[,2] +1)
aa = log2(edata[,1] +1) + log2(edata[,2] +1)
plot (aa, mm , col=2)
# X axis: aa is the sum of the samples, Y axis you take the difference between two samples. 
 edata = as.data.frame(edata) # make the dataset as the data frame so that you can use dplyr  filter command 
 filt_edata = filter(edata, rowMeans(edata) > 1) # filter the mean of each row, if the mean is greater than 1, it is going to keep it, otherwise
 dim (filt_edata)
 
 boxplot (as.matrix(log2(filt_edata +1)) , col=2 )
 