library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)
install.packages(c("devtools"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase","sva","bladderbatch","snpStats"))
data(bladderdata)
a
y
ls()
pheno = pData(bladderEset)
edata = exprs(bladderEset)
dim(edata)
dim(pheno)
head(pheno)

# build a matrix model  between the factor cancer and batch variables
mod = model.matrix(~as.factor(cancer) + as.factor(batch), data = pheno)
fit = lm.fit(mod, t(edata))
hist(fit$coefficients[2,], col = 2, breaks = 100)

table (pheno$cancer, pheno$batch)
