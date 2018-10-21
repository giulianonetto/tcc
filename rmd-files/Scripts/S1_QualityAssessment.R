library(madsim)
library(vsn)
library(affy)
library
system("mkdir ExampleMicroarrayData")
fparams <- data.frame(m1 = 7, m2 = 7, shape2 = 4, lb = 4, ub =
                        14, pde = 0.02, sym = 0.5)
dparams <- data.frame(lambda1 = 0.13, lambda2 = 2, muminde = 1,
                      sdde = 0.5)
sdn <- 0.4
rseed <- 50
n <- 35000
expdata1 <- madsim(mdata=NULL, n=35000, ratio=0, fparams, dparams,
                 sdn, rseed)
expdata <- normalizeBetweenArrays(expdata[[1]])
# Get Exp Data in the long format for boxplot and histograms/density plots
write.table(expdata, "ExampleMicroarrayData/expressionData.tsv", sep = "\t", quote = F)
system("sed -i s/,/./g ExampleMicroarrayData/expressionData.tsv")

