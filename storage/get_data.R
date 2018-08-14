############################# DATA GATHERING FOR TCC ##################################

########## 08/12/2018

##### REQUIREMENTS  #####
## SEARCH - GEO Datasets
"pulmonary fibrosis"
"bleomycin IPF"
"bleomycin pulmonary fibrosis"
"time course pulmonary fibrosis"
"silica pulmonary fibrosis"

# R version 3.4.1 (2017-06-30) - Ubuntu
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)

mydir <- getwd()
datadir <- "~/windows/tcc/Datasets"
setwd(datadir)
system("find . -type d -name '*GSE*' > list1.txt")
system("grep -o 'GSE[ ]*[0-9]\\+' list1.txt > list.txt")
system("rm list1.txt")
system("sed 's/ //g' list.txt")

