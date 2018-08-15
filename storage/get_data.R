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
system("find . -type d -name '*GSE*' > list.txt")
system("grep -o 'GSE[ ]*[0-9]\\+' list.txt > list1.txt")
system("sed -i 's/ //g' list1.txt")
system("grep -o '[A-Z][a-zA-Z]\\+, [0-9]\\+[a-z]*' list.txt > list2.txt")
system("cat list.txt")
system("cat list1.txt")
system("cat list2.txt")
c1 <- as.vector(read.csv('list1.txt', header = FALSE)[,1])
c2 <- as.vector(read.csv('list2.txt', header = FALSE, sep = '\t')[,1])
df <- data.frame(c1, c2)
write.table(df, file="datasets.tsv", sep = "\t", row.names = FALSE, col.names = c("geo", "paper"))

