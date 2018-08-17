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
mydir <- '/home/giulianonetto/windows/tcc/storage/'
datadir <- "~/windows/tcc/Datasets"
setwd(datadir)
system("find . -type d -name '*GSE*' > dataset_dirs.txt")
system("mv dataset_dirs.txt /home/giulianonetto/windows/tcc/storage/")
setwd(mydir)
system("grep -o 'GSE[ ]*[0-9]\\+' dataset_dirs.txt > geo.txt")
system("sed -i 's/ //g' geo.txt")
system("grep -o '[a-zA-Z]\\+, [0-9]\\+' dataset_dirs.txt > authors.txt")
system("grep -o -E '[0-9]{4}' authors.txt > years.txt")
(c1 <- read.csv('authors.txt')[,1])
(c2 <- read.csv('years.txt')[,1])
(c3 <- read.csv('geo.txt')[,1])
datasets <- data.frame(c1, c2, c3)
datasets <- datasets[order(datasets$c2),]
system("rm dataset_dirs.txt")
system("rm geo.txt")
system("rm authors.txt")
system("rm years.txt")
write.table(datasets, file="datasets.tsv", sep = "\t", row.names = FALSE, col.names = c("author", "year", "geo"))

system("cat datasets.tsv")
datasets <- read.csv("datasets.tsv", sep = "\t", header = T, stringsAsFactors = F)
plot(table(datasets$year), main="Datasets by year", ylab = "# of Datasets", ylim = c(0,15))
