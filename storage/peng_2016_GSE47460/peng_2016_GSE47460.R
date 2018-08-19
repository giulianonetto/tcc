library(dplyr)
library(GEOquery)
mydir <- getwd()
datasets <- read.csv("datasets.tsv", sep = "\t", header = T, stringsAsFactors = F)
geo <- datasets %>% filter(author == "Peng", year == 2016) %>% select(geo) %>% unlist()
system("mkdir peng_2016_GSE47460")
analysis_dir <- "peng_2016_GSE47460/"
setwd(analysis_dir)

filePaths <- getGEOSuppFiles(GEO = geo)
