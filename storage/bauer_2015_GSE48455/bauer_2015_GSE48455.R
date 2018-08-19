########### DOWNLOADS ############
library(dplyr)
library(GEOquery)
mydir <- getwd()
datasets <- read.csv("../datasets.tsv", sep = "\t", header = T, stringsAsFactors = F)
geo <- (datasets %>% filter(author == "Bauer", year == 2015) %>% select(geo) %>% unlist())[2]
system("mkdir bauer_2015_GSE48455")
analysis_dir <- "bauer_2015_GSE48455/"
setwd(analysis_dir)

#filePaths <- getGEOSuppFiles(GEO = geo)
data <- getGEO(geo)[[1]]

# vns package for normalizaton and glog transformation

########### METADATA ############
library(stringr)
metadata <- pData(data)
sample_names <- metadata$geo_accession
treatment <- metadata$`treatment:ch1`
times <- metadata$`time:ch1`
metadata <- data.frame(sample_names, treatment, times)
D <- vector()
N <- vector()
times_clean <- vector()
for (i in times){
  d <- str_extract(i, '[A-Z]{4}')
  n <- as.numeric(str_extract(i, '[0-9]'))
  if(d == "WEEK"){
    clean <- 7
  } else if(d == "DAYS"){
    clean <- 1
  } else {
    print(d, n)
  }
  times_clean <- c(times_clean, n*clean)
}
times <- times_clean
metadata$times <- NULL
metadata[,"times"] <- as.numeric(times)

metadata[metadata$treatment == "NAIVE (UNTREATED)",3] <- c(1,1,1)
levels(metadata$treatment)[2] <- "CONTROL"
levels(metadata$treatment)[3] <- "COLTROL"

metadata <- metadata[with(metadata, order(treatment, times)),]
write.table(metadata, file = "metadata.tsv", sep = "\t", row.names = F)

# Manual editing of data!
metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t")
metadata = metadata[with(metadata, order(bleomycin, times)),]
write.table(metadata, file = "metadata.tsv", sep = "\t", row.names = F)
# Manual editing of data!
metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t")
rownames(metadata) <- metadata$sample_names
gsm <- substr(metadata$sample_names, start = 8, stop = 10)
Names <- data.frame(gsm = rownames(metadata), sampleID = paste0("T",metadata$times,"_C",
                                                     metadata$control, "_B",
                                                     metadata$bleomycin,"_", gsm))
write.table(Names, file = "names.tsv", sep = "\t", row.names = F)
system("sed -i 's/_[CB]0//g' names.tsv")
Names <- read.table("names.tsv", header = TRUE, sep = "\t")
rownames(metadata) <- Names$sampleID
########## DESIGN MATRIX ############
library(maSigPro)

edesign <- metadata[,3:6]


