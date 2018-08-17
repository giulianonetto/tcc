########### DOWNLOADS ############
library(dplyr)
library(GEOquery)
mydir <- getwd()
datasets <- read.csv("datasets.tsv", sep = "\t", header = T, stringsAsFactors = F)
geo <- (datasets %>% filter(author == "Bauer", year == 2015) %>% select(geo) %>% unlist())[2]
system("mkdir bauer_2015_GSE48455")
analysis_dir <- "bauer_2015_GSE48455/"
setwd(analysis_dir)

#filePaths <- getGEOSuppFiles(GEO = geo)
data <- getGEO(geo)[[1]]


########### METADATA ############
library(stringr)
metadata <- pData(data)
sample_names <- metadata$geo_accession
treatment <- metadata$`treatment:ch1`
times <- metadata$`time:ch1`
data_processing <- metadata$data_processing[1]
organism <- metadata$organism_ch1[1]
scan <- metadata$scan_protocol[1]
times <- metadata$`time:ch1`
metadata_clean <- data.frame(sample_names, treatment, times)
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
metadata_clean$times <- NULL
metadata_clean[,"times"] <- times
########## DESIGN MATRIX ############
library(maSigPro)
levels(metadata_clean$treatment)[2] <- "NAIVE"
metadata <- metadata_clean[with(metadata_clean, order(times, treatment)),]
metadata[6:8,3] <- rep(1, 3)
metadata <- metadata[with(metadata, order(times, treatment)),]
rownames(metadata) <- c()
rownames(metadata)
write.table(metadata, file = "metadata.tsv", sep = "\t", row.names = F)
