library(tidyverse)
library(plyr)
library(stringr)
library(TcGSA)
library(biomaRt)
library(dplyr)
library(Biobase)
library(limma)

# Load Data
setwd("/home/giulianonetto/windows/tcc/rmd-files")
eset <- readRDS("data/GSE48455/suppdata/ExpressionSetClean.rds")

# Expression Matrix

genexp <- exprs(eset)

# Genesets - Map orthologs --> PYTHON

system("#python ~/windows/tcc/rmd-files/Scripts/S2.3_mapXueorthologs.py")

# Load mapped genesets
genesets <- read.csv("data/XueModules/genesets.rat.gmt", sep = "\t",
                     header = F, row.names = 1, stringsAsFactors = F,
                     na.strings = c("", " ", "NA"))
M1_M2_modules <- rownames(genesets[which(genesets$V2 %in% c("M1", "M2")),])

# Get only M1- and M2-associated modules
genesets <- genesets %>% dplyr::filter(V2 == "M1" | V2 == "M2") %>% as.data.frame()
rownames(genesets) <- M1_M2_modules

# remove genes from Xue modules that match no genes in expression matrix
genesets_copy = genesets
for (i in 1:nrow(genesets)){
  nas <- is.na(genesets[i,2:ncol(genesets)])
  query <- genesets[i,2:ncol(genesets)][!nas]
  print(str_glue("Query length: {length(query)}"))
  matches <- query %in% rownames(genexp)
  print(str_glue("Remove: {genesets[i,c(FALSE, !matches)]}"))
  genesets[i,c(FALSE, !matches)] <- NA
}
for (i in 1:nrow(genesets)){
  nas <- is.na(genesets[i,])
  genesets[i,] <- c(genesets[i,][!nas], genesets[i,][nas])
}
genesets$modules <- rownames(genesets)
x <- "modules"
genesets <- genesets[c(x, setdiff(names(genesets), x))]
for (i in 1:nrow(genesets)){
  line <- t(data.frame(genes = genesets[i,][!is.na(genesets[i,])]))
  write.table(line, "data/XueModules/genesets_subset.gmt", quote = F, 
              sep = "\t", col.names = F, row.names = F, append = T)
}


# Load new gmt file
Xue.gmt <- GSA::GSA.read.gmt("data/XueModules/genesets_subset.gmt")

