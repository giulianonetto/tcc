library(tidyverse)
library(stringr)
library(TcGSA)
library(biomaRt)

# Load Data
setwd("/home/giulianonetto/windows/tcc/rmd-files")
eset <- readRDS("data/GSE48455/suppdata/ExpressionSetClean.rds")
Xue.gmt <- GSA::GSA.read.gmt("data/XueModules/genesets.gmt")
genesets <- read.csv("data/XueModules/genesets.gmt", 
                     sep = "\t", header = F, row.names = 1)

is.na(genesets) <- genesets == ""

# Build grouping indeces
'
pheno <- pData(eset)
healthy <- pheno$Phase == "Healthy"
injury <-  pheno$Phase == "Injury"
fibrosis <- pheno$Phase == "Early Fibrosis" | pheno$Phase == "Late Fibrosis"
healing <- pheno$Phase == "Healing"
control <- pheno$Cy3 == "control"
bleomycin <- pheno$Cy3 == "bleomycin"
'


# Expression Matrix

genexp <- exprs(eset)
not_found = vector("list", length = 49)
for (i in 1:nrow(genesets)){
  nas <- is.na(genesets[i,2:ncol(genesets)])
  query <- genesets[i,2:ncol(genesets)][!nas]
  print(str_glue("Query length: {length(query)}"))
  matches <- query %in% rownames(genexp)
  unmatched <- query[!matches]
  perc = round(length(unmatched)/length(query) * 100)
  if (perc1 == perc2) {
  print(str_glue("Unmatched: {length(unmatched)} ({perc}%)"))
  }
  if (perc < 100){
    not_found[[i]] <- unmatched
  }
}
not_found <- ldply(not_found, rbind)

# BiomaRt Annotation fixes

