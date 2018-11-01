library(tidyverse)
library(plyr)
library(stringr)
library(TcGSA)
library(biomaRt)
library(dplyr)
library(limma)
# Load Data
setwd("/home/giulianonetto/windows/tcc/rmd-files")
eset <- readRDS("data/GSE48455/suppdata/ExpressionSetClean.rds")
Xue.gmt <- GSA::GSA.read.gmt("data/XueModules/genesets.gmt")


# Build grouping indices

pheno <- pData(eset)
pheno$Rows <- rownames(pheno)
healthy <- pheno$Phase == "Healthy"
injury <-  pheno$Phase == "Injury"
fibrosis <- pheno$Phase == "Early Fibrosis" | pheno$Phase == "Late Fibrosis"
healing <- pheno$Phase == "Healing"
control <- pheno$Cy3 == "control"
bleomycin <- pheno$Cy3 == "bleomycin"



# Expression Matrix

genexp <- exprs(eset)

# Genesets - Map orthologs --> PYTHON

system("python ~/windows/tcc/rmd-files/Scripts/S2.3_mapXueorthologs.py")

# Load mapped genesets
genesets <- read.csv("data/XueModules/genesets.rat.gmt", sep = "\t",
                     header = F, row.names = 1, stringsAsFactors = F,
                     na.strings = c("", " ", "NA"))

# Find genes in gmt file that have no match in gene expression matrix

not_found = vector("list", length = 49)
for (i in 1:nrow(genesets)){
  nas <- is.na(genesets[i,2:ncol(genesets)])
  query <- genesets[i,2:ncol(genesets)][!nas]
  print(str_glue("Query length: {length(query)}"))
  matches <- query %in% rownames(genexp)
  unmatched <- query[!matches]
  perc = round(length(unmatched)/length(query) * 100)
  print(str_glue("Unmatched: {length(unmatched)} ({perc}%)"))
  if (perc < 100){
    not_found[[i]] <- unmatched
  }
}
not_found <- ldply(not_found, rbind)
rownames(not_found) <- rownames(genesets)
# Use biomaRt!
convert_human_gene_list <- function(x){
  
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  synonims = 
  no_matches = setdiff(x, mouse_matches[,1]) 
  
  x = data.frame(MGI.symbol = x)
  df <- merge(x, mouse_matches, by = "HGNC.symbol")
  
  return(list(df, data.frame(no_matches)))
}
library(org.Rn.eg.db)
alias2Symbol(not_found[3,3], species="Hs")

# downloading gene_ingo from ncbi
system("python /home/giulianonetto/windows/tcc/storage/gene_info/get_aliases.py")


