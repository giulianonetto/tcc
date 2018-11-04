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
M1_M2_modules <- rownames(genesets[which(genesets$V2 %in% c("M1", "M2")),])

# Get only M1- and M2-associated modules
genesets <- genesets %>% dplyr::filter(V2 == "M1" | V2 == "M2") %>% as.data.frame()
rownames(genesets) <- M1_M2_modules

# Find genes in modules that have no match in gene expression matrix

not_found = vector("list", length = 6)
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
not_found_df <- ldply(not_found, rbind)
rownames(not_found_df) <- M1_M2_modules
write.table(not_found_df, "data/GSE48455/suppdata/not_found.gmt",
            row.names = T, col.names = F, sep = '\t', quote = F)

# downloading and sorting gene_ingo from ncbi
system("#python /home/giulianonetto/windows/tcc/storage/gene_info/get_aliases.py")

# get aliases for not found genes 
## output: /home/giulianonetto/windows/tcc/rmd-files/data/GSE48455/suppdata/not_found_map.gmt
system("#/home/giulianonetto/windows/tcc/rmd-files/Scripts/S2.4_not_found_Xue_modules.py")

# map aliases and Xue modules genes - delete genes not mapped

aliases <- read.csv("/home/giulianonetto/windows/tcc/rmd-files/data/GSE48455/suppdata/not_found_map.gmt",
                    sep = "\t", header = F, stringsAsFactors = F)
# should be something like below
not_found <- not_found %>% unlist()
counter = 0
for (gene in not_found){
  query_alias = grepl(str_glue("^{gene} "), aliases$V1)
  alias_list <- str_split(aliases[query_alias,1], " ") %>% unlist()
  for (alias.match in alias_list){
    if (alias.match %in% rownames(genexp)){
      counter = counter + 1
      pos <- which(genesets == gene, arr.ind = TRUE)
      print(pos)
      genesets[pos] <- alias.match
      not_found_df
    }
  }
}
coutner # got more 70 genes

# remove not_found genes_ref
for ()