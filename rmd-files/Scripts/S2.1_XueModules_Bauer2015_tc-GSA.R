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

# # Genesets - Map orthologs --> PYTHON
# 
# system("#python ~/windows/tcc/rmd-files/Scripts/S2.3_mapXueorthologs.py")
# 
# # Load mapped genesets
# genesets <- read.csv("data/XueModules/genesets.rat.gmt", sep = "\t",
#                      header = F, row.names = 1, stringsAsFactors = F,
#                      na.strings = c("", " ", "NA"))
# M1_M2_modules <- rownames(genesets[which(genesets$V2 %in% c("M1", "M2")),])
# 
# # Get only M1- and M2-associated modules
# genesets <- genesets %>% dplyr::filter(V2 == "M1" | V2 == "M2") %>% as.data.frame()
# rownames(genesets) <- M1_M2_modules
# 
# # remove genes from Xue modules that match no genes in expression matrix
# genesets_copy = genesets
# for (i in 1:nrow(genesets)){
#   nas <- is.na(genesets[i,2:ncol(genesets)])
#   query <- genesets[i,2:ncol(genesets)][!nas]
#   print(str_glue("Query length: {length(query)}"))
#   matches <- query %in% rownames(genexp)
#   print(str_glue("Remove: {genesets[i,c(FALSE, !matches)]}"))
#   genesets[i,c(FALSE, !matches)] <- NA
# }
# for (i in 1:nrow(genesets)){
#   nas <- is.na(genesets[i,])
#   genesets[i,] <- c(genesets[i,][!nas], genesets[i,][nas])
# }
# genesets$modules <- rownames(genesets)
# x <- "modules"
# genesets <- genesets[c(x, setdiff(names(genesets), x))]
# for (i in 1:nrow(genesets)){
#   line <- t(data.frame(genes = genesets[i,][!is.na(genesets[i,])]))
#   write.table(line, "data/XueModules/genesets_subset.gmt", quote = F, 
#               sep = "\t", col.names = F, row.names = F, append = T)
# }
# 

# Load new gmt file
Xue.gmt <- GSA::GSA.read.gmt("data/XueModules/genesets_subset.gmt")
design <- pData(eset)[,-c(3,5)]
design$Cy3 <- factor(design$Cy3)
design$Name <- design$Name %>% as.factor()

tcgsa <- TcGSA.LR(expr = genexp, gmt = Xue.gmt,
                  design = design,
                  subject_name = "Name", 
                  time_name = "Times",
                  group_name = "Cy3",
                  time_func = "cubic")
summary(tcgsa)
sigs <- signifLRT.TcGSA(tcgsa)
multi_test <- multtest.TcGSA(tcgsa)
plot1GS(tcgsa$Estimations, gmt = Xue.gmt, TimePoint = "Times",
        Subject_ID = "Name", clustering = F, 
        geneset.name = "Module15",
        gg.add=list(ggplot2::theme(legend.position="none")))

clust <- clustTrend(tcgsa, genexp, design$Name, 
                    design$Times, baseline = 0,
                    group_of_interest = "bleomycin")
plot(x = tcgsa, expr = genexp, Subject_ID = design$Name, 
     TimePoint = design$Times,
     group_of_interest = "bleomycin", 
     clust_trends = clust,
     subtitle="Bleomycin vs Ctrl", 
     cex.label.row=1, cex.label.col=1, cex.main=0.7,
     heatmap.width=0.2, dendrogram.size=0.1, margins=c(2,6),
     heatKey.size=0.8)

## VER https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html