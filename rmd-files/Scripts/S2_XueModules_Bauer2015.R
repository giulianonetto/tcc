library(GEOquery)
library(dplyr)
library(stringr)
library(limma)
# library(vsn)
library(rafalib)
# library(genefilter)
library(arrayQualityMetrics)
# -----> always make sure you know where you stand < ------- #
setwd("/home/giulianonetto/windows/tcc/rmd-files")
# these guys take a while (over 1GB of space for Supp Data)
geo <- getGEO("GSE48455")
# geoSup <- getGEOSuppFiles("GSE48455", makeDirectory = F, baseDir = "~/windows/tcc/rmd-files/data/GSE48455/suppdata/")

# Make .gmt file
xue_path <- "/home/giulianonetto/windows/tcc/checklist/genesigs/papers/Xue, 2014/modules.csv"
xue <- read.csv(xue_path, sep = "\t", stringsAsFactors = F)
head(xue)
gmtfile <- data.frame(cbind(names(xue),"Other",t(as.matrix(xue))[,2:884]), stringsAsFactors = F)
gmtfile[7:9, 2] <- "M1"; gmtfile[13:15,2] <- "M2"; gmtfile[c(30, 32, 33),2] <- "TPP"
write.table(gmtfile, "data/XueModules/genesets.gmt", col.names = F, row.names = F, quote = F)

# Make gene exp matrix suitable for TcGSA - data from Bauer, 2015 

## Load RAW data and annotate with gene names
targets <- readTargets(path = "data/GSE48455/suppdata/", row.names = "Name")
RG <- read.maimages(targets$FileName, source="agilent")
spottypes <- readSpotTypes(file="data/GSE48455/suppdata/SpotTypes.txt")
RG$genes$Status <- controlStatus(spottypes, RG)


# VSN normalization
RG_vsn <- normalizeVSN(RG)
meanSdPlot(RG_vsn$M, ranks = TRUE) # Not so cool...
plotDensities(RG_vsn$M, main = "RG_vsn", legend = F) # Awesome!

# Remove control probes and rename columns

colnames(RG_vsn) <- str_glue('{str_extract(colnames(RG_vsn$M), "GSM[0-9]+")} / {targets$Name}')
probes <- RG_vsn$genes$Status == "cDNA"
RG_vsn_cDNA <- RG_vsn[probes,]

# Calculate median between probes that map to the same gene
# ----> BE AWARE THIS IS NOT VERY COOL - SEE PLOTS IN ()

dupl_probes <- duplicated(RG_vsn_cDNA$genes$ProbeName) 
dupl_probes_count <- sum(dupl_probes) # 2367 repeats
dupl_probes <- unique(RG_vsn_cDNA$genes[dupl_probes, "ProbeName"])
dupl_genes <- duplicated(RG_vsn_cDNA$genes$GeneName)
dupl_genes <- unique(RG_vsn_cDNA$genes[dupl_genes, "GeneName"])

# See how messy the distribution is for only one gene/probe on 3 different arrays
plotDensities(RG_vsn_cDNA[which(RG_vsn_cDNA$genes$ProbeName == dupl_probes[1]),1:3]$M)
plotDensities(RG_vsn_cDNA[which(RG_vsn_cDNA$genes$GeneName == dupl_genes[1]),1:3]$M)

M <- data.frame(RG_vsn_cDNA$M)
A <- data.frame(RG_vsn_cDNA$A)

# From a list of duplicated probe names, take their expression median and
# return new M and A data.frames with updated values 
# i.e. all entries mapping to same probe will have same values
for (probe in dupl_probes){
  entries <- RG_vsn_cDNA[which(RG_vsn_cDNA$genes$ProbeName == probe),]
  print(str_glue("{probe}:  {dim(entries)[1]} repeats."))
  entries$M <- 
  M[which(RG_vsn_cDNA$genes$ProbeName == probe),] <- data.frame(entries$M) %>%
    summarise_all(median) %>% as.vector()
  A[which(RG_vsn_cDNA$genes$ProbeName == probe),] <- data.frame(entries$A) %>%
    summarise_all(median) %>% as.vector()  
}

RG_vsn_cDNA_dereplicated <- RG_vsn_cDNA
RG_vsn_cDNA_dereplicated$M <- as.matrix(M)
RG_vsn_cDNA_dereplicated$A <- as.matrix(A)
dimnames(RG_vsn_cDNA_dereplicated$M)[[2]] <- dimnames(RG_vsn_cDNA$M)[[2]]
# Make sure no other alterations were made:
all.equal(RG_vsn_cDNA_dereplicated$M[!which(RG_vsn_cDNA_dereplicated$genes$ProbeName %in% 
                                      dupl_probes),], 
          RG_vsn_cDNA$M[!which(RG_vsn_cDNA$genes$ProbeName %in% 
                                 dupl_probes),])

# Make sure the difference between 
# the copy and the unique(copy) matches the number of repeats:
dim(RG_vsn_cDNA_dereplicated$M)[1] - dim(unique(RG_vsn_cDNA_dereplicated$M))[1] == dupl_probes_count # Awesome!
RG_vsn_cDNA_dereplicated$M <- unique(RG_vsn_cDNA_dereplicated$M)
RG_vsn_cDNA_dereplicated$A <- unique(RG_vsn_cDNA_dereplicated$A)

saveRDS(RG_vsn_cDNA_dereplicated, file = "data/GSE48455/suppdata/norm_vsn/RG_vsn_cDNA_dereplicated.rds")
# STOPPED HERE!!!






rownames(RG_vsn_cDNA_dereplicated$genes) <- RG_vsn_cDNA$genes$ProbeName

meanSdPlot(RG_vsn_cDNA$M, ranks = T)
plotDensities(RG_vsn_cDNA$M, legend = F, 
              main = "Distribution of all 72 arrays after VSN Norm.")


## Get Gene annotations
features <- getGEO(filename = "data/GSE48455/GSE48455_series_matrix.txt.gz")@featureData@data
features <- features[,c("REFSEQ", "GB_ACC", "GENE_SYMBOL", "GENE_NAME")]

RG_vsn_cDNA$genes <- features # RG_norm_cDNA$genes <- merge(RG_norm_cDNA$genes, features, by="row.names",all.x=TRUE)
RG_vsn_cDNA$targets <- targets




# arrayQualityMetrics

aQM <- arrayQualityMetrics(RG_vsn_cDNA, outdir = "data/GSE48455/suppdata/arrayQM",
                           intgroup = targets$Cy3,reporttitle = "aQM-RG_vsn_cDNA")



