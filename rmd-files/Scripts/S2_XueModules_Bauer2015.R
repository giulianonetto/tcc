library(GEOquery)
library(stringr)
library(limma)
library(vsn)
library(rafalib)
library(genefilter)
library(arrayQualityMetrics)
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

# Calculate mean between probes that map to the same ge

dupl_probes <- duplicated(RG_vsn_cDNA$genes$ProbeName) # 2367 repeats
dupl_probes <- RG_vsn_cDNA$genes[dupl_probes, "ProbeName"]
M <- data.frame(RG_vsn_cDNA$M)
A <- data.frame(RG_vsn_cDNA$A)
# M$ProbeNames <-RG_vsn_cDNA$genes$ProbeName
for (probe in dupl_probes){
  entries <- RG_vsn_cDNA[which(RG_vsn_cDNA$genes$ProbeName == probe),]
  print(str_glue("{probe}:  {dim(entries)[1]} repeats."))
  M
  
}

rownames(RG_vsn_cDNA$genes) <- RG_vsn_cDNA$genes$ProbeName 
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



