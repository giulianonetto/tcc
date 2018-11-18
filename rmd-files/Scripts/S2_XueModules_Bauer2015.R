library(GEOquery)
library(tidyverse)
library(stringr)
library(limma)
library(vsn)
library(rafalib)
library(arrayQualityMetrics)
library(ggfortify)
library(ggpubr)
library(FactoMineR)
library(convert)

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
write.table(gmtfile, "data/XueModules/genesets.gmt", 
            col.names = F, row.names = F, quote = F, sep = "\t")

# Make gene exp matrix suitable for TcGSA - data from Bauer, 2015 

## Load RAW data and annotate with gene names
targets <- readTargets(path = "data/GSE48455/suppdata/", row.names = "Name")
RG <- read.maimages(targets$FileName, source="agilent")
spottypes <- readSpotTypes(file="data/GSE48455/suppdata/SpotTypes.txt")
RG$genes$Status <- controlStatus(spottypes, RG)
RG$genes$ProbeName <- make.names(RG$genes$ProbeName, unique = TRUE)
RG$genes$GeneName <- make.names(RG$genes$GeneName, unique = TRUE)
RG$targets <- targets

# VSN normalization
RG_vsn <- normalizeVSN(RG)
rownames(RG_vsn$M) <- RG_vsn$genes$ProbeName
meanSdPlot(RG_vsn$M, ranks = TRUE) # Not so cool...
plotDensities(RG_vsn$M, main = "RG_vsn", legend = F) # Awesome!

# Remove control probes, rename columns and
# make sure probes and genes have unique names

colnames(RG_vsn) <- str_glue('{str_extract(colnames(RG_vsn$M), "GSM[0-9]+")} / {targets$Name}')
probes <- RG_vsn$genes$Status == "cDNA"
RG_vsn_cDNA <- RG_vsn[probes,]


saveRDS(RG_vsn_cDNA, file = "data/GSE48455/suppdata/norm_vsn/RG_vsn_cDNA.rds")


## Get Gene annotations
geo <- getGEO(filename = "data/GSE48455/GSE48455_series_matrix.txt.gz")
features <- geo@featureData@data
features <- features[,c("REFSEQ", "GB_ACC", "GENE_SYMBOL", "GENE_NAME")]
features$GENE_SYMBOL <- make.names(features$GENE_SYMBOL, unique = TRUE)

M <- data.frame(RG_vsn_cDNA$M, stringsAsFactors = F)
write.table(M, "data/GSE48455/Buaer2015_vsn_cDNA_all.tsv", sep = "\t", quote = F)
write.table(features, "data/GSE48455/Bauer2015_features.tsv", sep = "\t",quote = F)
## Python processing...

M <- read.table("data/GSE48455/Bauer2015_vsn_cDNA_all_mapped.tsv", sep = "\t",
                 header = T)
M$Genes <- make.names(M$Genes, unique=T)

rownames(M) <- M$Genes
M$Genes <- NULL
RG_vsn_cDNA$M <- M

# Make ExpressionSet object
eset <- as(RG_vsn_cDNA, "ExpressionSet")
names <- sapply(pData(eset)$Name, str_split, pattern = "-")
pData(eset)$Times <- as.numeric(str_extract(unlist(names)[c(TRUE, FALSE)], "\\d+"))
pData(eset)$Phase <- cut(pData(eset)$Times, breaks = c(-1,0,7,15,29,56), include.lowest = F,
                         labels = c("Healthy", "Injury", "Early Fibrosis", "Late Fibrosis", "Healing"))
genexp <- exprs(eset)
rownames(genexp) <- toupper(rownames(genexp))
pheno <- new("AnnotatedDataFrame", data = pData(eset))
experiment <- new("MIAME", name="Bauer2015-normalized",
                  lab="Imuno Clinica, UFRGS", 
                  contact = "cruz_gnf@gmail.com",
                  title = "IPF time-series - rats", 
                  url = "github.com/giulianonetto/tcc")
eset <- new("ExpressionSet", exprs = genexp, phenoData = pheno,
            experimentData = experiment, annotation = "GPL7294")

saveRDS(eset, "data/GSE48455/ExpressionSet.rds")

eset <- readRDS("data/GSE48455/ExpressionSet.rds")


# Array Quality Metrics outlier detection

pheno <- pData(eset)
healthy <- pheno$Phase == "Healthy"
injury <-  pheno$Phase == "Injury"
fibrosis <- pheno$Phase == "Early Fibrosis" | pheno$Phase == "Late Fibrosis"
healing <- pheno$Phase == "Healing"
control <- pheno$Cy3 == "control"
bleomycin <- pheno$Cy3 == "bleomycin"

aQM.all <- arrayQualityMetrics(eset, outdir = "data/GSE48455/suppdata/arrayQM.all",
                            intgroup = c("Cy3", "Phase"),
                            reporttitle = "Bauer2015 - all", force = TRUE)

aQM.healty <- arrayQualityMetrics(eset[,healthy], outdir = "data/GSE48455/suppdata/arrayQM.healthy",
                           intgroup = c("Cy3", "Times"), reporttitle = "Bauer2015 - healthy rats")

aQM.injury <- arrayQualityMetrics(eset[,injury], outdir = "data/GSE48455/suppdata/arrayQM.injury",
                                  intgroup = c("Cy3", "Times"),
                                  reporttitle = "Bauer2015 - injury phase")
aQM.fibrosis <- arrayQualityMetrics(eset[,fibrosis], outdir = "data/GSE48455/suppdata/arrayQM.fibrosis",
                                    intgroup = c("Cy3", "Times", "Phase"),
                                    reporttitle = "Bauer2015 - fibrosis phase")
aQM.healing <- arrayQualityMetrics(eset[,healing], outdir = "data/GSE48455/suppdata/arrayQM.healing",
                                            intgroup = c("Cy3", "Times"),
                                            reporttitle = "Bauer2015 - healing phase")
aQM.control <- arrayQualityMetrics(eset[,control], outdir = "data/GSE48455/suppdata/arrayQM.control",
                                   intgroup = c("Cy3", "Times", "Phase"), reporttitle = "Bauer2015 - controls")
aQM.bleomycin <- arrayQualityMetrics(eset[,bleomycin], outdir = "data/GSE48455/suppdata/arrayQM.bleomycin",
                                     intgroup = c("Cy3", "Times", "Phase"), reporttitle = "Bauer2015 - bleomycin")


Outliers <- rownames(pheno) %in% c("GSM1179148.B3.5", "GSM1179104.V42.4", 
                                   "GSM1179124.B28.3", "GSM1179133.B7.5", "GSM1179085.B3.2")
PCA.all <- eset@assayData$exprs %>% t() %>% prcomp()
PCA.clean <- eset[,!Outliers]@assayData$exprs %>% t() %>% prcomp()
# ggplot stuff for pca
{ colours <- scale_color_manual(values = c("deeppink", "red", "black", "blue", "green"))
  colours_Cy3 <- scale_color_manual(values = c("black", "red"))
  scales_y <- scale_y_continuous(limits=c(-0.6, 0.60))
  margins <- theme(plot.margin = unit(c(1.5, 1, 0.5, 1), "cm"))}
# PCA by Disease Phase
pca.all.1.2 <- autoplot(PCA.all, x = 1, y = 2, 
                        data = pheno, colour = "Phase", 
                        shape = "Cy3", size = 3)+colours+scales_y+margins
pca.all.1.3 <- autoplot(PCA.all, x = 1, y = 3, 
                        data = pheno, colour = "Phase", 
                        shape = "Cy3", size = 3)+colours+scales_y+margins
pca.clean.1.2 <- autoplot(PCA.clean, x = 1, y = 2,
                        data = pheno[!Outliers,], colour = "Phase", 
                        shape = "Cy3", size = 3)+colours+scales_y+margins
pca.clean.1.3 <- autoplot(PCA.clean, x = 1, y = 3,
                        data = pheno[!Outliers,], colour = "Phase", 
                        shape = "Cy3", size = 3)+colours+scales_y+margins
gg_phase <- ggarrange(pca.all.1.2, pca.all.1.3, 
                      pca.clean.1.2, pca.clean.1.3, nrow = 2, ncol = 2,
                      common.legend = TRUE, legend = "bottom")
annotate_figure(gg_phase, 
                top = text_grob("Before outlier removal",
                                face = "bold", size = 20, 
                                vjust = 1),
                bottom = text_grob("After outlier removal",
                                   face = "bold", size = 20,
                                   vjust = -17.5))
ggsave("data/GSE48455/suppdata/pca.phase.jpg", 
       width = 20, height = 20, units = "cm")


# PCA by treatment condition
pca.all.1.2_Cy3 <- autoplot(PCA.all, x = 1, y = 2, 
                            data = pheno, colour = "Cy3", 
                            size = 3)+colours_Cy3+scales_y+margins
pca.all.1.3_Cy3 <- autoplot(PCA.all, x = 1, y = 3, 
                            data = pheno, colour = "Cy3", 
                            size = 3)+colours_Cy3+scales_y+margins
pca.clean.1.2_Cy3 <- autoplot(PCA.clean, x = 1, y = 2, 
                              data = pheno[!Outliers,], colour = "Cy3", 
                              size = 3)+colours_Cy3+scales_y+margins
pca.clean.1.3_Cy3 <- autoplot(PCA.clean, x = 1, y = 3, 
                              data = pheno[!Outliers,], colour = "Cy3", 
                              size = 3)+colours_Cy3+scales_y+margins

gg_Cy3 <- ggarrange(pca.all.1.2_Cy3, pca.all.1.3_Cy3, 
                    pca.clean.1.2_Cy3, pca.clean.1.3_Cy3, nrow = 2, ncol = 2,
                    common.legend = TRUE, legend = "bottom")

annotate_figure(gg_Cy3, 
                top = text_grob("Before outlier removal",
                                face = "bold", size = 20, 
                                vjust = 1),
                bottom = text_grob("After outlier removal",
                                   face = "bold", size = 20,
                                   vjust = -17.5))

ggsave("data/GSE48455/suppdata/pca.disease.jpg", 
       width = 20, height = 20, units = "cm")


# Quantile-Quantile plots

longEset <- as.data.frame(exprs(eset)) %>% gather()
ggplot(longEset)+
  geom_qq(aes(sample = value, color = key))+
  theme(legend.position="none")+ggtitle("Q-Q plot with outliers")

ggsave("data/GSE48455/suppdata/qq_before.png", 
       width = 20, height = 20, units = "cm")

longEset.clean <- as.data.frame(exprs(eset[,!Outliers])) %>% gather()
ggplot(longEset.clean)+
  geom_qq(aes(sample = value, color = key))+
  theme(legend.position="none")+ggtitle("Q-Q plot without outliers")
ggsave("data/GSE48455/suppdata/qq_after.png", 
       width = 20, height = 20, units = "cm")

saveRDS(eset[,!Outliers], "data/GSE48455/suppdata/ExpressionSetClean.rds")


# More PCAs colored by phase...

## Eset without outliers

eset <- readRDS("data/GSE48455/suppdata/ExpressionSetClean.rds")
pheno <- pData(eset); expt <- t(exprs(eset)); rm(eset)
{ pca <- prcomp(expt, scale. = T, center = T)
  pca$sdev %>% plot(type="b")
  pca %>% summary() %>% print()
  barplot(abs(colSums(pca$x)))
  autoplot(pca, x = 1, y = 2, 
           data = pheno, colour = "Phase", shape = "Cy3",
           size = 3)+colours+scales_y+margins}

# rever gene expression matrix!!!
