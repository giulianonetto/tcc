library(stringr)
library(dplyr)
library(ggplot2)
library(maSigPro)

setwd("/home/giulianonetto/windows/tcc/rmd-files")
filepath <- "/home/giulianonetto/windows/tcc/storage/bauer_2015_GSE48455/"

# Expression matrix and metadata - exclude time 56 and 42

metadata = read.csv(str_glue("{filepath}metadata_renamed.tsv"), sep = "\t", header = T)
metadata[1:3,3] = c(0,0,0) # untreated are taken as starting point
rownames(metadata) <- metadata$Geo
edes <- metadata[1:52,3:6]
edes <- edes[,c(1,4,2,3)]

genexp <- read.table(str_glue("{filepath}gene_expression_matrix_mapped"), sep = "\t", header = T, row.names = 1)

# Remove time points 42 and 56

bigtimes <- grepl("T56|T42", names(genexp))
genexp <- genexp[,!bigtimes]
genexp = genexp[, match(rownames(edes), colnames(genexp))]

# make design matrix
design <- make.design.matrix(edesign = edes, degree = 3) # cubic regression model

{
  ##### REGRESSION MODEL #####
  
  fit <- p.vector(genexp, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
  
  ##### FIND SIGNIFICANT GENES #####
  
  # step-wise regression
  tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
  
  # check goodness of fit - using R-squared
  (b <- mean(tstep$sol[,"R-squared"]>0.5)*100) # 31.7% for degree=3

  potential_outliers <- tstep$influ.info %>% colnames()  # check out potential outliers

}

# RUN AGAIN WITHOUT OUTLIERS

{
  genexp = genexp[!(rownames(genexp) %in% potential_outliers), ]
  
  ##### REGRESSION MODEL #####
  
  fit <- p.vector(genexp, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
  
  ##### FIND SIGNIFICANT GENES #####
  
  # step-wise regression
  tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
  
  # check goodness of fit - using R-squared
  (b <- mean(tstep$sol[,"R-squared"] > 0.5)*100) # 31.7% for degree=3
  
  (potential_outliers <- tstep$influ.info %>% colnames())  # check out potential outliers
  
  # No more outliers! (at least according to maSigPro base algorithm...)
  saveRDS(tstep, "data/GSE48455/suppdata/tstep_bauer2015.rds")
}

# Get and Plot Significant genes

sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")

clusters <- see.genes(sigs$sig.genes$BleomycinvsControl, show.fit = T, dis = design$dis,
                      cluster.method = "Mclust", k.mclust = TRUE)
clusters.df <- clusters$cut %>% as.data.frame()
clusters.df$genes <- rownames(clusters.df); rownames(clusters.df) <- NULL

hcluster <- see.genes(sigs$sig.genes$BleomycinvsControl, show.fit = T, dis = design$dis,
                      cluster.method = "hclust", legend = F)

# Plot Genes 
{
  genes.names <- c("Cd80", "Cd86", "Cd68", "Ccl2", "Ccl12", "Ccl17", "Ccl22", "Ccl22_1", 
                   "Ccl24", "Ccl7", "Ccl7_1", "Cxcl12", "Fcer1a", "Ms4a2_1",
                   "Tgfb3", "Il1b", "Il24", "Il15", "Il33_1", "Il13ra2", "Il13ra1_3", 
                   "Nostrin", "Mmp23", "Mmp7", "Mmp14", "Mmp19_1",  "Mmp12", 
                   "Mmp2_1", "Lbp", "Il1rl1_1", "C1qtnf9", "Mcpt2", "Mcpt8l3",
                   "Tgfbi_1")
  genes.groups <- data.frame(Clusters = hcluster$cut[which(names(hcluster$cut) %in% genes.names)])
  genes <- genexp[which(rownames(genexp) %in% genes.names),]
  
  for (g in genes.names){ 
    main = str_glue("{g} - cluster {genes.groups[g,]}")
    png(str_glue("Development_files/figure-docx/masigpro/{g}.png"))
    PlotGroups(genes[g,], edesign = design$edesign, show.fit = T,
               dis = design$dis, groups.vector = design$groups.vector, 
               main = main, legend = F)
    dev.off()
  }

}
