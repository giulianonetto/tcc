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
png("Development_files/figure-docx/hclust_bauer2015.png", res = 100)
hcluster <- see.genes(sigs$sig.genes$BleomycinvsControl, show.fit = T, dis = design$dis,
                      cluster.method = "hclust", legend = F, newX11 = F)
dev.off()
png("Development_files/figure-docx/hclust_bauer2015_better.png")
see.genes(sigs$sig.genes$BleomycinvsControl, show.fit = T, dis = design$dis,
          cluster.method = "hclust", legend = F, newX11 = F)
dev.off()
# Plot Genes 
{
  genes.names <- c("Cd80", "Cd86", "Cd68", "Ccl2", "Ccl12", "Ccl17", "Ccl22", "Ccl22_1", 
                   "Ccl24", "Ccl7", "Ccl7_1", "Cxcl12", "Fcer1a", "Ms4a2_1",
                   "Tgfb3", "Il1b", "Il24", "Il15", "Il33_1", "Il13ra2", "Il13ra1_3", 
                   "Nostrin", "Mmp23", "Mmp7", "Mmp14", "Mmp19_1",  "Mmp12", 
                   "Mmp2_1", "Lbp", "Il1rl1_1", "C1qtnf9", "Mcpt2", "Mcpt8l3",
                   "Tgfbi_1", "Ccr2")
  genes.groups <- data.frame(Clusters = hcluster$cut[which(names(hcluster$cut) %in% genes.names)])
  genes <- genexp[which(rownames(genexp) %in% genes.names),]
  
  for (g in genes.names){ 
    main = str_glue("{g} - cluster {genes.groups[g,]}")
    png(str_glue("Development_files/figure-docx/masigpro/{g}.png"))
    PlotGroups(genexp["Arg",], edesign = design$edesign, show.fit = T,
               dis = design$dis, groups.vector = design$groups.vector, 
               main = "Apoe", legend = F)
    dev.off()
  }

}


# TGFb1
library(ggsignif); library(car); library(ggpubr)
genesLong <- genexp["Tgfb1",] %>% reshape2::melt()
genesLong$time <- factor(as.numeric(str_extract(genesLong$variable, "[0-9]+")), 
                         levels = c(0, 3, 7, 14, 21, 28))
genesLong$treatment <- str_extract(genesLong$variable, "[^T][A-Z]+")
genesLong$treatment[1:3] <- "_CTRL"
genesLong <- rbind(genesLong, genesLong[1:3,])
genesLong[(nrow(genesLong) - 2):nrow(genesLong), "treatment"] <- "_BLEO"

fit <- lm(value ~ poly(as.numeric(as.character(time)), 3), genesLong[-c(10,13),])
outlierTest(fit)
qqPlot(fit, id.method = "identity", simulate = T)
qqnorm(genesLong[-c(10,13),"value"])
qqline(genesLong[-c(10,13),"value"])  # fairly normal
genesLong$treatment <- factor(genesLong$treatment, 
                              labels = c("Bleomycin", "Control"))
wilcox_test(value ~ treatment, data = genesLong, distribution = "exact")
p1 <- ggplot(genesLong[-c(10,13),], aes(x = as.numeric(as.character(time)), y = value, 
                      group = treatment, color = treatment)) +
        geom_point() + geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
        xlab("Days after treatment") + ylab("Expression") +
        theme(legend.title = element_blank())
p2 <- ggplot(genesLong, aes(x = treatment, y = value, group = treatment)) +
        geom_point(aes(color = treatment)) +
        geom_boxplot(aes(color = treatment)) +
        geom_signif(comparisons = list(c("Bleomycin", "Control")),
                    map_signif_level = T, tip_length = 0, test = "t.test") +
        ylab("Expression") +
        theme(legend.title = element_blank(), axis.title.x = element_blank())
ggarrange(p1, p2, ncol = 2, common.legend = T, legend = "right") 

# IL1b1
genesLong <- genexp["Il1b",] %>% reshape2::melt()
genesLong$time <- as.numeric(str_extract(genesLong$variable, "[0-9]+"))
genesLong$treatment <- factor(str_extract(genesLong$variable, "[^T][A-Z]+"), 
                              labels = c("Bleomycin", "Control", "Untreated"))
genesLong$treatment[1:3] <- "Control"
genesLong <- rbind(genesLong, genesLong[1:3,])
genesLong[(nrow(genesLong) - 2):nrow(genesLong), "treatment"] <- "Bleomycin"

fit <- lm(value ~ poly(time, 3), genesLong)
outlierTest(fit)
qqPlot(fit, id.method = "identity", simulate = T) # 13 and 12 should leave

genesLong$rats <- paste(genesLong$treatment, genesLong$time, sep = ".")
genesLong$treatment <- factor(genesLong$treatment, levels = c("Bleomycin", "Control")) 
genesLong$time <- factor(genesLong$time)
# multiple comparisons para il1b1
pairwise.wilcox.test(genesLong[which(genesLong$time != 0),]$value, 
                     genesLong[which(genesLong$time != 0),]$treatment, 
                     p.adjust.method = "fdr")


# Check IL13ra1
a1 <- melt(genexp[c("Il13ra1", "Il13ra1_1", "Il13ra1_2", "Il13ra1_3"), ])
a1 <- aggregate(value ~ variable, a1, FUN=median)
a1$treatment <- ifelse(grepl("BLEO",a1$variable), "Bleomycin", "Control")
a1[(nrow(a1)+1):(nrow(a1)+3),] <- a1[1:3,]
a1[(nrow(a1)-2):(nrow(a1)),3] <- "Bleomycin"
a1$time <- str_extract(a1$variable, "[0-9]+") %>% as.numeric()
ggplot(a1, aes(x = time, y = value, color = treatment)) + 
  geom_point()+
  geom_smooth(method = "glm", formula = y ~ poly(x, 3))
