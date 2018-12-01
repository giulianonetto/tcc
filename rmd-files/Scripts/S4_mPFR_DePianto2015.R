######  Load packages and set working dir ######
{
  library(limma)
  library(convert)
  library(vsn)
  library(car)
  library(arrayQualityMetrics)
  library(ggpubr)
  library(ggfortify)
  library(ggrepel)
  library(stringr)
  library(dplyr)
  library(GEOquery)
  library(reshape2)
  setwd("/home/giulianonetto/windows/tcc/rmd-files")
  
}



###### Preprocessing - YOU CAN SKIP THIS ######


  # Takes a while!
{

#   geoSup <- getGEOSuppFiles("GSE53845", 
#                             makeDirectory = F, 
#                             baseDir = "~/windows/tcc/rmd-files/data/GSE53845/suppdata/")
#   write.table(pheno@data, "~/windows/tcc/rmd-files/data/GSE53845/suppdata/phenodata.tsv",
#               row.names = F, col.names = T, quote = F, sep = "\t")
}
  # Load raw data
{
  # bora usar GSE53845 (40 ipf e 8 controls)
  geo <- getGEO("GSE53845")
  pheno <- phenoData(geo$GSE53845_series_matrix.txt.gz)
  
  ## Load RAW data and annotate with gene names
  targets <- readTargets(path = "data/GSE53845/suppdata", row.names = "Name")
  RG <- read.maimages(targets$FileName, source = "agilent")
  spottypes <- readSpotTypes(file = "data/GSE53845/suppdata/SpotTypes.txt")
  RG$genes$Status <- controlStatus(spottypes, RG)
  RG$targets <- targets

    # Set up gene names
  genes <- geo$GSE53845_series_matrix.txt.gz@featureData@data
  get_gene_name <- function(probe, gene.df){
    if (probe %in% rownames(gene.df)) {
      gene.name <- gene.df[probe,"GENE_SYMBOL"] %>% str_to_title()
      if (gene.name == "") {
        gene.name <- "to.remove"
      }
    }
    else {
      gene.name <- "to.remove"
    }
    return(gene.name)
  }
  GeneNames <- sapply(RG$genes$ProbeName, get_gene_name, genes)
  RG$genes$GeneName <- make.names(unlist(GeneNames), unique = TRUE)
  RG$genes$ProbeName <- make.names(RG$genes$ProbeName, unique = TRUE)
}
  # VSN Normalization
{
  
  # VSN normalization
  RG_vsn <- normalizeVSN(RG)
  rownames(RG_vsn$M) <- RG_vsn$genes$GeneName
  meanSdPlot(RG_vsn$M, ranks = TRUE) # Pretty cool...
  plotDensities(RG_vsn$M, main = "RG_vsn", legend = F) # Awesome!
}
  # Remove unused probes
{
  # Remove control probes, rename columns and
  # make sure probes and genes have unique names
  condition <- str_extract_all(targets$Name, "control|IPF|biopsy|transplant|necropsy") %>% 
    lapply(paste, collapse = " ") %>% unlist()
  colnames(RG_vsn) <- str_glue('{str_extract(colnames(RG_vsn$M), "GSM[0-9]+")} / {condition}')
  valid.probes <- RG_vsn$genes$Status == "cDNA" & 
    !grepl("to.remove*", RG_vsn$genes$GeneName)
  RG_vsn_cDNA <- RG_vsn[valid.probes,]
  saveRDS(RG_vsn_cDNA, file = "data/GSE53845/suppdata/RG_vsn_cDNA.rds")
}
  # Build eSet objct
{
  # Build ExpressionSet object
  phenodf <- data.frame(Array = targets$SlideNumber,
                        Treatment = str_extract(condition, "IPF|control"),
                        Source = str_extract(condition, 
                                    "biopsy|transplant|necropsy"),
                        row.names = colnames(RG_vsn_cDNA))
  pheno <- new("AnnotatedDataFrame", 
               data = phenodf)
  experiment <- new("MIAME", name = "DePianto2015-normalized",
                    lab = "Imuno Clinica, UFRGS", 
                    contact = "cruz_gnf@gmail.com",
                    title = "IPF patients", 
                    url = "github.com/giulianonetto/tcc")
  eset <- new("ExpressionSet", exprs = RG_vsn_cDNA$M, phenoData = pheno,
              experimentData = experiment, annotation = "GPL6480")
  
  saveRDS(eset, "data/GSE53845/suppdata/ExpressionSet.rds")
  
}


###### LOAD PREPROCESSED DATA ######

{
  eset <- readRDS("data/GSE53845/suppdata/ExpressionSet.rds")
  pheno <- pData(eset)
}




###### QUALITY CONTROL - TAKES FOREVER!!! ######
  # ArrayQualityMetrics - takes forever!!
{
  # Array Quality Metrics outlier detection
  aQM.all <- arrayQualityMetrics(eset, outdir = "data/GSE53845/suppdata/arrayQM.all",
                                 intgroup = c("Treatment", "Source"),
                                 reporttitle = "DePianto2015 - all", force = TRUE)
}
  


###### EXPLORATORY DATA ANALYSIS ######
outliers.vec <- colnames(eset) %in% c("GSM1302034 / IPF transplant", 
                                      "GSM1302054 / control necropsy",
                                      "GSM1302065 / IPF biopsy")

# Compute PCAs - with outliers
{

  PCA.all <- eset@assayData$exprs %>% t() %>% prcomp(scale. = T, center = F)
  PCA.all.sum <- summary(PCA.all)$importance
  loadings.df <- data.frame(VAR = PCA.all.sum[2,1:10] * 100, 
                            PCs = 1:10, 
                            stringsAsFactors = F)
  pp1 <- autoplot(PCA.all, x = 1, y = 2, 
                  data = pheno, colour = "Source", 
                  shape = "Treatment", size = 3)
  pp2 <- ggplot(loadings.df,
                aes(y = VAR, 
                    x = PCs)) + 
    geom_point() + geom_line(color = "blue", group = 1) +
    xlab("Principal Components") +
    ylab("Variance Explained (%)") +
    scale_x_continuous(breaks = 1:10,
                     labels = paste0("PC", 1:10))
  plots <- ggarrange(pp1, pp2, ncol = 2)
  annotate_figure(plots, 
                  top = text_grob("PCA - IPF patients - with outliers\n",
                                  face = "bold", size = 20, 
                                  vjust = 1))
}

  # Compute PCAs - without outliers
{
  PCA.clean <- eset[,!outliers.vec]@assayData$exprs %>% t() %>% prcomp(scale. = T, center = F)
  PCA.clean.sum <- summary(PCA.clean)$importance
  loadings.df <- data.frame(VAR = PCA.clean.sum[2,1:10] * 100, 
                            PCs = 1:10, 
                            stringsAsFactors = F)
  pp1 <- autoplot(PCA.clean, x = 1, y = 2, 
                  data = pheno[!outliers.vec,],
                  colour = "Treatment", size = 3)
  pp2 <- autoplot(PCA.clean, x = 3, y = 2, 
                  data = pheno[!outliers.vec,],
                  colour = "Treatment", size = 3) 
  pp3 <- autoplot(PCA.clean, x = 1, y = 3, 
                  data = pheno[!outliers.vec,],  
                  colour = "Treatment", size = 3)
  
  pp4 <- ggplot(loadings.df,
                  aes(y = VAR, x = PCs)) + 
    geom_point() + geom_line(color = "blue", group = 1) +
    xlab("Principal Components") +
    ylab("Variance Explained (%)") +
    scale_x_continuous(breaks = 1:10,
                       labels = paste0("PC", 1:10))
  plots <- ggarrange(pp1, pp2, pp3, pp4, ncol = 2, nrow = 2, 
                     common.legend = T)
  annotate_figure(plots, 
                  top = text_grob("PCA - IPF and Control patients\n",
                                  face = "bold", size = 20, 
                                  vjust = 1))
}

  # Save ExpressionSetClean
{
  saveRDS(eset[,!outliers.vec], "data/GSE53845/suppdata/ExpressionSetClean.rds")
}




####### mPFR Calculations #######

# Load Clean Data
{
  eset <- readRDS("data/GSE53845/suppdata/ExpressionSetClean.rds")
  pheno <- pData(eset)
  genexp <- exprs(eset)
}

# Define functions

{
  
  # Calculate mPFR
  
  mPFR <- function(sets_data.frame, exprs_matrix, Source, Treat) {
    mpfr_data.frame <- data.frame()
    for (i in 1:nrow(sets_data.frame)) {
      set1.vec <- sets_data.frame[i,1] %>% unlist()
      set2.vec <- sets_data.frame[i,2] %>% unlist()
      set1 <- apply(exprs_matrix[set1.vec, ,drop = F], 2, mean) # means across genesets
      set2 <- apply(exprs_matrix[set2.vec, ,drop = F], 2, mean)
      set1_all <- exprs_matrix[set1.vec,] %>% unlist() %>% mean() # baseline means
      set2_all <- exprs_matrix[set2.vec,] %>% unlist() %>% mean()
      mpfr <- (set1 - set2) * (set2_all / set1_all)
      mpfr_data.frame[1:length(mpfr),i] <- mpfr
    }
    colnames(mpfr_data.frame) <- rownames(sets_data.frame)
    mpfr_data.frame$source <- Source
    mpfr_data.frame$treatment <- Treat
    return(mpfr_data.frame)
  }
  
  # Plot mPFR comparisons
  
  plot_signif <- function(df, formula, filters, 
                          titles = NULL, outTest = TRUE,
                          method = "wilcox.test", point = F) {
    
    df <- df[which(filters),]
    if (outTest == TRUE) {
      while (TRUE) {
        t <- outlierTest(lm(formula, data = df))
        if (sum(t$bonf.p < 0.05) < 1) {
          break
        }
        else {
          print(t)
          to.remove <- names(t$bonf.p) %>% as.numeric()
          df <- df[!(rownames(df) %in% to.remove),]
        }
      }
      print("Outliers check done!")
    }
    else print("Skip outlier tests...")
    p <- ggplot(df, aes(x = treatment, y = value, color = treatment))
    stat <- stat_compare_means(comparisons = list(c("IPF", "control")),
                               method = method,
                               aes(label = ..p.signif..))
    fac <- facet_wrap(~ variable)
    if (point == FALSE) {
      p + stat + fac + geom_boxplot()
    }
    else {
      p + stat + fac + geom_point() + geom_boxplot(alpha = 0.25)
    }
  }
  
  # Fit linear models to check 
  # correlation btw mPFR candidates and IPF conditions
  
  fit_mPFR <- function(df, Nvars = 28){
    fitdf <- data.frame(Rsquared = vector('numeric'),
                        pval = vector('numeric'),
                        Estimate = vector('numeric'), 
                        stringsAsFactors = F)
    for (i in 1:Nvars) {
      var.name <- colnames(df)[i]
      fit <- lm(df[,i] ~ treatment, data = df)
      s <- summary(fit)
      rs <- s$r.squared %>% as.numeric()
      p <- s$coefficients[2,4] %>% as.numeric()
      effect <- s$coefficients[2,1] %>% as.numeric()
      fitdf[i,1:3] <- c(rs, p, effect)
      rownames(fitdf)[i] <- var.name
      print(str_glue("{var.name}:\t{rs}\t{effect}\t{p}\n"))
    }
    return(fitdf)
  }
}

# Define clusters from Bauer2015

"

        note that ~ 70 - 80% of cluster's genes are NOT
        present in this dataset

"

{
  # define genesets from masigpro clusters

  hcluster <- readRDS("data/GSE48455/suppdata/hcluster_masigpro_bauer2015.rds")
  c1a = names(hcluster$cut)[which(hcluster$cut == 1)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c1 <- c1a[c1a %in% rownames(genexp)]; print(length(c1)/length(c1a) * 100)
  c1not <- c1a[!c1a %in% c1]
  c2a = names(hcluster$cut)[which(hcluster$cut == 2)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c2 <- c2a[c2a %in% rownames(genexp)]; print(length(c2)/length(c2a) * 100)
  c3a = names(hcluster$cut)[which(hcluster$cut == 3)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c3 <- c3a[c3a %in% rownames(genexp)]; print(length(c3)/length(c3a) * 100)
  c4a = names(hcluster$cut)[which(hcluster$cut == 4)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c4 <- c4a[c4a %in% rownames(genexp)]; print(length(c4)/length(c4a) * 100)
  c5a = names(hcluster$cut)[which(hcluster$cut == 5)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c5 <- c5a[c5a %in% rownames(genexp)]; print(length(c5)/length(c5a) * 100)
  c6a = names(hcluster$cut)[which(hcluster$cut == 6)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c6 <- c6a[c6a %in% rownames(genexp)]; print(length(c6)/length(c6a) * 100)
  c7a = names(hcluster$cut)[which(hcluster$cut == 7)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c7 <- c7a[c7a %in% rownames(genexp)]; print(length(c7)/length(c7a) * 100)
  c8a = names(hcluster$cut)[which(hcluster$cut == 8)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c8 <- c8a[c8a %in% rownames(genexp)]; print(length(c8)/length(c8a) * 100)
  c9a = names(hcluster$cut)[which(hcluster$cut == 9)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c9 <- c9a[c9a %in% rownames(genexp)]; print(length(c9)/length(c9a) * 100)

} 


# Define dataframe of comparisons

" check expression profiles for those
  in the S1-3 clusters (Il13ra1.3 is not present,
  but Il13ra1.0-2 are)
"

{
  # we take the mean expression for 3 probes that map to IL13ra1
  # REMOVED: Ccl12 (not in the dataset)
  genexp["Il13ra1",] <- colMeans(genexp[grepl("Il13ra1.*", rownames(genexp)),])
  
  s1l <- list("S1" = c("Il12b"), 
              "S2" = c("Ccr7", "Ccl22", "Ccl24"),
              "S3" = c("Ccl2", "Ccr2", "Il1b"), 
              "C1" = c1, "C2" = c2, "C3" = c4, "C4" = c5, 
              "C5" = c1, "C6" = c2, "C7" = c4, "C8" = c5, 
              "C9" = c1, "C10" = c2, "C11" = c4, "C12" = c5, 
              "C13" = c1, "C14" = c2, "C15" = c4, "C16" = c5,
              "C17" = c1, "C18" = c2, "C19" = c4, "C20" = c5,
              "C21" = c1, "C22" = c1, "C23" = c1,
              "C24" = c2, "C25" = c2
  )
  s2l <- list("S1" = c("Arg1"), 
              "S2" = c("Il13ra1", "Il13ra2", "Cd36", "Cxcr1", "St7"),
              "S3" = c("Apoe", "Mmp14", "Lbp", "Il1rl1"), 
              "C1" = c3, "C2" = c3, "C3" = c3, "C4" = c3,
              "C5" = c6, "C6" = c6, "C7" = c6, "C8" = c6,               
              "C9" = c7, "C10" = c7, "C11" = c7, "C12" = c7,
              "C13" = c8, "C14" = c8, "C15" = c8, "C16" = c8,
              "C17" = c9, "C18" = c9, "C19" = c9, "C20" = c9,
              "C21" = c2, "C22" = c4, "C23" = c5,
              "C24" = c4, "C25" = c5
  )
  df <- data.frame(Sets1 = I(s1l), Sets2 = I(s2l))
}


# Calculate mPFR

"
  Adapt functions for the current dataset metadata
"
{
  mpfr.all <- mPFR(df, genexp, pheno$Source, pheno$Treatment)
  mpfr.all.long <- melt(mpfr.all, id.vars = c("source", "treatment"))
}


# EDA

{
  # QQ plots
  ggplot(mpfr.all.long, aes(sample = value)) +
    facet_wrap(~ variable) +
    ggplot2::stat_qq() + ggplot2::stat_qq_line(color = "red")
  # Overall Density
  plotDensities(t(mpfr.all[,1:28])) # maybe you should normalize again
  ggplot(mpfr.all.long, aes(value)) +
    geom_density() + xlim(-10,10)
  ggplot(mpfr.all.long, aes(sample = value)) +
    ggplot2::stat_qq() + ggplot2::stat_qq_line(color = 'red')
  "
  NORMALIZE AGAIN???? - normalizeVSN, quantile normalization...
  library(preprocessCore)
  mpfr.norm <- normalize.quantiles(t(mpfr.all[,1:28]))
  qqnorm(mpfr.norm)
  ggplot(mpfr.all.long, aes(sample = value)) +
    ggplot2::stat_qq() + ggplot2::stat_qq_line(color = 'red')
  "
}


# PCA

{

  pca <- prcomp(mpfr.all[,1:28], scale. = T, center = F)
  PCA.sum <- summary(pca)$importance
  loadings.df <- data.frame(VAR = PCA.sum[2,1:10] * 100, 
                            PCs = 1:10, 
                            stringsAsFactors = F)
  
  ppp1 <- autoplot(pca, data = mpfr.all, colour = "treatment", 
                   size = 3)
  ppp2 <- autoplot(pca, x = 3, y = 2, data = mpfr.all, 
                   size = 3, colour = "treatment")
  ppp3 <- autoplot(pca, x = 1, y = 3, data = mpfr.all, 
                   size = 3, colour = "treatment")
  ppp4 <- ggplot(loadings.df, aes(y = VAR, x = PCs)) + 
                 geom_point() + geom_line(color = "blue", group = 1) +
                 xlab("Principal Components") +
                 ylab("Variance Explained (%)") +
                 scale_x_continuous(breaks = 1:10,
                                    labels = paste0("PC", 1:10))
  plots <- ggarrange(ppp1, ppp2, ppp3, ppp4, ncol = 2, nrow = 2, 
                     common.legend = T)
  annotate_figure(plots, 
                  top = text_grob("PCA - IPF and Control patients\n",
                                  face = "bold", size = 20, 
                                  vjust = 1))
}

# Boxplots
"
    SERIA LEGAL MOSTRAR O N EM CADA FACET!!!
    TAMBEM MODELAR CORRELACAO ENTRE VARIABLES DE HUMANOS E RATOS!!!!
"
{
  ivars <- 13; nvars <- 24
  filtro = mpfr.all.long$variable %in% unique(mpfr.all.long$variable)[ivars:nvars]
  form <- value ~ treatment
  plot_signif(mpfr.all.long, form,
              filtro, outTest = T,  method = "wilcox.test",
              point = T)
}

# Modelos

{
  dat <- fit_mPFR(mpfr.all) 
  dat$variables <- rownames(dat)
}

# volcano

{
ggplot(dat, aes(x = Estimate, y = -log10(pval), 
                label = variables, size = Rsquared)) +
  geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) + 
  labs(size = "R squared", color = NULL) +
  geom_text_repel(force = T) +
  xlab("IPF effect") + xlim(-6, 6)
}

# pval vs rsquared
{
  ggplot(dat, aes(x = -log10(pval), y = Rsquared, 
                  label = variables)) +
    geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) + 
    geom_text_repel(nudge_x = 0.02, segment.alpha = 0.3) + 
    theme(legend.title = element_blank())
}

# rsquared vs bleo effect (Estimate)

{
  ggplot(dat, aes(x = Rsquared, y = Estimate, 
                  label = variables)) +
    geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) + 
    geom_text_repel(
      data = dat[dat$Rsquared >= 0.5, ],
      segment.alpha =  0.3) + 
    geom_vline(xintercept = 0.95, alpha = 0.5, linetype = "dashed", color = "blue") +
    theme(legend.title = element_blank()) +
    annotate("text", x = 0.85, y = 5.5, color = "blue", alpha = 0.5,
             label = "R^2 == 0.95", parse = TRUE) +
    ylab("IPF Effect")
}



############# XUE MODULES ###########
{
  xue <- read.csv("data/XueModules/genesets_subset.gmt", sep = "\t", 
                  na.strings = "", header = F, stringsAsFactors = F)
  
  m1.1 <- xue[1, 3:ncol(xue)] %>% str_to_title(); m1.1 <- m1.1[!is.na(m1.1)]
  m1.1 <- m1.1[m1.1 %in% rownames(genexp)]
  m1.2 <- xue[2, 3:ncol(xue)] %>% str_to_title(); m1.2 <- m1.2[!is.na(m1.2)]
  m1.2 <- m1.2[m1.2 %in% rownames(genexp)]
  m1.3 <- xue[3, 3:ncol(xue)] %>% str_to_title(); m1.3 <- m1.3[!is.na(m1.3)]
  m1.3 <- m1.3[m1.3 %in% rownames(genexp)]
  m2.1 <- xue[4, 3:ncol(xue)] %>% str_to_title(); m2.1 <- m2.1[!is.na(m2.1)]
  m2.1 <- m2.1[m2.1 %in% rownames(genexp)]
  m2.2 <- xue[5, 3:ncol(xue)] %>% str_to_title(); m2.2 <- m2.2[!is.na(m2.2)]
  m2.2 <- m2.2[m2.2 %in% rownames(genexp)]
  m2.3 <- xue[6, 3:ncol(xue)] %>% str_to_title(); m2.3 <- m2.3[!is.na(m2.3)]
  m2.3 <- m2.3[m2.3 %in% rownames(genexp)]
  
  # Define dataframe of comparisons
  set1 <- list("Module7" = m1.1, "Module7" = m1.1, "Module7" = m1.1,
               "Module8" = m1.2, "Module8" = m1.2, "Module8" = m1.2,
               "Module9" = m1.3, "Module9" = m1.3, "Module9" = m1.3)
  set2 <- list("Module13" = m2.1, "Module14" = m2.2, "Module15" = m2.3,
               "Module13" = m2.1, "Module14" = m2.2, "Module15" = m2.3,
               "Module13" = m2.1, "Module14" = m2.2, "Module15" = m2.3)
  df <- data.frame(Set1 = I(set1), Set2 = I(set2))
  rownames(df) <- c("m7m13", "m7m14", "m7m15",
                    "m8m13", "m8m14", "m8m15",
                    "m9m13", "m9m14", "m9m15")
}


# calculate mPFR

{
  
  mpfr.all <- mPFR(df, genexp, pheno$Source, pheno$Treatment)
  mpfr.all.long <- melt(mpfr.all, id.vars = c("source", "treatment"))
}


#### EDA


# QQ plots
ggplot(mpfr.all.long, aes(sample = value)) +
  facet_wrap(~ variable) +
  ggplot2::stat_qq() + ggplot2::stat_qq_line(color = "red")
# Overall Density
ggplot(mpfr.all.long, aes(value)) +
  geom_density() + xlim(-10,10)


# PCA

"
  Falha completamente em separar controles de IPF
"
autoplot(prcomp(mpfr.all[,1:9], scale. = T), 
         data = mpfr.all, colour = "treatment", 
         size = 3)


# Boxplots
{
  ivars <- 1; nvars <- 9
  filtro = mpfr.all.long$variable %in% unique(mpfr.all.long$variable)[ivars:nvars]
  form <- value ~ treatment
  plot_signif(mpfr.all.long, form,
              filtro, outTest = T,  method = "wilcox.test",
              point = T)
  
}

# Modelos

{
  dat <- fit_mPFR(mpfr.all, Nvars = 9) 
  dat$variables <- rownames(dat)
}

# volcano

{
  ggplot(dat, aes(x = Estimate, y = -log10(pval), 
                  label = variables, size = Rsquared)) +
    geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) + 
    labs(size = "R squared", color = NULL) +
    geom_text_repel(force = T) +
    xlab("IPF effect") + xlim(-6, 6)
}

# pval vs rsquared
{
  ggplot(dat, aes(x = -log10(pval), y = Rsquared, 
                  label = variables)) +
    geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) + 
    geom_text_repel(nudge_x = 0.02, segment.alpha = 0.3) + 
    theme(legend.title = element_blank())
}

# rsquared vs bleo effect (Estimate)

{
  ggplot(dat, aes(x = Rsquared, y = Estimate, 
                  label = variables)) +
    geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) + 
    geom_text_repel(
      data = dat[dat$Rsquared >= 0.5, ],
      segment.alpha =  0.3) + 
    geom_vline(xintercept = 0.95, alpha = 0.5, linetype = "dashed", color = "blue") +
    theme(legend.title = element_blank()) +
    annotate("text", x = 0.85, y = 5.5, color = "blue", alpha = 0.5,
             label = "R^2 == 0.95", parse = TRUE) +
    ylab("IPF Effect")
}



###### TO-DO LIST ########

"
  1) Buscar ortologos dos genes dos ratos para humano
  2) Testar outros moduloes do Xue2014
  3) Mostrar o N nas facets dos boxplots das variables
  4) Modelar correlacao entre variables no humanos e nos ratos
"
{
  # Convert rat to human gene names
  
  convert_rat_gene_list <- function(x){
    
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    mouse_matches = getLDS(attributes = c("mgi_symbol"), 
                           filters = "mgi_symbol", 
                           values = x , mart = rat, 
                           attributesL = c("hgnc_symbol"), 
                           martL = human, uniqueRows=T)
    no_matches = setdiff(x, mouse_matches[,1]) 
    
    x = data.frame(HGNC.symbol = x)
    df <- merge(x, mouse_matches, by = "HGNC.symbol")
    
    return(list(df, data.frame(no_matches)))
  }
  
  genes <- convert_rat_gene_list(c1not)
  mymatches <- genes[[1]]
  nomathes <- genes[[2]]
  nrow(nomathes) # over 2 thousand genes with no matches
  nrow(mymatches)
}

}