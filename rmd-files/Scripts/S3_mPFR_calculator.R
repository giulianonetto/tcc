{
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(reshape2)
  library(ggpubr)
  library(ggloop)
  library(ggrepel)
  library(qqplotr)
  library(rafalib)
  library(maSigPro)
  library(Biobase)
  library(car)
  library(yhat)
  
  
  setwd("/home/giulianonetto/windows/tcc/rmd-files")
  eset <- readRDS("data/GSE48455/suppdata/ExpressionSetClean.rds")
  pheno <- pData(eset)
  genexp <- exprs(eset)
  rownames(genexp) <- rownames(genexp) %>% str_to_title()
}

mPFR <- function(sets_data.frame, exprs_matrix, Time, Treat) {
  # rownames(exprs_matrix) <- toupper(rownames(exprs_matrix)) # normalize cases
  mpfr_data.frame <- data.frame()
  # rownames(mpfr_data.frame) <- colnames(exprs_matrix)
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
  mpfr_data.frame$times <- Time
  mpfr_data.frame$treatment <- Treat
  return(mpfr_data.frame)
}


# define genesets from masigpro clusters
{
  hcluster <- readRDS("data/GSE48455/suppdata/hcluster_masigpro_bauer2015.rds")
  c1a = names(hcluster$cut)[which(hcluster$cut == 1)] %>% gsub(pattern = "_", replacement = ".") %>% str_to_title()
  c1 <- c1a[c1a %in% rownames(genexp)]; print(length(c1)/length(c1a) * 100)
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

  # build metadata df from genexp colnames

  times <- str_extract(colnames(genexp), "[0-9]+") %>% as.numeric()
  treatment <- ifelse(grepl("BLEO", colnames(genexp)), "Bleomycin", "Control")
  colors <- ifelse(treatment == "Bleomycin", "red", "blue")
  meta <- data.frame(arrays = colnames(genexp), 
                     times = times, 
                     treatment = treatment, 
                     colors = colors, 
                     stringsAsFactors = F) %>% arrange(times, treatment)
  # meta[which(meta$times == 0), 3] <- "Untreated"
  # meta[which(meta$times == 0), 4] <- "gray40"
}

{
  # Define dataframe of comparisons
  s1l <- list("S1" = c("Il12b"), 
              "S2" = c("Ccl12", "Ccr7", "Ccl22", "Ccl24"),
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
              "S2" = c("Il13ra1.3", "Il13ra2", "Cd36", "Cxcr1", "St7"),
              "S3" = c("Apoe", "Mmp14", "Lbp", "Il1rl1.1"), 
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
  p <- ggplot(df, aes(x = treatment, y = value))
  stat <- stat_compare_means(comparisons = list(c("bleomycin", "control")),
                             method = method,
                             aes(label = ..p.signif..))
  fac <- facet_wrap(~ variable)
  title <- ggtitle(titles)
  if (point == FALSE) {
    p + stat + fac + geom_boxplot() + title
  }
  else {
    p + stat + fac + geom_point() + geom_boxplot(alpha = 0.25) + title
  }
}



# Calculate mPFR

{

  mpfr.all <- mPFR(df, genexp, pheno$Times, pheno$Cy3)
  mpfr.all.long <- melt(mpfr.all, id.vars = c("times", "treatment"))
  mpfr.all$times <- factor(mpfr.all$times)
  mpfr.all$treatment <- factor(mpfr.all$treatment, levels = c("bleomycin", "control"))
  contrasts(mpfr.all$treatment) <- c("bleomycin" = 1, "control" = 0)
}


# EDA


# QQ plots
ggplot(mpfr.all.long, aes(sample = value)) +
  facet_wrap(~ variable) +
  ggplot2::stat_qq() + ggplot2::stat_qq_line(color = "red")
# Overall Density
ggplot(mpfr.all.long, aes(value)) +
  geom_density() + xlim(-10,10)


# PCA
autoplot(prcomp(mpfr.all[,1:28], scale. = T), 
         data = mpfr.all, colour = "treatment", 
         size = 2)
autoplot(prcomp(mpfr.all[,1:28], scale. = T), 
         data = mpfr.all, colour = "times", 
         size = 2, shape = "treatment")


# Boxplots
{
  day = 7
  ivars <- 25; nvars <- 30
  filtro = mpfr.all.long$variable %in% unique(mpfr.all.long$variable)[ivars:nvars] &
    mpfr.all.long$times %in% day
  form <- value ~ treatment
  
  plot_signif(mpfr.all.long, form,
              filtro, outTest = T,  method = "wilcox.test",
              titles = str_glue("Day {day}"), point = T)
  
}


# Modelos

fit_mPFR <- function(df, day = 3, Nvars = 28){
  fitdf <- data.frame(Rsquared = vector('numeric'),
                      pval = vector('numeric'),
                      Estimate = vector('numeric'), 
                      stringsAsFactors = F)
  df = df[df$times == day, ]
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

{
day = 7
dat <- fit_mPFR(mpfr.all, day = day) 
dat$variables <- rownames(dat)
}

# volcano
ggplot(dat, aes(x = Estimate, y = -log10(pval), label = variables)) +
  geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) + 
  theme(legend.title = element_blank()) +
  geom_text_repel(force = T) +
  xlab("Bleomycin Effect")

# pval vs rsquared
ggplot(dat, aes(x = -log10(pval), y = Rsquared, 
                label = variables)) +
  geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) + 
  geom_text_repel(nudge_x = 0.02, segment.alpha = 0.3) + 
  theme(legend.title = element_blank())

# rsquared vs bleo effect (Estimate)
ggplot(dat, aes(x = Rsquared, y = Estimate, 
                label = variables)) +
  geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) + 
  geom_text_repel(
    data = dat[dat$Rsquared >= 0.5, ],
    segment.alpha =  0.3) + 
  geom_vline(xintercept = 0.95, alpha = 0.5, linetype = "dashed", color = "blue") +
  theme(legend.title = element_blank()) +
  annotate("text", x = 0.91, y = 5.5, color = "blue", alpha = 0.5,
           label = "R^2 == 0.95", parse = TRUE) +
  ylab("Bleomycin Effect")



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
  
  mpfr.all <- mPFR(df, genexp, pheno$Times, pheno$Cy3)
  mpfr.all.long <- melt(mpfr.all, id.vars = c("times", "treatment"))
  mpfr.all$times <- factor(mpfr.all$times)
  mpfr.all$treatment <- factor(mpfr.all$treatment, levels = c("bleomycin", "control"))
  contrasts(mpfr.all$treatment) <- c("bleomycin" = 1, "control" = 0)
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
autoplot(prcomp(mpfr.all[,1:9], scale. = T), 
         data = mpfr.all, colour = "treatment", 
         size = 2)
autoplot(prcomp(mpfr.all[,1:9], scale. = T), 
         data = mpfr.all, colour = "times", 
         size = 2, shape = "treatment")

# Boxplots
{
  day = 3
  ivars <- 1; nvars <- 9
  filtro = mpfr.all.long$variable %in% unique(mpfr.all.long$variable)[ivars:nvars] &
    mpfr.all.long$times %in% day
  form <- value ~ treatment
  
  plot_signif(mpfr.all.long, form,
              filtro, outTest = T,  method = "t.test",
              titles = str_glue("Day {day}"), point = T)
  
}

# longituginal profile of module8 vs module13
ggplot(mpfr.all, aes(x = as.numeric(as.character(times)), y = m8m13, 
                     group = treatment, color = treatment)) +
  geom_point() + geom_smooth(method = "lm", formula = y ~ poly(x, 3))



# Modelos


{
  day = 3
  dat <- fit_mPFR(mpfr.all, day = day, Nvars = 9) 
  dat$variables <- rownames(dat)
}

# volcano
ggplot(dat, aes(x = Estimate, y = -log10(pval), label = variables)) +
  geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) + 
  theme(legend.title = element_blank()) +
  geom_text_repel(force = T) +
  xlab("Bleomycin Effect")

# pval vs rsquared
ggplot(dat, aes(x = -log10(pval), y = Rsquared, 
                label = variables)) +
  geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) + 
  geom_text_repel(nudge_x = 0.02, segment.alpha = 0.3) + 
  theme(legend.title = element_blank())

# rsquared vs bleo effect (Estimate)
ggplot(dat, aes(x = Rsquared, y = Estimate, 
                label = variables)) +
  geom_point(aes(color = ifelse(dat$pval < 0.05, "Significant", "Not significant"))) + 
  geom_text_repel(
    data = dat[dat$Rsquared >= 0.5, ],
    segment.alpha =  0.3) + 
  geom_vline(xintercept = 0.95, alpha = 0.5, linetype = "dashed", color = "blue") +
  theme(legend.title = element_blank()) +
  annotate("text", x = 0.88, y = 5.5, color = "blue", alpha = 0.5,
           label = "R^2 == 0.95", parse = TRUE) +
  ylab("Bleomycin Effect")
