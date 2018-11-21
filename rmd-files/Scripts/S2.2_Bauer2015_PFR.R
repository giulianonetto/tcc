
{
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(stringr)
  library(nlme)
  library(car)
  library(coin)
  library(multcomp)
  library(MASS)
  library(ggpubr)
  library(ggsignif)
  setwd("/home/giulianonetto/windows/tcc/rmd-files")
  eset <- readRDS("data/GSE48455/suppdata/ExpressionSetClean.rds")
  pheno <- pData(eset)
  genexp <- exprs(eset)
  
  # Ploting the Polarization Factor Ratio (PFR) by [@Buscher2017]
  pheno$Rows <- rownames(pheno)
  args <- grep("^ARG.*[0-9]*", rownames(genexp))
  arg1 <- args[3]
  il12b <- grep("IL12.*", rownames(genexp))[1]
}

{
    # -> Outliers identification functions!
    "
    Studentized residuals with Bonferonni p < 0,05:
    The car::outlierTest() function reports the Bonferroni adjusted p-value for the largest
    absolute studentized residual from your fit [@RobertKabacoff - R in Action].
    "
  ResidualsPlot <- function(Fit, title = ""){
    res <- rstudent(Fit)
    Residuals <- data.frame(obs = names(res), values = res)
    h = 2*sd(res)
    ggplot(Residuals, aes(reorder(obs, values), values))+
      geom_bar(stat = "identity")+
      theme(text = element_text(family = "Decima WE", color = "grey20"),
            axis.text.x = element_text(angle = 90, hjust = 1))+
      ggtitle(str_glue("Studenized Residuals {title}"))+
      geom_hline(yintercept=c(h, -h),
                 linetype="dashed", color = "red")+
      geom_text(aes(0,h,label = str_glue("2 SD's ({round(h, 2)})"),
                    hjust = -0.1, vjust = -1.15,
                    family = "Courier"), color = "grey40")
  }
  outlierDetec <- function(df, Main = ""){
    df <- as.data.frame(df)
    rownames(df) <- paste0("Obs ", rownames(df))
    fit <- lm(PFR ~ poly(as.numeric(as.character(Times)), 3, raw = T) * Treatment, data = df)
    outlierTest(fit) %>% print()
    qqPlot(fit, labels = row.names(mydat),
           id.method = "identify", simulate = T,
           main = str_glue("Q-Q Plot {Main}"))
    ResidualsPlot(fit, Main)
  }
}

" O PFR parece nao ter distribuicao normal, alem dos dados terem muitos
outliers. Tem que comparar o controle com a bleo e fazer wilcoxon.
Tem que plotar com o teste. Tem que fazer o mesmo para os set de genes
do Xue. Tem que avaliar a possibilidade de subset do dataset e de nova
normalizacao dos PFRs. Decidir formula para PFR modular. Fazer os 
mesmos plots com os dados tratados do bauer, so pra checar
congruencia.

"
# ANOVA fitting to PFR

# build long-formatted df to compare PFR
{
df <- data.frame(IL12B = genexp[il12b,], ARG1 = genexp[arg1,])
df$Arrays <- rownames(df)
df$Times <- pheno$Times
df$Treatment <- pheno$Cy3
df <- df[,c(3,1,2,4,5)]; rownames(df) <- 1:nrow(df)
# we replicate 0 timepoints for both treatment groups
untreated <- df %>% filter(Times == 0)
untreated$Treatment <- "bleomycin"
df <- rbind(df, untreated)
str(df)
# convert times and treatment to factors
df$Times <- factor(df$Times); df$Treatment <- factor(df$Treatment)
# calculate PFR 
df <- df %>% arrange(Treatment, Times) %>%
  mutate(Ratio = IL12B / ARG1,
         Correction =  mean(df$ARG1) / mean(df$IL12B), 
         Correction.Robust = median(df$ARG1) / median(df$IL12B)) %>%
  mutate(PFR = Ratio * Correction, 
         PFR.Robust = Ratio * Correction.Robust) %>%
  group_by(Times, Treatment) %>%
  mutate(PFRm = mean(PFR), PFRm.Robust = median(PFR.Robust))
table(df$Times)
saveRDS(df, "data/GSE48455/suppdata/anova_table.rds")
df <- readRDS("data/GSE48455/suppdata/anova_table.rds")

outlierDetec(df, Main = "Poly(Times) * Treatment") # 41 is outlier
df <- df[-41,]
outlierDetec(df, Main = "Poly(Times) * Treatment") # 44 is outlier
df <- df[-44,]
outlierDetec(df, Main = "Poly(Times) * Treatment") # 61 is outlier
df <- df[-61,]
outlierDetec(df, Main = "Poly(Times) * Treatment") # 55 is outlier
df <- df[-55,]
outlierDetec(df, Main = "Poly(Times) * Treatment") # 62 is outlier
df <- df[-62,]
outlierDetec(df, Main = "Poly(Times) * Treatment") # 54 is outlier
df <- df[-54,]
outlierDetec(df, Main = "Poly(Times) * Treatment") # That's it!

} 


{
  df$Rat <- paste0(toupper(str_extract(df$Treatment, "\\w")), 
                   str_extract(str_split(df$Arrays, "\\.",simplify = T)[,2],"[0-9]+"))
  BleovsCtrl <- c(1,0)
  contrasts(df$Treatment) <- BleovsCtrl
  contrasts(df$Times) <- contrasts(df$Times)
  
  baseline <- lme(PFR ~ 1, random = ~1|Rat/Times, 
                  data = df, method = "ML")
  treats <- update(baseline, .~. + Treatment)
  treats_and_times <- update(treats, .~. + Times)
  final_fit <- update(treats_and_times, .~. + Treatment:Times)
  
  # only treatment significant
  anova(baseline, treats, treats_and_times, final_fit)
  
  bartlett.test(PFR ~ Treatment, data = df)
  t.test(PFR ~ Treatment, data = df, var.equal = TRUE, conf.level=0.95)
  wilcox.test(PFR ~ Treatment, data = df, conf.level = 0.95)
  wilcox_test(PFR ~ Treatment, data = df, distribution = "exact")
}

{
  dftoplot <- df
  dftoplot$Times <- dftoplot$Times %>% levels %>% as.numeric()

  p1 <- dftoplot %>% ggplot(aes(x = Times, y = PFR, group = Treatment)) +
          geom_point(aes(color = Treatment), size = 3, alpha = 0.3) +
          geom_smooth(aes(color = Treatment), method = "glm")+
          xlab("Days after treatment")+
          theme(text = element_text(size = 20))
  
  p2 <- dftoplot %>% ggplot(aes(x = Treatment, y = PFR, group = Treatment)) +
          geom_boxplot(aes(color = Treatment)) +
          geom_smooth(aes(color = Treatment), method = "lm") +
          # stat_compare_means(comparisons = list(c("bleomycin", "control")))
          geom_signif(comparisons = list(c("bleomycin", "control")),
                      map_signif_level = TRUE, tip_length = 0) +
          theme(text = element_text(size = 20))
  ggarrange(p1, p2, nrow = 1, ncol = 2, 
            common.legend = T, legend = "right") }


# Build PFR for other genes

get_pfr_table <- function(exp, pheno, m1, m2){
  mydat <- data.frame(M1 = exp[m1,], M2 = exp[m2,])
  mydat$Arrays <- rownames(mydat)
  mydat$Times <- pheno$Times
  mydat$Treatment <- pheno$Cy3
  mydat <- mydat[,c(3,1,2,4,5)]; rownames(mydat) <- 1:nrow(mydat)
  mydat <- mydat %>% arrange(Treatment, Times) %>%
    mutate(Ratio = M1 / M2,
           Correction =  mean(mydat$M2) / mean(mydat$M1), 
           Correction.Robust = median(mydat$M2) / median(mydat$M1)) %>%
    mutate(PFR = Ratio * Correction, 
           PFR.Robust = Ratio * Correction.Robust) %>%
    group_by(Times, Treatment) %>%
    mutate(PFRm = mean(PFR), PFRm.Robust = median(PFR.Robust))
  mydat$Times <- factor(mydat$Times); mydat$Treatment <- factor(mydat$Treatment)
  return(mydat)
}

outlierDetec <- function(mydat, Main = ""){
  mydat <- as.data.frame(mydat)
  rownames(mydat) <- paste0("Obs ", rownames(mydat))
  fit <- lm(PFR ~ poly(as.numeric(as.character(Times)), 3, raw = T) * Treatment,
            data = mydat)
  outlierTest(fit) %>% print()
  qqPlot(fit, id.method = "identify", simulate = T,
         main = str_glue("Q-Q Plot {Main}"))
  ResidualsPlot(fit, Main)
}

# relevant genes - test for outliers
{
inos <- grep("^NOS2.*", rownames(genexp))
genexp[inos[1],] <- apply(genexp[inos,], 2, median) 
mydat <- get_pfr_table(genexp, pheno = pheno, inos[1], arg1)
outlierDetec(mydat = mydat, Main = "all")
mydat <- mydat[-38,]
outlierDetec(mydat = mydat, Main = "-38")
mydat <- mydat[-41,]
outlierDetec(mydat = mydat, Main = "-41")
mydat <- mydat[-c(58,52),]
outlierDetec(mydat = mydat, Main = "-58 e 52")
mydat <- mydat[-c(53,39),]
outlierDetec(mydat = mydat, Main = "-53 e 39")
mydat <- mydat[-57,]
outlierDetec(mydat = mydat, Main = "-57")
}

mydat$Rat <- paste0(toupper(str_extract(mydat$Treatment, "\\w")), 
                 str_extract(str_split(mydat$Arrays, "\\.",simplify = T)[,2],"[0-9]+"))
BleovsCtrl <- c(1,0)
contrasts(mydat$Treatment) <- BleovsCtrl
contrasts(mydat$Times) <- contrasts(mydat$Times)
print(contrasts(mydat$Treatment))
print(contrasts(mydat$Times))
baseline <- lme(PFR ~ 1, random = ~1|Rat/Times, 
                data = mydat, method = "ML")
treats <- update(baseline, .~. + Treatment)
treats_and_times <- update(treats, .~. + Times)
interaction_treat_and_times <- update(treats_and_times, .~. + Treatment:Times)

anova(baseline, treats, treats_and_times)  # only treatment significant
bartlett.test(PFR ~ Treatment, data = mydat)
t.test(PFR ~ Treatment, data = mydat, var.equal = TRUE, conf.level=0.95)
wilcox.test(PFR ~ Treatment, data = mydat, conf.level = 0.95)
wilcox_test(PFR ~ Treatment, data = mydat, distribution = "exact") 
# p = 4,335e-07 for PFR = NOS2/arg1 3,932e-07

dftoplot <- mydat
dftoplot$Times <- dftoplot$Times %>% as.character() %>% as.numeric()
{
  p3 <- dftoplot %>% ggplot(aes(x = Times, y = PFR, group = Treatment)) +
    geom_point(aes(color = Treatment), size = 3, alpha = 0.3) +
    geom_smooth(aes(color = Treatment), method = "glm")+
    xlab("Days after treatment")+
    theme(text = element_text(size = 20))
  
  p4 <- dftoplot %>% ggplot(aes(x = Treatment, y = PFR, group = Treatment)) +
    geom_boxplot(aes(color = Treatment)) +
    geom_smooth(aes(color = Treatment), method = "lm") +
    # stat_compare_means(comparisons = list(c("bleomycin", "control")))
    geom_signif(comparisons = list(c("bleomycin", "control")),
                map_signif_level = TRUE, tip_length = 0) +
    theme(text = element_text(size = 20))

}

# Set common legends
{
p1b <- p1 + theme(axis.title.x = element_blank()) + 
            ylab("PFR (IL12b / Arg1)")
p2b <- p2 + theme(axis.title = element_blank())
p3b <- p3 + ylab("PFR (iNOS / Arg1)")
p4b <- p4 + theme(axis.title.y = element_blank())

ggarrange(p1b, p2b, p3b, p4b, nrow = 2, ncol = 2, common.legend = T, 
          legend = "right") + xlab("")
}
# ** is p=0.007983 and *** is p=4,335e-07

