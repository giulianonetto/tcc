library(dplyr)
library(reshape2)
library(ggplot2)
library(nlme)
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
arg2 <- args[1]

il12b <- grep("IL12.*", rownames(genexp))[1]

bleo_arrays <- pheno %>% filter(Cy3 == "bleomycin" | Times == 0) %>%
  arrange(Times) %>% dplyr::select(Rows, Times)
control_arrays <- pheno %>% filter(Cy3 == "control") %>%
  arrange(Times) %>% dplyr::select(Rows, Times)


# Arg1

bleo_arg1 <- data.frame(ARG1 = genexp[arg1,bleo_arrays[,1]])
bleo_arg1 <- cbind(bleo_arg1, bleo_arrays)
bleo_arg1 %>% group_by(Times) %>% summarise(ARG1 = median(ARG1)) %>% 
  ggplot(aes(x = Times, y = ARG1))+geom_point()+
  stat_smooth(method="lm", se = TRUE,
              formula=y ~ poly(x, 3, raw=TRUE),colour="red")+
  labs(title="Arg1 (bleo) - Bauer2015")

control_arg1 <- data.frame(ARG1 = genexp[arg1,control_arrays[,1]])
control_arg1 <- cbind(control_arg1, control_arrays)
control_arg1 %>% group_by(Times) %>% summarise(ARG1 = mean(ARG1)) %>% 
  ggplot(aes(x = Times, y = ARG1))+geom_point()+
  stat_smooth(method="lm", se=TRUE, fill=NA,
              formula=y ~ poly(x, 3, raw=TRUE),colour="blue")+
  labs(title="Arg1 (cntrl) - Bauer2015")


# IL12b

bleo_il12b <- data.frame(IL12B = genexp[il12b,bleo_arrays[,1]])
bleo_il12b <- cbind(bleo_il12b, bleo_arrays)
bleo_il12b %>% group_by(Times) %>% summarise(IL12B = median(IL12B)) %>% 
  ggplot(aes(x = Times, y = IL12B))+geom_point()+
  stat_smooth(method="lm", se=TRUE,
              formula=y ~ poly(x, 3, raw=TRUE),colour="red")+
  labs(title="IL12b (bleo) - Bauer2015")

control_il12b <- data.frame(IL12B = genexp[il12b,control_arrays[,1]])
control_il12b <- cbind(control_il12b, control_arrays)
control_il12b %>% group_by(Times) %>% summarise(IL12B = mean(IL12B)) %>% 
  ggplot(aes(x = Times, y = IL12B))+geom_point()+
  stat_smooth(method="lm", se=TRUE,
              formula=y ~ poly(x, 3, raw=TRUE),colour="blue")+
  labs(title="IL12b (cntrl) - Bauer2015")

# (IL12B/ARG1) * (median(ARG1)/median(IL12B)) BLEOMYCIN GROUP
df <- cbind(bleo_arg1$ARG1,bleo_il12b)
df <- df %>%
  mutate(Ratio = IL12B / `bleo_arg1$ARG1`, 
         Correction =  median(df$`bleo_arg1$ARG1`) / median(df$IL12B)) %>%
  mutate(PFR = Ratio * Correction) %>% group_by(Times) %>% mutate(PFRm = median(PFR)) %>% as.data.frame()
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
    fit <- lm(PFR ~ poly(as.numeric(Times), 3, raw = T) * Treatment, data = df)
    outlierTest(fit) %>% print()
    qqPlot(fit, labels=row.names(df), 
           id.method = "identify", simulate = T, 
           main = str_glue("Q-Q Plot {Main}"))
    ResidualsPlot(fit, Main)
  }

}
outlierDetec(df, Main = "All") # 29 is outlier!
df <- df[-grep("29", rownames(df)),] 
outlierDetec(df, Main = "Without 29") # That's it!

df %>% ggplot(aes(x = Times, y = PFR))+
  geom_point(color = "grey40", alpha = 0.6)+
  stat_smooth(aes(x = Times, y = PFR),
              method = "lm", formula = y ~ poly(x, 3, raw = T))+
  ggtitle("Median PFR")

ggsave("data/GSE48455/suppdata/PFR_median_with_obs.png", width = 25, height = 20, units = "cm")

# (IL12B/ARG1) * (median(ARG1)/median(IL12B)) - CONTROLS

dfc <- cbind(control_arg1$ARG1, control_il12b)
dfc <- dfc %>%
  mutate(Ratio = IL12B / `control_arg1$ARG1`,
         Correction =  median(dfc$`control_arg1$ARG1`) / median(dfc$IL12B)) %>%
  mutate(PFR = Ratio * Correction) %>% group_by(Times) %>% mutate(PFRm = median(PFR))

outlierDetec(dfc) # 7 is outlier!
dfc <- dfc[-7,]
outlierDetec(dfc) # 10 is outlier!
dfc <- dfc[-10,]
outlierDetec(dfc) # 27 is outlier!
dfc <- dfc[-27,]
outlierDetec(dfc) # that's it!

dfc %>% ggplot(aes(x = Times, y = PFR))+geom_point()+
  stat_smooth(aes(x = Times, y = PFRm), method = "lm", 
              formula = y ~ poly(x, 3, raw = T))+ggtitle("Median PFR - control")






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
# 
# df.bleo <- df %>% filter(Treatment == "bleomycin")
# df.ctrl <- df %>% filter(Treatment == "control")
# df.bleo %>% as.data.frame() %>% filter(Times == 14) %>% select(PFR) %>% unlist()
# df.bleo <- df.bleo[-c(15,18,29,33),] %>% as.data.frame()
# df.ctrl <- df.ctrl %>% as.data.frame()

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

anova(baseline, treats, treats_and_times, final_fit) # only treatment significant

bartlett.test(PFR ~ Treatment, data = df)
t.test(PFR ~ Treatment, data = df, var.equal = TRUE, conf.level=0.95)
wilcox.test(PFR ~ Treatment, data = df, conf.level = 0.95)
wilcox_test(PFR ~ Treatment, data = df, distribution = "exact")

dftoplot <- df
dftoplot$Times <- dftoplot$Times %>% levels %>% as.numeric()
{
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
