library(tidyverse)
library(stringr)
library(TcGSA)
library(biomaRt)
library(dplyr)
library(GEOquery)
# Load Data
setwd("/home/giulianonetto/windows/tcc/rmd-files")
eset <- readRDS("data/GSE48455/suppdata/ExpressionSetClean.rds")
Xue.gmt <- GSA::GSA.read.gmt("data/XueModules/genesets.gmt")
genesets <- read.csv("data/XueModules/genesets.gmt", 
                     sep = "\t", header = F, row.names = 1)

is.na(genesets) <- genesets == ""

# Build grouping indeces

pheno <- pData(eset)
pheno$Rows <- rownames(pheno)
healthy <- pheno$Phase == "Healthy"
injury <-  pheno$Phase == "Injury"
fibrosis <- pheno$Phase == "Early Fibrosis" | pheno$Phase == "Late Fibrosis"
healing <- pheno$Phase == "Healing"
control <- pheno$Cy3 == "control"
bleomycin <- pheno$Cy3 == "bleomycin"



# Expression Matrix

genexp <- exprs(eset)

## BiomaRt Annotation fixes

not_found = vector("list", length = 49)
for (i in 1:nrow(genesets)){
  nas <- is.na(genesets[i,2:ncol(genesets)])
  query <- genesets[i,2:ncol(genesets)][!nas]
  print(str_glue("Query length: {length(query)}"))
  matches <- query %in% rownames(genexp)
  unmatched <- query[!matches]
  perc = round(length(unmatched)/length(query) * 100)
  if (perc1 == perc2) {
  print(str_glue("Unmatched: {length(unmatched)} ({perc}%)"))
  }
  if (perc < 100){
    not_found[[i]] <- unmatched
  }
}
not_found <- ldply(not_found, rbind)



# Ploting the Polarization Factor Ratio (PFR) by [@Buscher2017]

args <- grep("^ARG.*[0-9]*", rownames(genexp))
arg1 <- args[3]
arg2 <- args[1]

il12b <- grep("IL12.*", rownames(genexp))[1]

bleo_arrays <- pheno %>% filter(Cy3 == "bleomycin" | Times == 0) %>%
  arrange(Times) %>% select(Rows, Times)
control_arrays <- pheno %>% filter(Cy3 == "control") %>%
  arrange(Times) %>% select(Rows, Times)


# Arg1

bleo_arg1 <- data.frame(ARG1 = genexp[arg1,bleo_arrays[,1]])
bleo_arg1 <- cbind(bleo_arg1, bleo_arrays)
bleo_arg1 %>% group_by(Times) %>% summarise(ARG1 = median(ARG1)) %>% 
  ggplot(aes(x = Times, y = ARG1))+geom_point()+
  stat_smooth(method="lm", se=TRUE, fill=NA,
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

# (IL12B/ARG1) * (median(ARG1)/median(IL12B))
df <- cbind(bleo_arg1$ARG1,bleo_il12b)
df <- df %>%
  mutate(Ratio = IL12B / `bleo_arg1$ARG1`, 
         Correction =  median(df$`bleo_arg1$ARG1`) / median(df$IL12B)) %>%
  mutate(PFR = Ratio * Correction)

df %>% group_by(Times) %>% summarise(PFR = median(PFR)) %>% ggplot(aes(x = Times, y = PFR))+
  geom_point()+stat_smooth(method = "lm", formula = y ~ poly(x, 3, raw = T))+ggtitle("Median PFR")

ggsave("data/GSE48455/suppdata/PFR_median.png", width = 25, height = 20, units = "cm")

# (IL12B/ARG1) * (median(ARG1)/median(IL12B)) - CONTROLS

dfc <- cbind(control_arg1$ARG1, control_il2b)
# dfc <- dfc %>%
#   mutate(Ratio = IL12B / `control_arg1$ARG1`, 
#          Correction =  median(dfc$`control_arg1$ARG1`) / median(dfc$IL12B)) %>%
#   mutate(PFR = Ratio * Correction) %>% group_by(Times) %>% mutate(PFRmedian = median(PFR))
# 
# dfc[-c(7,11,23,29),] %>% ggplot()+
#   geom_point(aes(x = Times, y = PFR))+
#   stat_smooth(aes(x = Times, y = PFRmedian), 
#               method = "lm", formula = y ~ poly(x, 3, raw = T))+
#   ggtitle("Median PFR - control")
dfc <- dfc %>%
  mutate(Ratio = IL12B / `control_arg1$ARG1`,
         Correction =  median(dfc$`control_arg1$ARG1`) / median(dfc$IL12B)) %>%
  mutate(PFR = Ratio * Correction)

dfc[-c(7,11,23,29),] %>% group_by(Times) %>% summarise(PFR = median(PFR)) %>% ggplot(aes(x = Times, y = PFR))+
  geom_point()+stat_smooth(method = "lm", formula = y ~ poly(x, 3, raw = T))+ggtitle("Median PFR - control")



" O PFR parece nao ter distribuicao normal, alem dos dados terem muitos
outliers. Tem que comparar o controle com a bleo e fazer wilcoxon.
Tem que plotar com o teste. Tem que fazer o mesmo para os set de genes
do Xue. Tem que avaliar a possibilidade de subset do dataset e de nova
normalizacao dos PFRs. Decidir formula para PFR modular. Fazer os 
mesmos plots com os dados tratados do bauer, so pra checar
congruencia.

"

# DOING THE SAME FOR BAUER-NORMALIZED DATA


bauer <- getGEO("GSE48455")
bauer <- bauer$GSE48455_series_matrix.txt.gz
exp.bauer <- exprs(bauer)
pheno.bauer <- pData(bauer)
