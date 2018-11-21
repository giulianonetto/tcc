setwd("/home/giulianonetto/windows/tcc/storage/bauer_2015_GSE48455")
load(".RData")
library(maSigPro)
library(rafalib)

pdf('il1b.pdf')
il1b <- gene_exp[rownames(gene_exp)=='Il1b',]
counter = 2
for (des in list(design2, design3, design4)){
  PlotGroups(il1b, edesign = edes, show.fit = T,
             dis = des$dis, main = counter, 
             groups.vector = des$groups.vector,
            widt  )
  counter = counter +1
}
dev.off()


###### RPUBS BY NEIL SAUNDERS 2015 - http://rpubs.com/neilfws/83863

sigs3 <- get.siggenes(tstep3, rsq = 0.6, vars = "groups")

names(sigs3)
names(sigs3$sig.genes)
names(sigs3$sig.genes$Control)

#We can get the data frames with p-values for control, celecoxib- and rofecoxib-treated cells.
control <- sigs3$sig.genes$Control$sig.pvalues
bleomycin <- sigs3$sig.genes$BleomycinvsControl$sig.pvalues

##### Plotting timecourses for genes of interest #####


## 1 Control samples

head(control[order(control$`p-value`, decreasing = FALSE), ])

plotGenes <- function(e, g, md) {
  require(ggplot2)
  d <- e[g, ] %>% 
    as.matrix() %>% t() %>% as.data.frame() %>% 
    setNames("value") %>% 
    mutate(Rep = md$Replicate,
           time = md$Time,
           agent = md$Bleomycin)
  gg <- ggplot(d, aes(time, value)) + geom_boxplot(aes(position = factor(time)), 
                                                   outlier.shape = NA) + theme_bw() + scale_x_continuous(breaks = unique(d$time)) + 
    geom_jitter(aes(color = factor(agent))) + geom_smooth() + labs(title = g, x = "time (hours)", y = "RMA value") + scale_color_discrete(name = "treatment")
  return(gg)
}


g <- rownames(control[order(control$`p-value`, decreasing = FALSE), ])[1]

plotGenes(gene_exp, g, edes) 

## 4.2 Treated samples

g1 <- rownames(bleomycin[order(bleomycin$`p-value`, decreasing = FALSE), ])[1]
plotGenes(gene_exp, g1, edes) 

# IL1b

plotGenes(gene_exp, "inos", edes)
