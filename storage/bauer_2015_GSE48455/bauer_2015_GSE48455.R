library(stringr)
library(maSigPro)
"
Test with dregees 2, 3 and 4. min.obs = 20. step.method=backward
"



# Expression matrix and metadata - exclude time 56 and 42

metadata = read.csv("metadata_renamed.tsv", sep = "\t", header = T)
metadata[1:3,3] = c(0,0,0) # untreated are taken as starting point
rownames(metadata) <- metadata$Geo
edes <- metadata[1:52,3:6]
edes <- edes[,c(1,4,2,3)]

gene_exp <- read.table("gene_expression_matrix_mapped", sep = "\t", header = T, row.names = 1)
bigtimes <- grepl("T56|T42", names(gene_exp))
gene_exp <- gene_exp[,!bigtimes]
gene_exp = gene_exp[, match(rownames(edes), colnames(gene_exp))]


# make design matrix
design2 <- make.design.matrix(edesign = edes, degree = 2) # quadratic regression model
design3 <- make.design.matrix(edesign = edes, degree = 3) # cubic regression model
design4 <- make.design.matrix(edesign = edes, degree = 4) # cubic regression model

##### REGRESSION MODEL #####

fit2 <- p.vector(gene_exp, design2, Q = 0.05, MT.adjust = "BH", min.obs = 20) # quadratic
fit3 <- p.vector(gene_exp, design3, Q = 0.05, MT.adjust = "BH", min.obs = 20) # cubic
fit4 <- p.vector(gene_exp, design4, Q = 0.05, MT.adjust = "BH", min.obs = 20) # fourth degree

##### FIND SIGNIFICANT GENES #####

# step-wise regression
tstep2 <- T.fit(fit2, step.method = "backward", alfa = 0.05) # quadratic
tstep3 <- T.fit(fit3, step.method = "backward", alfa = 0.05) # cubic
tstep4 <- T.fit(fit4, step.method = "backward", alfa = 0.05) # fourth degree

# check how well fit the models are - using R-squared
(a <- mean(tstep2$sol[,"R-squared"]>0.5)*100) #  15.44% for degree=2
(b <- mean(tstep3$sol[,"R-squared"]>0.5)*100) # 31.7% for degree=3
(c <- mean(tstep4$sol[,"R-squared"]>0.5)*100) # 41.7% for degree=4
pdf("out_min.obs_20_back_bauer2.pdf")
barplot(c(15.44, 31.7, 41.7), main = c("Percentage of models with R^2 > 0.5"),
        names.arg = c("Dregee 2", "Dregee 3", "Dregee 4"), ylim = c(0, 60))
text(0.75, 10, as.character(a))
text(1.9, 20, as.character(b))
text(3.15, 30, as.character(c))
dev.off()

save.image()
