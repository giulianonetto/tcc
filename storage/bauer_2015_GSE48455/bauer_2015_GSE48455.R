########### DOWNLOADS ############
library(GEOquery)
library(maSigPro)
geo <- "GSE48455"
#filePaths <- getGEOSuppFiles(GEO = geo)
data <- getGEO(geo)[[1]]

# vns package for normalizaton and glog transformation

########### METADATA ############
metadata1 <- pData(data)
metadata1 = metadata1[,c("title", "geo_accession")]
write.table(metadata1, "meta-tmp.tsv", sep = "\t")
# Processed with python script: bauer_2015_GSE48455.py
metadata = read.csv("metadata.tsv", sep = "\t", header = T)
metadata[1:3,3] = c(0,0,0) # untreated are taken as starting point
########## DESIGN MATRIX ############
tags = vector()
counter = 1
#  create tags for GSM (sample IDs)
for(i in 1:nrow(metadata)){
  row = metadata[i, 3:6]
  if(row[1] == 0){tag = paste0("T",0,"_UNTR_","R",row[4],"-", counter)}
  else{
    if(row[2] == 0){tag = paste0("T",row[1],"_BLEO_","R",row[4],"-", counter)}
    if(row[2] == 1){tag = paste0("T",row[1],"_CTRL_","R",row[4],"-", counter)}
  }
  counter = counter + 1
  tags = c(tags, tag)
}
print(tags)
rownames(metadata) <- tags
#  create names file
names_file = data.frame(tags = rownames(metadata), geo = metadata$Geo)
write.table(names_file, "names.tsv", sep = "\t", row.names = F, col.names = F)
edes <- metadata[,3:6]
edes <- edes[,c(1,4,2,3)]

# make design matrix
design2 <- make.design.matrix(edesign = edes, degree = 2) # quadratic regression model
design3 <- make.design.matrix(edesign = edes, degree = 3) # cubic regression model
design4 <- make.design.matrix(edesign = edes, degree = 4) # cubic regression model

######### EXPRESSION MATRIX PREP ###########
gene_exp <- exprs(data)
gene_exp1 <- gene_exp
#  rename expression matrix from names file
for(i in c(1:length(colnames(gene_exp)))){
  for(j in c(1:length(names_file[,"tags"]))){
    expname = colnames(gene_exp)[i]
    mygeo = as.character(names_file[j,"geo"])
    if(expname == mygeo){
      colnames(gene_exp1)[i] = as.character(names_file[j,"tags"])
    }
  }
}

# check if renaming worked
counter = 1
for(i in colnames(gene_exp)){
  print(names_file[which(names_file$tags == i),2] == colnames(gene_exp1)[counter])
  counter = counter + 1
}

##### REGRESSION MODEL #####

fit2 <- p.vector(gene_exp, design2, Q = 0.05, MT.adjust = "BH", min.obs = 20) # quadratic
fit3 <- p.vector(gene_exp, design3, Q = 0.05, MT.adjust = "BH", min.obs = 20) # cubic
fit4 <- p.vector(gene_exp, design4, Q = 0.05, MT.adjust = "BH", min.obs = 20) # fourth degree

##### FIND SIGNIFICANT GENES #####

# step-wise regression
tstep2 <- T.fit(fit2, step.method = "backward", alfa = 0.05) # quadratic
tstep3 <- T.fit(fit3, step.method = "backward", alfa = 0.05) # cubic
tstep4 <- T.fit(fit4, step.method = "backward", alfa = 0.05) # cubic

# check how well fit the model is - using R-squared
mean(tstep2$sol[,"R-squared"]>0.5) #  6.9% for degree=2
mean(tstep3$sol[,"R-squared"]>0.5) # 12% for degree=3
mean(tstep4$sol[,"R-squared"]>0.5) # 21% for degree=4

sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
sigs$sig.genes$bleomycinvscontrol[1]
random_probe = expression[rownames(expression)=="A_42_P467630",]
PlotGroups(random_probe, edesign = edesign, 
           show.fit = T, dis = design$dis, groups.vector = design$groups.vector)
see.genes(sigs$sig.genes$bleomycinvscontrol, main = "bleomycinvscontrol", show.fit = T,
          dis =design$dis, cluster.method="kmeans" ,cluster.data = 1, k = 9)
