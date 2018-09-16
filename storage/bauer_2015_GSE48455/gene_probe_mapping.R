 
expression = read.csv("~/windows/tcc/storage/bauer_2015_GSE48455/gene_expression_matrix", sep = "\t",
                      stringsAsFactors = FALSE)
agilentids = as.vector(rownames(expression))

# FROM GPL AT GEO website
glp = read.csv("GPL7294-9579.txt", header = TRUE, skip = 17, sep = "\t", stringsAsFactors = FALSE)[,c(1, 7)]
mapped_genes = vector()
mapped_df = data.frame()
for(i in agilentids){
  gene_name = glp[which(glp$ID == i), 2]
  if(nchar(gene_name) == 0){
    mapped_genes = c(mapped_genes, "Nothing")
  }
  else{
    print(gene_name)
    mapped_genes = c(mapped_genes, gene_name)
  }
  
}

sum(mapped_genes == "Nothing") # 18509 unmapped probes!!
write.table(data.frame(probes = agilentids, genes = mapped_genes),
            file = "probe_mapping.tsv", sep = "\t", row.names = F)
# Try biomaRt!

library(biomaRt)
ensembl=useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
geneids = getBM(attributes = c('agilent_wholegenome_4x44k_v3', 'external_gene_name'),
                filters = 'agilent_wholegenome_4x44k_v3', values = agilentids, mart = ensembl)

dim(geneids) # only returns 9548 entries... not as good then

