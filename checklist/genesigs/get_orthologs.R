library(tidyverse, lib ="~/R/x86_64-pc-linux-gnu-library/3.4/")
library(stringr, lib ="~/R/x86_64-pc-linux-gnu-library/3.4/")

# load gene symbols
dat <- read.csv("/home/giulianonetto/windows/tcc/checklist/genesigs/papers/Xue, 2014/modules.csv", 
                        sep = "\t", header = T, stringsAsFactors = F)
# Rename columns
names(dat) <- sapply(names(dat), function(x) paste0("Module", str_extract(x, "[0-9]+")))

# get unique gene symbols from Xue, 2014
valid_elem=c()
for(i in 1:ncol(dat)){
  valid_elem=c(valid_elem, dat[which(dat[,i] != ""), i])
}
genes_list <- data.frame('GeneSymbols' = unique(valid_elem), stringsAsFactors = F)
rm(valid_elem)

# Convert human to mouse gene names

convert_human_gene_list <- function(x){

  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  mouse_matches = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = x , mart = human, 
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  no_matches = setdiff(x, mouse_matches[,1]) 
  
  x = data.frame(HGNC.symbol = x)
  df <- merge(x, mouse_matches, by = "HGNC.symbol")

    return(list(df, data.frame(no_matches)))
}

genes <- convert_human_gene_list(genes_list$GeneSymbols)
mymatches <- genes[[1]]
nomathes <- genes[[2]]
nrow(nomathes) # over 2 thousand genes with no matches
nrow(mymatches)


# write multiple tables for manual bioGPS queries

# write.table(nomathes[1:500,1], "nomatches1.500.txt", quote = F, row.names = F, col.names = F)
# write.table(nomathes[501:1000,1], "nomatches501.1000.txt", quote = F, row.names = F, col.names = F)
# write.table(nomathes[1001:1500,1], "nomatches1001.1500.txt", quote = F, row.names = F, col.names = F)
# write.table(nomathes[1501:2000,1], "nomatches1501.2000.txt", quote = F, row.names = F, col.names = F)
# write.table(nomathes[2000:nrow(nomathes),1], "nomatches2000.end.txt", quote = F, row.names = F, col.names = F)

"
first we made five queries at http://biogps.org/#goto=welcome for the total of
2331 gene HGCN symbols from which we got no matches from biomaRt. Then we
got the respective gene symbols for human, mouse and rat into .csv tables
(one for each query). Then we made five files with all 'matchless' 
queries from bioGPS.
"

"
it appears that 431 human genes had no mouse orthology mathes through biomaRt
neither any matcc at bioGPS (not even for human genes).
Further validation needed?
"

biogps <- read.csv("nomatchesfrombiomart/biogps.allmatches.csv", 
                   sep = "\t", header = T, stringsAsFactors = F)
biogps <- biogps %>% filter(SPECIES != "SPECIES") # take out repeated headers
table(biogps$SPECIES) 
names(biogps) <- c("Species", "HGNC.symbol", "MGI.symbol")
"
From biogps: 1917 human genes, 924 mouse orthologs and 843 rat orthologs.
"

# Try getting MGI symbols for bioGPS results

biogps.humans.biomart = biogps %>% filter(Species == "human") %>% 
  dplyr::select(MGI.symbol) %>% unlist() %>% convert_human_gene_list()
biogps.humans.biomart.MGI <- biogps.humans.biomart[[1]]
nrow(biogps.humans.biomart.MGI) # 1232 mouse matches! ~ 300 more...

matches_all = biogps %>% filter(Species == "mouse") %>% 
              dplyr::select(HGNC.symbol, MGI.symbol) %>%
              rbind(mymatches, biogps.humans.biomart.MGI) %>% unique()
matches_all%>%nrow()

matches_duplicated = matches_all[duplicated(matches_all$HGNC.symbol),]
matches_unique = matches_all[!duplicated(matches_all$HGNC.symbol),] 
# 9117 unique HGNC entries

matches_duplicated_mgi = matches_all[duplicated(matches_all$MGI.symbol),]
matches_unique_mgi = matches_all[!duplicated(matches_all$MGI.symbol),]
# 9073 unique MGI entries

no_mouse_orthologs = data.frame(HGNC.symbol = biogps.humans.biomart[[2]][[1]],
                                MGI.symbol = rep(NA, 
                                    length(biogps.humans.biomart[[2]][[1]])))

final.table = data.frame(rbind(matches_unique, no_mouse_orthologs))

final.table = final.table[!duplicated(final.table$HGNC.symbol),]
write.table(final.table, 
            "/home/giulianonetto/windows/tcc/checklist/genesigs/human_mouse.map",
            sep = "\t", quote = F, row.names = F)
