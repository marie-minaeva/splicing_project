setwd("~/splicing_project/")

if (file.exists("Data/output_non_sQTL_full.csv")){
                 input1 = read.csv("Data/output_non_sQTL_full.csv")
         } else {
                         input1 = data.frame()
                 }
View(input1)

library(dplyr)
library(stringr)

pdb <- read.pdb( "UP000005640_9606_HUMAN/AF-A0A0A0MRZ7-F1-model_v1.pdb" )
pdbseq(pdb)
paste(pdbseq(pdb), collapse = "")

View(read.table("~/Downloads/uniclust_uniprot_mapping.tsv")[1:1000,])
View(listAttributes(mart))


data = read.table("Data/cross_tissue_nonsignificant_genes.tsv", header = T)
data_out = sample_n(data, 1000)

colnames(data_out) = colnames(data)
colnames(data_out)
write.table(data_out, "Data/cross_tissue_nonsignificant_genes.tsv", row.names = F, sep = '\t', col.names = T)


View(d1)

load('myEnvironment.RData')
d1 = read.table("Data/top_sQTLs_top_coloc.tsv", header = T, sep="\t")
d2 = read.csv("Data/output_top_sQTL.csv")
d2 = rbind(d2[1:2234,], c(2234, "ABCA2", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), d2[2235:(nrow(d2)-1),])
uni_name = annotLookup[annotLookup$ensembl_gene_id %in% d1$group,]$external_gene_name
uni_name = unique(uni_name)
setdiff(uni_name, d2$GENE.NAME)
annotLookup[annotLookup$external_gene_name == setdiff(uni_name, d2$GENE.NAME),]
View(d2)
View(cbind(d1, d2[2:769,]))
View(d1)
nrow(d2)
nrow(d1)
diff = setdiff(d1$gene_symbol, d2$GENE.NAME)
diff
d1 = d1[ !(d1$gene_symbol %in% diff), ]
write.csv(cbind(d1, d2[, 2:ncol(d2)]), "Data/combined_sQTL_data.csv", row.names = F)


d1 = read.table("Data/cross_tissue_nonsignificant_genes.tsv", header = T)
d1 = distinct(d1, top_pid, .keep_all= TRUE)
View(d1)
write.table(d1, "Data/non_colocalizing_sQTLs.tsv", sep='\t', row.names = F)
nrow(d1)



library(UniprotR)

write(GetSequences('A0A0A0MRZ7')$Sequence, "A0A0A0MRZ7.fa")
load('myEnvironment.RData')



prots = readLines("prots_in_AlphaFold.txt")
prots
for (i in 1:length(prots)){
        prots[i] = unlist(strsplit(prots[i], split='-'))[2]
}
write(prots, "prots_in_AlphaFold.txt")
View(prots)
for (j in 1:10){
        data = read.table("Data/non_colocalizing_sQTLs.tsv", header=T)
        gene_name = data$gene_id[j]
        gene_name
        uni_name = annotLookup[annotLookup$ensembl_gene_id == gene_name,]$uniprot_gn_id
        uni_name
        gene = intersect(uni_name, prots)
        print(gene)
        print(GetSequences(gene)$Sequence)
}


try = read.csv("Data/top_sQTLs_top_coloc.tsv", header = T, sep = '\t')
colnames(try)
table(try$has_coloc)
has_coloc = try[try$has_coloc == "At least one\ncolocalized trait", ]$top_pid
View(has_coloc)
write(has_coloc, "Data/true_coloc.txt")

new_data = read.table("Data/top_sQTLs_top_coloc.tsv", header = T, sep="\t")
View(new_data)
nrow(new_data)
new_data = new_data[new_data$PP.H4.abf > 0.5, ]
nrow(new_data)

#write.table(new_data, "Data/new_coloc_data.tsv", row.names = F, sep = '\t', col.names = T)
write(new_data$top_pid, "Data/true_coloc_new_new.txt")



data = read.csv("Data/combined_sQTL_data.csv")
data2 = data[data$has_coloc == "No colocalization", ]
nrow(data2)
write.csv(data2, "Data/output_data_non_coloc.csv", row.names=F)
