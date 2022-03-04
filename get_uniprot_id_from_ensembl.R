library(UniprotR)
setwd("~/splicing_project/")
## SIGNIFICANT GENES

# load('myEnvironment.RData')
# args <- commandArgs(trailingOnly = TRUE)
# j = as.numeric(args[1])
# 
# prots = readLines("prots_in_AlphaFold.txt")
# data = read.table("Data/significant_sQTL_GWAS_coloc_events.tsv", header=T)
# gene_name = data$GeneID[j]
# gene_name
# uni_name = annotLookup[annotLookup$ensembl_gene_id == gene_name,]$uniprot_gn_id
# gene = intersect(uni_name, prots)
# GetSequences(gene)$Sequence
# 
# write(GetSequences(gene)$Sequence, "ref_seq.txt")
# print(gene)
# if (identical(gene, character(0)) ){
#         print("GEEEEEEENEEEEEEE1111111111")
#         write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
#                 annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1], ""), "gene_metadata.txt")
# }else{
#         print("GEEEEEEENEEEEEEE")
#         write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
#                 annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1], gene), "gene_metadata.txt")
# }

## NON SIGNIFICANT GENES

# load('myEnvironment.RData')
# args <- commandArgs(trailingOnly = TRUE)
# j = as.numeric(args[1])
# 
# prots = readLines("prots_in_AlphaFold.txt")
# data = read.table("Data/non_colocalizing_sQTLs.tsv", header=T)
# gene_name = data$gene_id[j]
# uni_name = annotLookup[annotLookup$ensembl_gene_id == gene_name,]$uniprot_gn_id
# gene = intersect(uni_name, prots)
# 
# 
# write(GetSequences(gene)$Sequence, "ref_seq.txt")
# print(gene)
# if (identical(gene, character(0)) ){
#         print("GEEEEEEENEEEEEEE1111111111")
#         write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
#                 annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1], ""), "gene_metadata.txt")
# }else{
#         print("GEEEEEEENEEEEEEE")
#         write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
#                 annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1], gene), "gene_metadata.txt")
# }





## NON sQTL GENES

load('myEnvironment.RData')
args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])

prots = readLines("prots_in_AlphaFold.txt")
data = read.table("Data/cross_tissue_nonsignificant_genes.tsv", header=T)
gene_name = data$group[j]
uni_name = annotLookup[annotLookup$ensembl_gene_id == gene_name,]$uniprot_gn_id
gene = intersect(uni_name, prots)


write(GetSequences(gene)$Sequence, "ref_seq.txt")
print(gene)
if (identical(gene, character(0)) ){
        print("GEEEEEEENEEEEEEE1111111111")
        write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
                annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1], ""), "gene_metadata.txt")
}else{
        print("GEEEEEEENEEEEEEE")
        write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
                annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1], gene), "gene_metadata.txt")
}

#All top sQTLs

# load('myEnvironment.RData')
# args <- commandArgs(trailingOnly = TRUE)
# j = as.numeric(args[1])
# 
# prots = readLines("prots_in_AlphaFold.txt")
# data = read.table("Data/top_sQTLs_top_coloc.tsv", header=T, sep="\t")
# gene_name = data$group[j]
# uni_name = annotLookup[annotLookup$ensembl_gene_id == gene_name,]$uniprot_gn_id
# gene = intersect(uni_name, prots)
# 
# 
# write(GetSequences(gene)$Sequence, "ref_seq.txt")
# print(gene)
# if (identical(gene, character(0)) ){
#         print("GEEEEEEENEEEEEEE1111111111")
#         write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
#                 annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1], ""), "gene_metadata.txt")
# }else{
#         print("GEEEEEEENEEEEEEE")
#         write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
#                 annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1], gene), "gene_metadata.txt")
# }
