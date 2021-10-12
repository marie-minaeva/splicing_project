library(UniprotR)

load('myEnvironment.RData')
args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])

data = read.table("non_colocalizing_sQTLs.tsv", header=T)
gene_name = data$gene_id[j]
gene_name
uni_name = annotLookup[annotLookup$ensembl_gene_id == gene_name,]$uniprot_gn_id
GetSequences(uni_name[1])$Sequence
write(GetSequences(uni_name[1])$Sequence, "ref_seq.txt")

write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
              annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1]), 
      "gene_metadata.txt")
