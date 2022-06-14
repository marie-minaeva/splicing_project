setwd("~/splicing_project/")
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])

data = read.table("Data/top_sQTLs_top_coloc.tsv", header=T, sep='\t')
# Gene name extraction
#gene_name = unlist(strsplit(data$phenotype_id[j],split = "_"))[1]
# for non sQTL
gene_name = unlist(str_split(data$top_pid[j],pattern = "_"))[1]
gene_name
bash = c("bash retrieve_nucleotide_sequence_from_genename.sh -g", gene_name)
bash =  paste(bash, collapse=" ")
print(bash)
## all gene exons extraction
system(bash)
