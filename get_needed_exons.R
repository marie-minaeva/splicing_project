setwd("~/splicing_project/")

args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])

data = read.table("significant_sQTL_GWAS_coloc_events.tsv", header=T)
# Gene name extraction
gene_name = unlist(strsplit(data$phenotype_id[j],split = "_"))[1]
gene_name
bash = c("bash retrieve_nucleotide_sequence_from_genename.sh -g", gene_name)
bash =  paste(bash, collapse=" ")
print(bash)
## all gene exons extraction
system(bash)