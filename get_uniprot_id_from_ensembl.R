library(UniprotR)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])

data = read.table("significant_sQTL_GWAS_coloc_events.tsv", header=T)
gene_name = data$GeneID[j]
gene_name
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
        mart = mart,
        attributes = c('ensembl_gene_id',
                'uniprotswissprot',
                'description',
                'external_gene_name'), uniqueRows=TRUE)

uni_name = annotLookup[annotLookup$ensembl_gene_id == gene_name,]$uniprotswissprot
write(GetSequences(uni_name[1])$Sequence, "ref_seq.txt")

write(c(annotLookup[annotLookup$ensembl_gene_id == gene_name,]$external_gene_name[1],
              annotLookup[annotLookup$ensembl_gene_id == gene_name,]$description[1]), 
      "gene_metadata.txt")
