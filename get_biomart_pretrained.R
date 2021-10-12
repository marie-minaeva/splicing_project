library(biomaRt)


mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
        mart = mart,
        attributes = c('ensembl_gene_id',
                       'uniprot_gn_id',
                       'description',
                       'external_gene_name'), uniqueRows=TRUE)
save.image(file='myEnvironment.RData')