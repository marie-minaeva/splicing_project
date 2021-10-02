setwd("~/splicing_project/")

gene_meta = readLines("gene_metadata.txt")

align_meta = readLines("best_align_meta.txt")

protein = readLines("best_aligned_protein.txt")

protein_len = nchar(protein)

exon = 0
prot_meta = readLines("protein_meta.txt")


if (file.exists("output_data.csv")){
        input = read.csv("output_data.csv")
} else {
        input = data.frame()
}

line = c(gene_meta, exon, align_meta, protein, protein_len, prot_meta)
input = rbind(input, line)
colnames(input) = c("GENE NAME", "PROTEIN", "EXON", "READ FRAME", "STRAND", "DIRECTION", "ALIGNED SEQ", "LENGT", "ASN #", "CYS #", "HYDROPATHICITY")
write.csv(input, "output_data.csv", row.names = F)
View(input)
