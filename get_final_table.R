setwd("~/splicing_project/")

gene_meta = readLines("gene_metadata.txt")

align_meta = readLines("best_align_meta.txt")

protein = readLines("best_aligned_protein.txt")
prot_seq = readLines("protein_seq.txt")

if (protein == "NOT IN ALPHAFOLD" || protein == "EXON NOT FOUND"){
        protein_len_ali = ''
} else {
        protein_len_ali = nchar(protein)
}
prot_len =  nchar(prot_seq)
exon = readLines("exon_meta.txt")
prot_meta = readLines("protein_meta.txt")
align_coords = readLines("align_coords.txt")

if (file.exists("Data/output_data_non_sQTL_try.csv")){
        input = read.csv("Data/output_data_non_sQTL_try.csv")
} else {
        input = data.frame()
}

line = c(gene_meta, exon, align_meta, prot_seq, protein, prot_len, protein_len_ali, align_coords, prot_meta)
input = rbind(input, line)
colnames(input) = c("GENE NAME", "PROTEIN", "ALPHAFOLD NAME", "EXON", "READ FRAME", "STRAND", "DIRECTION", "SEQ", "ALIGNED SEQ", "LENGTH", "LENGTH ALIGN", "ALIGN COORDS", "ASN #", "CYS #", "HYDROPATHICITY")
write.csv(input, "Data/output_data_non_sQTL_try.csv", row.names = F)
View(input)
