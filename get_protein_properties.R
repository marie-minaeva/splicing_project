library(Peptides)
library(stringr)

setwd("~/splicing_project/")
seq = readLines("protein_seq.txt")
if (seq == "NOT FOUND"){
        hyd = ''
        asn_n = ''
        cys_n = ''
} else {
        hyd = as.character(hydrophobicity(seq))
        asn_n = as.character(str_count(seq, "N"))
        cys_n = as.character(str_count(seq, "C"))
}
write(c(asn_n, cys_n, hyd), "protein_meta.txt")
