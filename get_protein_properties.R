library(Peptides)
library(stringr)

setwd("~/splicing_project/")
seq = readLines("best_aligned_protein.txt")
hyd = as.character(hydrophobicity(seq))
asn_n = as.character(str_count(seq, "N"))
cys_n = as.character(str_count(seq, "C"))
write(c(asn_n, cys_n, hyd), "protein_meta.txt")
