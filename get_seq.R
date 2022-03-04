setwd("~/splicing_project/")


args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])
data = read.csv("Data/output_non_sQTL_full.csv")

prot = data$SEQ[j]
write(prot, "temp.fa")
