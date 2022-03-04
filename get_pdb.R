args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])
data = read.csv("Data/output_non_sQTL_full.csv")
name = data$ALPHAFOLD.NAME[j]
name
file = paste("UP000005640_9606_HUMAN/AF-", name, "-F1-model_v1.pdb", sep='')
write(file, file = "name.txt")
