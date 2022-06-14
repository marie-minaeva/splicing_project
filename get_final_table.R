library(data.table)
library(stringr)


setwd("~/splicing_project/")

data = fread("try_blastx.txt", header=F)
print(head(data))

if (file.exists("Data/try_new_pipeline.csv")){
        input = read.csv("Data/try_new_pipeline.csv")
	colnames(input) = c("GENE ID", "PROT REFSEQ", "START", "END", "SEQ", "ALIGN LENGTH", "EVAL", "SCORE", "ASN #", "CYS #")
} else {
        input = data.frame()
}

line = data.table(data[1,], str_count(data[1, "V5"], pattern="C"), str_count(data[1, "V5"], pattern="N"))
colnames(line) = c("GENE ID", "PROT REFSEQ", "START", "END", "SEQ", "ALIGN LENGTH", "EVAL", "SCORE", "ASN #", "CYS #")


line = rbind(input, line)
write.csv(line, "Data/try_new_pipeline.csv", row.names = F)
