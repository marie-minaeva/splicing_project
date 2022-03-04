setwd("~/splicing_project/")

library(stringr)
library(diptest)

args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])
data = read.csv("Data/combined_sQTL_data.csv")
#filt = readLines("Data/true_coloc.txt")
#inp = read.table("Data/significant_sQTL_GWAS_coloc_events.tsv", header=T)
#data = data[inp$phenotype_id %in% filt, ]
poss = data$ALIGN.COORDS[j]

# data$MIN_pLLDT[j] = NA
# data$Q1_pLLDT[j] = NA
# data$Q2_pLLDT[j] = NA
# data$Q3_pLLDT[j] = NA
# data$MAX_pLLDT[j] = NA
# data$DIP[j] = NA
# data$DIP_P[j] = NA

start = as.numeric(unlist(str_split(poss, pattern = '-'))[1])
start
stop = as.numeric(unlist(str_split(poss, pattern = '-'))[2])
file = paste("UP000005640_9606_HUMAN/AF-", data$ALPHAFOLD.NAME[j], "-F1-model_v1.cif", sep = '')
#file = paste("UP000005640_9606_HUMAN/AF-", "Q13530", "-F1-model_v1.cif", sep = '')

pLDDT = data.frame()
flag=FALSE
if (file.exists(file) & !is.na(start) & !is.na(stop)){
        input = readLines(file)
        k=length(input)
        for (i in 1:length(input)){
        
                if (unlist(str_split(input[i], pattern = " "))[1] == "A" & i>k & flag == TRUE){
                        stri = unlist(str_split(input[i], pattern = ' '))
                        stri = stri[stri != ""]
                        pLDDT = rbind(pLDDT, stri)
                } else {
                        flag=FALSE
                }
                if (input[i] == "_ma_qa_metric_local.ordinal_id"){
                        k = i
                        flag = TRUE
                }
        }
        pLDDT = pLDDT[,2:ncol(pLDDT)]
        for (i in 2:ncol(pLDDT)){
                pLDDT[, i] = as.numeric(pLDDT[, i]) 
        }
        colnames(pLDDT) = c("AA", "NUM", "SMT", "pLDDT", "SMTH", "SOMTH")
        print(quantile(pLDDT$pLDDT))
        #write(pLDDT$pLDDT, paste("Data/pLDDTs/", "Q13530", "_all_pLDDT.txt", sep=""))
        write(pLDDT$pLDDT[start:stop], paste("Data/pLDDTs_sQTL/", data$ALPHAFOLD.NAME[j], "_all_pLDDT.txt", sep=""))
        # data$MIN_pLLDT[j] = quantile(pLDDT$pLDDT[start:stop])[1]
        # data$Q1_pLLDT[j] = quantile(pLDDT$pLDDT[start:stop])[2]
        # data$Q2_pLLDT[j] = quantile(pLDDT$pLDDT[start:stop])[3]
        # data$Q3_pLLDT[j] = quantile(pLDDT$pLDDT[start:stop])[4]
        # data$MAX_pLLDT[j] = quantile(pLDDT$pLDDT[start:stop])[5]
        # data$DIP[j] = dip.test(pLDDT$pLDDT[start:stop])$statistic
        # data$DIP_P[j] = dip.test(pLDDT$pLDDT[start:stop])$p.value
}
#write.csv(data, "Data/output_data_non_coloc.csv", row.names = F)



