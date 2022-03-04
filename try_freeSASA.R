library(stringr)

args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])
data = read.csv("Data/output_non_sQTL_full.csv")
poss = data$ALIGN.COORDS[j]
data$MIN[j] = NA
data$Q1[j] = NA
data$Q2[j] = NA
data$Q3[j] = NA
data$MAX[j] = NA

start = as.numeric(unlist(str_split(poss, pattern = '-'))[1])
start
stop = as.numeric(unlist(str_split(poss, pattern = '-'))[2])
file = "try_freesas.txt"
if (file.exists(file) & !is.na(start) & !is.na(stop)){
        try = readLines(file)
        parsed = data.frame()
        for (i in 2:length(try)){
                parser = unlist(str_split(try[i], pattern = " "))
                parser = parser[parser != ""]
                if (length(parser) == 14){
                        parser[3] = paste(parser[3], parser[4], "")
                        parser = parser[-4]
                }
                parsed = rbind(parsed, parser)
        }
        
        parsed = parsed[, 2:ncol(parsed)]
        colnames(parsed) = c('RES', 'NUM', 'All-atoms_ABS', 'All-atoms_REL', 'Total-Side_ABS', 'Total-Side_REL', 'Main-Chain_ABS', 'Main-Chain_REL', 'Non-polar_ABS', 'Non-polar_REL', 'All-polar_ABS', 'All-polar_REL')
        parsed$`All-atoms_ABS` = as.numeric(parsed$`All-atoms_ABS`)
        for (i in 3:ncol(parsed)){
                parsed[, i] = as.numeric(parsed[, i]) 
        }
        data$MIN[j] = quantile(parsed$`All-atoms_REL`[start:stop])[1]
        data$Q1[j] = quantile(parsed$`All-atoms_REL`[start:stop])[2]
        data$Q2[j] = quantile(parsed$`All-atoms_REL`[start:stop])[3]
        data$Q3[j] = quantile(parsed$`All-atoms_REL`[start:stop])[4]
        data$MAX[j] = quantile(parsed$`All-atoms_REL`[start:stop])[5]
}
write.csv(data, "Data/output_non_sQTL_full.csv", row.names = F)
