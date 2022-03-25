setwd("~/splicing_project/")

library(stringr)


out_data = read.csv("Data/output_non_sQTL_full.csv")
out_data$mean_01_psi = NA


files = list.files("Data/median_psi/")
files
for (file in files){
        data = read.table(paste0("Data/median_psi/", file), header = T, sep = '\t')
        tis = unlist(str_split(file, pattern = '_median'))[1]
        print(tis)
        if (file == "README"){
                next
        }
        data = data[data$exon_id %in% out_data[out_data$tiss == tis, ]$top_pid, ]
        out_data[out_data$tiss == tis, ]$mean_01_psi = data$median_psi
        
}
View(out_data)
write.csv(out_data, "Data/output_non_sQTL_full.csv", row.names = F)