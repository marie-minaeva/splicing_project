library(stringr)


args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])
data = read.csv("Data/output_non_sQTL_full.csv")
#View(data)
domains = readLines("domain.txt")
domains = unlist(str_split(domains, " "))
domains = domains[which(startsWith(domains, "alihmmname"))]
domains = unlist(str_split(domains, "\""))[2]
print(domains)
data$DOMAIN[j] = domains
write.csv(data, "Data/output_non_sQTL_full.csv", row.names = F)

