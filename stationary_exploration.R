setwd("~/splicing_project/")
data = read.csv("Data/output_data_non_coloc.csv")
stationary_non = data[data$KPSS_P >= 0.05 & data$ADF_P <= 0.05, ]
stationary_non = stationary_non[!is.na(stationary_non$GENE.NAME), ]


data1 = read.csv("Data/output_data_coloc.csv")
stationary = data1[data1$KPSS_P >= 0.05 & data1$ADF_P <= 0.05, ]
stationary = stationary[!is.na(stationary$GENE.NAME), ]

station = data.frame("sQTL" = c(length(stationary$GENE.NAME), length(data1$GENE.NAME) - length(stationary$GENE.NAME)),
                             "non_sQTL" = c(length(stationary_non$GENE.NAME), length(data$GENE.NAME) - length(stationary_non$GENE.NAME)),
                             #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             row.names = c("stat", "non_stat"),
                             stringsAsFactors = FALSE)
station
test = fisher.test(station)
test