setwd("~/splicing_project/")
data = read.csv("Data/output_data_coloc_try.csv")

data2 = read.csv("Data/output_data_non_coloc_try.csv")
data3 = read.csv("Data/output_data_non_sQTL_try.csv")

data = unique.data.frame(data)
data2 = unique.data.frame(data2)
nrow(data)
nrow(data2)
library(ggplot2)
library(stats)

flag_outliers <- function(x){
        quants <- quantile(x, na.rm=TRUE)[c(2,4)]
        iqr <- IQR(x, na.rm = T)
        return((x > quants[2] + 1.5 * iqr) | (x < quants[1] - 1.5 * iqr))
}

data = data[!flag_outliers(data$LENGTH),]
data2 = data2[!flag_outliers(data2$LENGTH),]
data3 = data3[!flag_outliers(data3$LENGTH),]

#View(data[!flag_outliers(data$LENGTH) & data$LENGTH!=0 ,])
ggplot() + geom_density(data=data[data$LENGTH!=0,], aes(LENGTH, fill="coloc"), alpha = 0.2) + geom_density(data=data2[data2$LENGTH!=0,], aes(LENGTH, fill="non coloc"), alpha = 0.2) #+ geom_density(data=data3, aes(LENGTH, fill="non sQTL"), alpha = 0.6)
ggplot() + geom_density(data=data[!flag_outliers(data$HYDROPATHICITY),], aes(HYDROPATHICITY, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$HYDROPATHICITY),], aes(HYDROPATHICITY, fill="non coloc"), alpha = 0.2) #+ geom_density(data=data3, aes(HYDROPATHICITY, fill="non sQTL"), alpha = 0.6)
ggplot() + geom_density(data=data[!flag_outliers(data$ASN../data$LENGTH *100),], aes(ASN../LENGTH *100, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$ASN../data2$LENGTH *100),], aes(ASN../LENGTH * 100, fill="non coloc"), alpha = 0.2)# + geom_density(data=data3, aes(CYS.., fill="non sQTL"), alpha = 0.6)
ggplot() + geom_density(data=data[!flag_outliers(data$CYS../data$LENGTH *100),], aes(CYS../LENGTH *100, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$CYS../data2$LENGTH *100),], aes(CYS../LENGTH * 100, fill="non coloc"), alpha = 0.2)# + geom_density(data=data3, aes(CYS.., fill="non sQTL"), alpha = 0.6)
ggplot() + geom_density(data=data, aes(CYS../LENGTH *100, fill="coloc"), alpha = 0.2) + geom_density(data=data2, aes(CYS../LENGTH * 100, fill="non coloc"), alpha = 0.2)# + geom_density(data=data3, aes(CYS.., fill="non sQTL"), alpha = 0.6)

data = cbind(data, data$LENGTH %% 3)
data2 = cbind(data2, data2$LENGTH %% 3)
data3 = cbind(data3, data3$LENGTH %% 3)

ggplot() + geom_density(data=data[!flag_outliers(data$LENGTH %% 3),], aes(LENGTH %% 3, fill="coloc"), alpha=0.3) + geom_density(data=data2[!flag_outliers(data2$LENGTH %% 3),], aes(LENGTH %% 3, fill="non coloc"), alpha=0.3)# + geom_density(data=data3, aes(LENGTH %% 3, fill="non sQTL"), alpha = 0.6)



ggplot() + geom_density(data=data[!flag_outliers(data$MEAN.ACC),], aes(MEAN.ACC, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$MEAN.ACC),], aes(MEAN.ACC, fill="non coloc"), alpha = 0.2) #+ geom_density(data=data3, aes(LENGTH, fill="non sQTL"), alpha = 0.6)




ggplot() + geom_density(data=data, aes(HELIX, fill="coloc"), alpha = 0.2) + geom_density(data=data2, aes(HELIX, fill="non coloc"), alpha = 0.2) #+ geom_density(data=data3, aes(LENGTH, fill="non sQTL"), alpha = 0.6)
ggplot() + geom_density(data=data, aes(SHEET, fill="coloc"), alpha = 0.2) + geom_density(data=data2, aes(SHEET, fill="non coloc"), alpha = 0.2) #+ geom_density(data=data3, aes(LENGTH, fill="non sQTL"), alpha = 0.6)
ggplot() + geom_density(data=data, aes(TURN, fill="coloc"), alpha = 0.2) + geom_density(data=data2, aes(TURN, fill="non coloc"), alpha = 0.2) #+ geom_density(data=data3, aes(LENGTH, fill="non sQTL"), alpha = 0.6)





ks.test(data[!flag_outliers(data$HYDROPATHICITY),]$HYDROPATHICITY, data2[!flag_outliers(data2$HYDROPATHICITY),]$HYDROPATHICITY)
ks.test(data$LENGTH, data2$LENGTH)
ks.test(data[!flag_outliers(data$ASN../data$LENGTH *100),]$ASN../data$LENGTH*100, data2[!flag_outliers(data2$ASN../data2$LENGTH *100),]$ASN../data2$LENGTH*100)
ks.test(data[!flag_outliers(data$CYS../data$LENGTH *100),]$CYS../data$LENGTH*100, data2[!flag_outliers(data2$CYS../data2$LENGTH *100),]$CYS../data2$LENGTH*100)
ks.test(data[!flag_outliers(data$LENGTH %% 3),]$LENGTH %% 3, data2[!flag_outliers(data2$LENGTH %% 3),]$LENGTH %% 3.)
ks.test(data[!flag_outliers(data$MEAN.ACC),]$MEAN.ACC, data2[!flag_outliers(data2$MEAN.ACC),]$MEAN.ACC)





nrow(data)
nrow(data2)

symet = data.frame("coloc" = c(length(data[data$LENGTH %% 3 == 0,]$LENGTH), length(data[data$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_coloc" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
test

symet_not_found = data.frame("coloc" = c(length(data[data$LENGTH %% 3 == 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data[data$LENGTH %% 3 != 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                   "non_coloc" = c(length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data2[data2$LENGTH %% 3 != 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)
symet_not_found
test = fisher.test(symet_not_found)
test


hel = data.frame("coloc" = c(length(data[data$HELIX == T,]$HELIX), length(data[data$HELIX != T,]$HELIX)),
                   "non_coloc" = c(length(data2[data2$HELIX == T,]$HELIX), length(data2[data2$HELIX != T,]$HELIX)),
                   row.names = c("helix", "no_helix"),
                   stringsAsFactors = FALSE)

hel
test = fisher.test(hel)
test

she = data.frame("coloc" = c(length(data[data$SHEET == T,]$SHEET), length(data[data$SHEET != T,]$SHEET)),
                 "non_coloc" = c(length(data2[data2$SHEET == T,]$SHEET), length(data2[data2$SHEET != T,]$SHEET)),
                 row.names = c("sheet", "no_sheet"),
                 stringsAsFactors = FALSE)

she
test = fisher.test(she)
test

turn = data.frame("coloc" = c(length(data[data$TURN == T,]$TURN), length(data[data$TURN != T,]$TURN)),
                 "non_coloc" = c(length(data2[data2$TURN == T,]$TURN), length(data2[data2$TURN != T,]$TURN)),
                 row.names = c("turn", "no_turn"),
                 stringsAsFactors = FALSE)

turn
test = fisher.test(turn)
test


ggplot() + geom_bar(data=data, aes(DOMAIN, fill="coloc"), stat = "count", alpha=0.3) + geom_bar(data=data2, aes(DOMAIN, fill="non coloc"), stat = "count", alpha=0.3)# + geom_density(data=data3, aes(LENGTH %% 3, fill="non sQTL"), alpha = 0.6)
View(data)
table(data$DOMAIN)[table(data$DOMAIN) >= 2]
table(data2$DOMAIN)[table(data2$DOMAIN) >= 2]
