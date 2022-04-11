setwd("~/splicing_project/")
library(ggplot2)
library(stats)
library(dplyr)
library(mltools)
library(foreach)
library(doParallel)
library(doSNOW)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

flag_outliers <- function(x){
  quants <- quantile(x, na.rm=TRUE)[c(2,4)]
  iqr <- IQR(x, na.rm = T)
  return((x > quants[2] + 1.5 * iqr) | (x < quants[1] - 1.5 * iqr))
}

meanDiff = function(data, func){
  data = func(data)
  m1 <- mean(data[data$V2=='sQTL', "totalData"])
  m2 <- mean(data[data$V2=='non_sQTL', "totalData"])
  return(m1-m2)
}

inv_norm <- function(x) qnorm((rank(x, na.last='keep') - 0.5)/sum(!is.na(x)))

my_bootstrap = function(x, R, func, func1){
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
    dist<- foreach(
      i = 1:R, 
      .combine = c,
      .packages = 'foreach'
    ) %dopar% (func(x, func1))
  parallel::stopCluster(cl)
  return(dist)
}

select_samples = function(x){
  # d <- foreach(
  #   j = 1:length(x[, 1]), 
  #   .combine = rbind,
  #   .packages = 'foreach'
  # ) %dopar% (
  #   x[sample(1:length(x[,1]), 1), ]
  # )
  d = x[x[,1] %in% sample(x[,1], replace = TRUE),]
  colnames(d) = c("totalData", "V2")
  print(nrow(d))
  return(d)
}


ks_stat = function(x, y){
  x = x[!is.na(x)]
  y = y[!is.na(y)]
  pval =  ks.test(x, y)$p.value
  totalData = c(x, y)
  totalData = cbind(totalData, c(replicate(length(x), "sQTL"), replicate(length(y), "non_sQTL")))
  totalData[,1] = inv_norm(as.numeric(totalData[,1]))
  totalData = data.frame(totalData)
  totalData = totalData[!is.na(totalData[,1]), ]
  totalData[,1] = as.numeric(totalData[,1])
  print(mean(totalData[totalData$V2 == "sQTL", 'totalData']) - mean(totalData[totalData$V2 == "non_sQTL", 'totalData']))
  
  totalBoot = my_bootstrap(totalData, R=1000, meanDiff, select_samples)
  totalBootCI = quantile(totalBoot, probs=c(0.05, 0.95))
  conf_low = totalBootCI[1]
  conf_top = totalBootCI[2]
  stat = mean(totalData[totalData$V2 == "sQTL", 'totalData']) - mean(totalData[totalData$V2 == "non_sQTL", 'totalData'])
  print(c(pval, stat, conf_low, conf_top))
  return(c(pval, stat, conf_low, conf_top))
}





data = read.csv("Data/output_top_sQTL.csv")
data2 = read.csv("Data/output_non_sQTL_full.csv")

ncol(data)
ncol(data2)
data = data[,2:ncol(data)]
data = unique.data.frame(data)
data2 = unique.data.frame(data2)
nrow(data)
nrow(data2)



data = data[!flag_outliers(data$LENGTH) & data$LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$LENGTH) & data2$LENGTH >= 15,]
nrow(data)
nrow(data2)
View(data)
dev.off()

#View(data)

# 
# ggplot() + geom_density(data=data[data$LENGTH!=0,], aes(LENGTH, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[data2$LENGTH!=0,], aes(LENGTH, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$HYDROPATHICITY),], aes(HYDROPATHICITY, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$HYDROPATHICITY),], aes(HYDROPATHICITY, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$ASN../data$LENGTH *100),], aes(ASN../LENGTH *100, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$ASN../data2$LENGTH *100),], aes(ASN../LENGTH * 100, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$CYS../data$LENGTH *100),], aes(CYS../LENGTH *100, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$CYS../data2$LENGTH *100),], aes(CYS../LENGTH * 100, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data, aes(CYS../LENGTH *100, fill="sQTL"), alpha = 0.2) + geom_density(data=data2, aes(CYS../LENGTH * 100, fill="non sQTL"), alpha = 0.2)
# data = cbind(data, data$LENGTH %% 3)
# data2 = cbind(data2, data2$LENGTH %% 3)
# data3 = cbind(data3, data3$LENGTH %% 3)
# 
# ggplot() + geom_density(data=data[!flag_outliers(data$LENGTH %% 3),], aes(LENGTH %% 3, fill="sQTL"), alpha=0.3) + geom_density(data=data2[!flag_outliers(data2$LENGTH %% 3),], aes(LENGTH %% 3, fill="non sQTL"), alpha=0.3)
# 
# 
# ggplot() + geom_density(data=data[!flag_outliers(data$MIN),], aes(MIN, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$MIN),], aes(MIN, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$Q1),], aes(Q1, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q1),], aes(Q1, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$Q2),], aes(Q2, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q2),], aes(Q2, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$Q3),], aes(Q3, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q3),], aes(Q3, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$MAX),], aes(MAX, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$MAX),], aes(MAX, fill="non sQTL"), alpha = 0.2)
# 
# 
# ggplot() + geom_density(data=data[!flag_outliers(data$MIN_pLLDT),], aes(MIN_pLLDT, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$MIN_pLLDT),], aes(MIN_pLLDT, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$Q1_pLLDT),], aes(Q1_pLLDT, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q1_pLLDT),], aes(Q1_pLLDT, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$Q2_pLLDT),], aes(Q2_pLLDT, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q2_pLLDT),], aes(Q2_pLLDT, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$Q3_pLLDT),], aes(Q3_pLLDT, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q3_pLLDT),], aes(Q3_pLLDT, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data[!flag_outliers(data$MAX_pLLDT),], aes(MAX_pLLDT, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$MAX_pLLDT),], aes(MAX_pLLDT, fill="non sQTL"), alpha = 0.2)
# 
# 
# 
# 
# ggplot() + geom_density(data=data, aes(HELIX, fill="sQTL"), alpha = 0.2) + geom_density(data=data2, aes(HELIX, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data, aes(SHEET, fill="sQTL"), alpha = 0.2) + geom_density(data=data2, aes(SHEET, fill="non sQTL"), alpha = 0.2)
# ggplot() + geom_density(data=data, aes(TURN, fill="sQTL"), alpha = 0.2) + geom_density(data=data2, aes(TURN, fill="non sQTL"), alpha = 0.2)

data_unstruc = data[data$Q2 > 40.0 & data$Q2_pLLDT < 50.0, ]
data_unstruc = data_unstruc[!is.na(data_unstruc$GENE.NAME),] 
nrow(data_unstruc)
data2_unstruc = data2[data2$Q2 > 40.0 & data2$Q2_pLLDT < 50.0, ]
data2_unstruc = data2_unstruc[!is.na(data2_unstruc$GENE.NAME),] 
nrow(data2_unstruc)






data_struc = data[data$Q2 < 25.0 & data$Q2_pLLDT >= 90.0, ]
data_struc = data_struc[!is.na(data_struc$GENE.NAME),] 
data2_struc = data2[data2$Q2 < 25.0 & data2$Q2_pLLDT >= 90.0, ]
data2_struc = data2_struc[!is.na(data2_struc$GENE.NAME),] 


data_cons = data[data$Q2_pLLDT >= 90.0, ]
data_cons = data_cons[!is.na(data_cons$GENE.NAME),]
data2_cons = data2[data2$Q2_pLLDT >= 90.0, ]
data2_cons = data2_cons[!is.na(data2_cons$GENE.NAME),] 


data_loop = data[data$TURN == "True" & data$HELIX == "False" & data$SHEET == "False", ]
data_loop = data_loop[!is.na(data_loop$GENE.NAME),] # 59/571 = 10,3 % 37 signal

data2_loop = data2[data2$TURN == "True" & data2$HELIX == "False" & data2$SHEET == "False", ]
data2_loop = data2_loop[!is.na(data2_loop$GENE.NAME),] # 54/686 = 7,9 % 32 signal



to_draw = data.frame()
round(nrow(data_unstruc) / nrow(data) * 100)
pa = replicate(round(nrow(data_unstruc) / nrow(data) * 100), "sQTL")
pa = c(pa, replicate(round(nrow(data2_unstruc) / nrow(data2) * 100), "non_sQTL"))
pa = c(pa, replicate(round(nrow(data_struc) / nrow(data) * 100), "sQTL"))
pa = c(pa, replicate(round(nrow(data2_struc) / nrow(data2) * 100), "non_sQTL"))
pa = c(pa, replicate(round(nrow(data_cons) / nrow(data) * 100), "sQTL"))
pa = c(pa, replicate(round(nrow(data2_cons) / nrow(data2) * 100), "non_sQTL"))
pa = c(pa, replicate(round(nrow(data_loop) / nrow(data) * 100), "sQTL"))
pa = c(pa, replicate(round(nrow(data2_loop) / nrow(data2) * 100), "non_sQTL"))
trait = replicate(round(nrow(data_unstruc) / nrow(data)*100) + round(nrow(data2_unstruc) / nrow(data2)*100), "Unstructured")
trait = c(trait, replicate(round(nrow(data_struc) / nrow(data)*100) + round(nrow(data2_struc) / nrow(data2)*100), "Structured"))
trait = c(trait, replicate(round(nrow(data_cons) / nrow(data)*100) + round(nrow(data2_cons) / nrow(data2)*100), "Conserved"))
trait = c(trait, replicate(round(nrow(data_loop) / nrow(data)*100) + round(nrow(data2_loop) / nrow(data2)*100), "Loop"))

to_draw = cbind(trait, pa)

colnames(to_draw) = c("trait", "pair")
to_draw = data.frame(to_draw)


pval = c()
pair = c()
trait = c()
conf_low = c()
conf_top = c()
stat = c()
test_type = c()

# #HYDROPATHICITY TESTS
# trait = c(trait, "Hydropathicity")
# pair = c(pair, "sQTL vs non_sQTL")
# test_type = c(test_type, "ks")
# out = ks_stat(data$HYDROPATHICITY,
#               data2$HYDROPATHICITY)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])

#LENGTH TESTS
pval = c(pval, test$p.value)
trait = c(trait, "Length")
pair = c(pair, "sQTL vs non_sQTL")
test_type = c(test_type, "ks")
out = ks_stat(data$LENGTH, data2$LENGTH)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

#ASN... TESTS

trait = c(trait, "% ASN")
pair = c(pair, "sQTL vs non_sQTL")
test_type = c(test_type, "ks")
out = ks_stat(data$ASN../data$LENGTH*100, 
              data2$ASN../data2$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


#CYS.. TESTS

trait = c(trait, "% CYS")
pair = c(pair, "sQTL vs non_sQTL")
test_type = c(test_type, "ks")
out = ks_stat(data$CYS../data$LENGTH*100, 
               data2$CYS../data2$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


# SYM TESTS

trait = c(trait, "Symm_dist")
pair = c(pair, "sQTL vs non_sQTL")
test_type = c(test_type, "ks")
out = ks_stat(data$LENGTH %% 3, 
              data2$LENGTH %% 3.)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

# #MIN TESTS
# 
# trait = c(trait, "MIN_RSA")
# pair = c(pair, "sQTL vs non_sQTL")
# test_type = c(test_type, "ks")
# out = ks_stat(data$MIN,
#               data2$MIN)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])

#Q1 TESTS

trait = c(trait, "Q1_RSA")
pair = c(pair, "sQTL vs non_sQTL")
test_type = c(test_type, "ks")
out = ks_stat(data$Q1,
              data2$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])
#
#
# #Q2 TESTS
# 
# trait = c(trait, "Q2_RSA")
# pair = c(pair, "sQTL vs non_sQTL")
# test_type = c(test_type, "ks")
# out = ks_stat(data$Q2,
#               data2$Q2)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# #Q3 TESTS
# 
# trait = c(trait, "Q3_RSA")
# pair = c(pair, "sQTL vs non_sQTL")
# test_type = c(test_type, "ks")
# out = ks_stat(data$Q3,
#               data2$Q3)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# #MAX TESTS
# 
# trait = c(trait, "MAX_RSA")
# pair = c(pair, "sQTL vs non_sQTL")
# test_type = c(test_type, "ks")
# out = ks_stat(data$MAX,
#               data2$MAX)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# #MIN_pLDDT TESTS
# 
# trait = c(trait, "MIN_pLDDT")
# pair = c(pair, "sQTL vs non_sQTL")
# test_type = c(test_type, "ks")
# out = ks_stat(data$MIN_pLLDT,
#               data2$MIN_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# #Q1_pLDDT TESTS
# 
# trait = c(trait, "Q1_pLDDT")
# pair = c(pair, "sQTL vs non_sQTL")
# test_type = c(test_type, "ks")
# out = ks_stat(data$Q1_pLLDT,
#               data2$Q1_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# #Q2_pLDDT TESTS
# 
# trait = c(trait, "Q2_pLDDT")
# pair = c(pair, "sQTL vs non_sQTL")
# test_type = c(test_type, "ks")
# out = ks_stat(data$Q2_pLLDT, data2$Q2_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])

#Q3_pLDDT TESTS

trait = c(trait, "Q3_pLDDT")
pair = c(pair, "sQTL vs non_sQTL")
test_type = c(test_type, "ks")
out = ks_stat(data$Q3_pLLDT, data2$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

# #MAX_pLDDT TESTS
# 
# trait = c(trait, "MAX_pLDDT")
# pair = c(pair, "sQTL vs non_sQTL")
# test_type = c(test_type, "ks")
# out = ks_stat(data$MAX_pLLDT, data2$MAX_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])

# #UNSTRUC TESTS
# unstruc = data.frame("sQTL" = c(nrow(data_unstruc), nrow(data) - nrow(data_unstruc)),
#                      "non_sQTL" = c(nrow(data2_unstruc), nrow(data2) - nrow(data2_unstruc)),
#                      #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                      row.names = c("unstruc", "other"),
#                      stringsAsFactors = FALSE)
#
# unstruc
# test = fisher.test(unstruc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Unstructured")
# pair = c(pair, "sQTL vs non_sQTL")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
#
# # STRUC TESTS
# struc = data.frame("sQTL" = c(nrow(data_struc), nrow(data)-nrow(data_struc)),
#                    "non_sQTL" = c(nrow(data2_struc), nrow(data2)-nrow(data2_struc)),
#                    #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                    row.names = c("struc", "other"),
#                    stringsAsFactors = FALSE)
#
# struc
# test = fisher.test(struc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Structured")
# pair = c(pair, "sQTL vs non_sQTL")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
#
# # CONSERVED TESTS
# cons = data.frame("sQTL" = c(nrow(data_cons), nrow(data) - nrow(data_cons)),
#                   "non_sQTL" = c(nrow(data2_cons), nrow(data2) - nrow(data2_cons)),
#                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                   row.names = c("cons", "other"),
#                   stringsAsFactors = FALSE)
#
# cons
# test = fisher.test(cons)
# pval = c(pval, test$p.value)
# trait = c(trait, "Conserved")
# pair = c(pair, "sQTL vs non_sQTL")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")

# SYM TESTS
symet = data.frame("sQTL" = c(length(data[data$LENGTH %% 3 == 0,]$LENGTH), length(data[data$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_sQTL" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
pval = c(pval, test$p.value)
trait = c(trait, "Symmetric")
pair = c(pair, "sQTL vs non_sQTL")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

pa = c("sQTL", "non_sQTL", "sQTL", "non_sQTL")
tra = c("Symmetric", "Symmetric", "Non symmetric", "Non symmetric")
count = c(length(data[data$LENGTH %% 3 == 0,]$LENGTH)/length(data$LENGTH),length(data2[data2$LENGTH %% 3 == 0,]$LENGTH)/length(data2$LENGTH), length(data[data$LENGTH %% 3 != 0,]$LENGTH)/length(data$LENGTH),  length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)/length(data2$LENGTH))
to_draw = cbind(tra, count)
to_draw = cbind(to_draw, pa)
to_draw
to_draw = data.frame(to_draw)
colnames(to_draw) = c("trait", "proportion", "Group")
to_draw$proportion = as.numeric(to_draw$proportion)
to_draw 
# ggplot(to_draw, aes(fill = as.factor(Group), x=trait, y=proportion )) + geom_col(position = position_dodge())  +  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("proportion") + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
# ggsave("symmetric_bars_der_all.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)


# SYM NOT FOUND TESTS
symet_not_found = data.frame("sQTL" = c(length(data[data$LENGTH %% 3 == 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data[data$LENGTH %% 3 != 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             "non_sQTL" = c(length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data2[data2$LENGTH %% 3 != 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             row.names = c("sym", "non_sym"),
                             stringsAsFactors = FALSE)
symet_not_found
test = fisher.test(symet_not_found)
pval = c(pval, test$p.value)
trait = c(trait, "Symm_not_found")
pair = c(pair, "sQTL vs non_sQTL")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

pa = c("sQTL", "non_sQTL", "sQTL", "non_sQTL")
tra = c("Symmetric", "Symmetric", "Non symmetric", "Non symmetric")
count = c(length(data[data$LENGTH %% 3 == 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data[data$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data2[data2$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH), length(data[data$LENGTH %% 3 != 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data[data$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),  length(data2[data2$LENGTH %% 3 !=  0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data2[data2$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH))
to_draw = data.frame()
to_draw = cbind(tra, count)
to_draw = cbind(to_draw, pa)
to_draw
to_draw = data.frame(to_draw)
colnames(to_draw) = c("trait", "proportion", "Group")
to_draw$proportion = as.numeric(to_draw$proportion)
to_draw 
#ggplot(to_draw, aes(fill = as.factor(Group), x=trait, y=proportion )) + geom_col(position = position_dodge())  +  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("proportion") + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
#ggsave("symmetric_not_found_bars_sQTL.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)



# HELIX TESTS
hel = data.frame("sQTL" = c(length(data[data$HELIX == "True",]$HELIX), length(data[data$HELIX != "True",]$HELIX)),
                 "non_sQTL" = c(length(data2[data2$HELIX == "True",]$HELIX), length(data2[data2$HELIX != "True",]$HELIX)),
                 #"non_sQTL" = c(length(data3[data3$HELIX ="True",]$HELIX), length(data3[data3$HELIX !"True",]$HELIX)),
                 row.names = c("helix", "no_helix"),
                 stringsAsFactors = FALSE)

hel
test = fisher.test(hel)
pval = c(pval, test$p.value)
trait = c(trait, "Helix")
pair = c(pair, "sQTL vs non_sQTL")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

# SHEETS TESTS
she = data.frame("sQTL" = c(length(data[data$SHEET == "True",]$SHEET), length(data[data$SHEET != "True",]$SHEET)),
                 "non_sQTL" = c(length(data2[data2$SHEET == "True",]$SHEET), length(data2[data2$SHEET != "True",]$SHEET)),
                 #"non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 row.names = c("sheet", "no_sheet"),
                 stringsAsFactors = FALSE)

she
test = fisher.test(she)
pval = c(pval, test$p.value)
trait = c(trait, "Sheet")
pair = c(pair, "sQTL vs non_sQTL")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

# TURNS TESTS
turn = data.frame("sQTL" = c(length(data[data$TURN == "True",]$TURN), length(data[data$TURN != "True",]$TURN)),
                  "non_sQTL" = c(length(data2[data2$TURN == "True",]$TURN), length(data2[data2$TURN != "True",]$TURN)),
                  #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("turn", "no_turn"),
                  stringsAsFactors = FALSE)

turn
test = fisher.test(turn)
pval = c(pval, test$p.value)
trait = c(trait, "Turn")
pair = c(pair, "sQTL vs non_sQTL")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

#TRANSMEMBRANE TEST

trans= data.frame("sQTL" = c(length(data[data$TRANSMEMBRANE == "True",]$LENGTH), length(data[data$TRANSMEMBRANE != "True",]$LENGTH)),
                  "non_sQTL" = c(length(data2[data2$TRANSMEMBRANE == "True",]$LENGTH), length(data2[data2$TRANSMEMBRANE != "True",]$LENGTH)),
                  #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                  row.names = c("trans", "non_trans"),
                  stringsAsFactors = FALSE)
trans
test = fisher.test(trans)
pval = c(pval, test$p.value)
trait = c(trait, "Transmembrane")
pair = c(pair, "sQTL vs non_sQTL")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

# #DOMAINS ANALYSIS
# ggplot() + geom_bar(data=data[data$DOMAIN != "['Disordered']" & data$DOMAIN != "[] ", ], aes(DOMAIN, fill="sQTL"), stat = "count", alpha=0.3) + geom_bar(data=data2[data2$DOMAIN != "['Disordered']" & data2$DOMAIN != "[]", ], aes(DOMAIN, fill="non sQTL"), stat = "count", alpha=0.3)
# table(data[data$DOMAIN != "['Disordered']" & data$DOMAIN != "[]", ]$DOMAIN)[table(data$DOMAIN) >= 2]
# table(data2[data2$DOMAIN != "Disordered  " & data2$DOMAIN != "", ]$DOMAIN)[table(data2$DOMAIN) >= 2]



# for (i in seq(1, length(pval), by=3)){
#           print(pval[i:(i+2)])
#           pval[i:(i+2)] = p.adjust(pval[i:(i+2)], method="hochberg")
# }

pval = p.adjust(pval, method="hochberg")

to_draw = data.frame()
to_draw = cbind(trait, log10(pval))
to_draw = cbind(to_draw, pair)
to_draw = cbind(to_draw, conf_low)
to_draw = cbind(to_draw, conf_top)
to_draw = cbind(to_draw, stat)
to_draw = cbind(to_draw, test_type)
#View(to_draw)
colnames(to_draw) = c("trait", "pval", "pair", "conf_low", "conf_top", "statistics", "test_type")
to_draw = data.frame(to_draw)
to_draw$pval = -as.numeric(to_draw$pval)
#to_draw
line = replicate(length(pval), 1.5)
ggplot(to_draw, aes(fill = as.factor(pair), x=trait, y=pval )) + geom_col(position = position_dodge())  + geom_hline(yintercept = line) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("-log10(p_value)") + labs(fill='Pairs') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5))) + gtex_v8_figure_theme()
ggsave('statistical_summary_without_correction_sQTL_non_sQTL.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
to_draw$conf_low = as.numeric(to_draw$conf_low)
to_draw$conf_top = as.numeric(to_draw$conf_top)
#to_draw$conf_low = ifelse(to_draw$conf_low < 0, 0, to_draw$conf_low)
#to_draw$conf_top = ifelse(to_draw$conf_top > 2, 2, to_draw$conf_top)
#to_draw$conf_low[is.nan(to_draw$conf_low)] = 0.0
write.csv(to_draw, "to_draw.csv", row.names = F)
line1 = replicate(length(pval), 0.0)
to_draw$alpha = ifelse((to_draw$test_type == "ks"  & ((to_draw$conf_low <= 0 & to_draw$conf_top >= 0) | 
                                                        (to_draw$conf_low >= 0 & to_draw$conf_top <= 0))) | 
                         (to_draw$test_type == "fisher"  & ((log(to_draw$conf_low) <= 0 & log(to_draw$conf_top) >= 0) | 
                                                              (log(to_draw$conf_low) >= 0 & log(to_draw$conf_top) <= 0))), yes = 0.2, 1.0)


ggplot(to_draw[to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_pointrange(aes(x = trait, y = log(as.numeric(statistics)), ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) +
        ylab("Enrichment in non sQTLs                                                                  Enrichment in sQTLs\n log(odd_ratio)") + 
        #geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) + coord_flip()  + 
  geom_hline(yintercept = line1, color="red") + guides(alpha = FALSE) + gtex_v8_figure_theme()
ggsave('statistical_summary_fisher_sQTL_non_sQTL.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)

line1 = replicate(length(pval), 0.0)
ggplot(to_draw[to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics))) +
        geom_pointrange(aes(x = trait, y = as.numeric(statistics), ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) +
        ylab("Enrichment in non sQTLs                                                                  Enrichment in sQTLs\n m1-m2") + 
        # geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + coord_flip() + 
  geom_hline(yintercept = line1, color="red") + gtex_v8_figure_theme() + guides(alpha = FALSE)


ggsave('statistical_summary_ks_sQTL_non_sQTL.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)


View(to_draw)

ggplot() + geom_density(data=data[data$DIP_P <= 0.05, ], aes(DIP, fill="sQTL"), alpha = 0.2) + geom_density(data=data2[data2$DIP_P <= 0.05, ], aes(DIP, fill="non sQTL"), alpha = 0.2)
ggplot() + geom_density(data=data[data$DIP_P <= 0.05, ], aes(-log10(DIP_P), fill="sQTL"), alpha = 0.2) + geom_density(data=data2[data2$DIP_P <= 0.05, ], aes(-log10(DIP_P), fill="non sQTL"), alpha = 0.2)
min(data2$DIP_P, na.rm=T)
View(data[data$DIP_P <= 0.001 & !is.na(data$DIP_P), ])



write.csv(to_draw, "Data/to_draw_sQTL.csv", quote=F, row.names = F)
# 
#
# tr = "_pLDDT"
# pa= "disordered vs structured"
# te_ty = "ks"
# out = ks_stat(data[grepl("['Disordered']", data$DOMAIN, fixed = TRUE) & data$HELIX == 'False' & data$SHEET == 'False' & data$TURN == 'False',]$MIN_pLLDT,
#               data[!grepl("['Disordered']", data$DOMAIN, fixed = TRUE) & (data$HELIX == 'True' | data$SHEET == 'True' | data$TURN == 'True'),]$MIN_pLLDT)
# pv = out[1]
# c_l = out[3]
# c_t =  out[4]
# st = out[2]
# st
# #mean(inv_norm(data[grepl("['Disordered']", data$DOMAIN, fixed = TRUE) & data$HELIX == 'False' & data$SHEET == 'False' & data$TURN == 'False',]$MIN_pLLDT), na.rm = T) - mean(inv_norm(data[!grepl("['Disordered']", data$DOMAIN, fixed = TRUE) & (data$HELIX == 'True' | data$SHEET == 'True' | data$TURN == 'True'),]$MIN_pLLDT), na.rm=T)
# to_draw = data.frame()
# to_draw = cbind(tr, log10(pval))
# to_draw = cbind(to_draw, pa)
# to_draw = cbind(to_draw, c_l)
# to_draw = cbind(to_draw, c_t)
# to_draw = cbind(to_draw, st)
# to_draw = cbind(to_draw, te_ty)
# to_draw
# colnames(to_draw) = c("trait", "pval", "pair", "conf_low", "conf_top", "statistics", "test_type")
# to_draw = data.frame(to_draw)
# to_draw$pval = -as.numeric(to_draw$pval)
# to_draw
#
#
# ggplot(to_draw[to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics))) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + ylab("m1-m2") + coord_flip() + geom_hline(yintercept = line1, color="red")
#
# #ggsave('statistical_summary_ks_sQTL_non_sQTL.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
#




data = read.csv("Data/output_top_sQTL.csv")
data2 = read.csv("Data/output_non_sQTL_full.csv")

#
pre_data_1 = data[, 19:36]
pre_data_1$LENGTH = data$LENGTH
#pre_data_1$LENGTH.ALIGN = data$LENGTH.ALIGN
#pre_data_1$STRUC = as.numeric(data$GeneID %in% data_struc$GeneID)
#pre_data_1$UNSTRUC = as.numeric(data$GeneID %in% data_unstruc$GeneID)
#pre_data_1$TRANSMEMBRANE = as.numeric(data$TRANSMEMBRANE == "True")
# pre_data_1$HELIX = as.numeric(data$HELIX == "True")
# pre_data_1$SHEET = as.numeric(data$SHEET == "True")
# pre_data_1$TURN = as.numeric(data$TURN == "True")
# pre_data_1$ASN = data$ASN../data$LENGTH *100
# pre_data_1$CYS = data$CYS../data$LENGTH *100
# pre_data_1$SIGNAL = as.numeric(data$SIGNAL != "[]" & data$SIGNAL != "['']" & data$SIGNAL != "['', '']")
pre_data_1$GROUP = replicate(nrow(pre_data_1), "sQTL")
#
#
pre_data_2 = data2[, 21:32]
pre_data_2$LENGTH = data2$LENGTH
#pre_data_2$LENGTH.ALIGN = data2$LENGTH.ALIGN
#pre_data_2$STRUC = as.numeric(data2$gene_id %in% data2_struc$gene_id)
#pre_data_2$UNSTRUC = as.numeric(data2$gene_id %in% data2_unstruc$gene_id)
# pre_data_2$TRANSMEMBRANE = as.numeric(data2$TRANSMEMBRANE == "True")
# pre_data_2$HELIX = as.numeric(data2$HELIX == "True")
# pre_data_2$SHEET = as.numeric(data2$SHEET == "True")
# pre_data_2$TURN = as.numeric(data2$TURN == "True")
# pre_data_2$ASN = data2$ASN../data2$LENGTH *100
# pre_data_2$CYS = data2$CYS../data2$LENGTH *100
# pre_data_2$SIGNAL = as.numeric(data2$SIGNAL != "[]" & data2$SIGNAL != "['']" & data2$SIGNAL != "['', '']")
pre_data_2$GROUP = replicate(nrow(pre_data_2), "non_sQTL")
#pre_data_2 = sample_n(pre_data_2, 200)
#
colnames(pre_data_1)
colnames(pre_data_2)
data_for_pca = rbind(pre_data_1[,colnames(pre_data_2)], pre_data_2)
data_for_pca = na.omit(data_for_pca)
data_for_pca.pca = prcomp(data_for_pca[, 1:13], center = T, scale. = T)
data_for_pca.pca
#View(data_for_pca)
#
library(ggfortify)
autoplot(data_for_pca.pca, data = data_for_pca, colour = 'GROUP')
autoplot(data_for_pca.pca, data = data_for_pca, colour = 'DIP')
#
wss <- (nrow(data_for_pca[, 1:10])-1)*sum(apply(data_for_pca[, 1:10],2,var))
for (i in 2:15) wss[i] <- sum(kmeans(data_for_pca[, 1:10],
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
comp <- data.frame(data_for_pca.pca$x[,1:4])
k <- kmeans(comp, 2, nstart=25, iter.max=1000)
library(RColorBrewer)
library(scales)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(comp, col=k$clust, pch=16)
comp$DIP = data_for_pca$DIP
comp$GROUP = data_for_pca$GROUP
ggplot(comp, aes(x=PC1, y=PC2)) +
        geom_point(aes(alpha=as.numeric(DIP), color=as.factor(k$cluster)), size=2.5) + scale_shape_manual(values=c(3, 25, 16))
ggsave("PCA_try.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
#
#
# ggplot(comp, aes(x=PC1, y=PC2)) +
#         geom_point(aes(size=as.factor(UNSTRUC), color=as.factor(k$cluster), shape=as.factor(DIP))) + scale_shape_manual(values=c(3, 13, 17))
sort(table(k$clust))
clust <- names(sort(table(k$clust)))
#First cluster
table(data_for_pca[k$clust==clust[1], ]$GROUP)
#Second Cluster
table(data_for_pca[k$clust==clust[2], ]$GROUP)
table(data_for_pca[k$clust==clust[3], ]$GROUP)
# table(data_for_pca[k$clust==clust[4], ]$STRUC)
# table(data_for_pca[k$clust==clust[1], ]$UNSTRUC)
# table(data_for_pca[k$clust==clust[2], ]$UNSTRUC)
# table(data_for_pca[k$clust==clust[3], ]$UNSTRUC)
# #
# table(data_for_pca[k$clust==clust[1], ]$STRUC)
# table(data_for_pca[k$clust==clust[2], ]$STRUC)
# table(data_for_pca[k$clust==clust[3], ]$STRUC)
# table(data_for_pca[k$clust==clust[4], ]$UNSTRUC)
# #
mean1 = numeric(0)
max1 = numeric(0)
for (i in 1:(ncol(pre_data_1)-1)){
        mean1 = c(mean1, mean(pre_data_1[,i], na.rm = T))
        max1 = c(max1, max(pre_data_1[,i], na.rm = T))
}
mean1
#
mean2 = numeric(0)
max2 = numeric(0)
for (i in 1:(ncol(pre_data_2)-1)){
        mean2 = c(mean2, mean(pre_data_2[,i], na.rm = T))
        max2 = c(max2, max(pre_data_2[,i], na.rm = T))
}
#
#
mean1 = mean1/max1
mean2 = mean2/max2
data_mean = rbind(mean1, mean2, mean3)
data_mean
colnames(data_mean) = colnames(pre_data_1)[1:(ncol(pre_data_1)-1)]
data_mean
ggplot(data_mean, aes(x=replicate(3, "MIN"), y=MIN)) + geom_col(fill=c(1,2,3), position = position_dodge())  +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("-log10(p_value)") + labs(fill='Pairs') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
# #ggsave("statistical_summary_with_correction.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
#
