setwd("~/splicing_project/")
###___________________________________________________CHANGED______________________________________________________

library(ggplot2)
library(stats)
library(dplyr)
library(boot)
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

medianDiff = function(data, func){
  data = func(data)
  m1 <- median(data[data$V2=='sQTL', "totalData"])
  m2 <- median(data[data$V2=='non_sQTL', "totalData"])
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
  # dist = c()
  # for (i in 1:R){
  #   print(i)
  #   dist = c(dist, func(x, func1))
  # }
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
  #x = data[grepl("['Disordered']", data$DOMAIN, fixed = TRUE),]$MIN_pLLDT
  x = x[!is.na(x)]
  #y = data[!grepl("['Disordered']", data$DOMAIN, fixed = TRUE),]$MIN_pLLDT
  y = y[!is.na(y)]
  pval =  ks.test(x, y)$p.value
  totalData = c(x, y)
  totalData = cbind(totalData, c(replicate(length(x), "sQTL"), replicate(length(y), "non_sQTL")))
  #View(totalData)
  totalData[,1] = inv_norm(as.numeric(totalData[,1]))
  totalData = data.frame(totalData)
  totalData = totalData[!is.na(totalData[,1]), ]
  totalData[,1] = as.numeric(totalData[,1])
  print(median(totalData[totalData$V2 == "sQTL", 'totalData']) - median(totalData[totalData$V2 == "non_sQTL", 'totalData']))
  
  totalBoot = my_bootstrap(totalData, R=1000, medianDiff, select_samples)
  #totalBoot = boot(totalData, statistic = medianDiff, R = 10000)
  totalBootCI = quantile(totalBoot, probs=c(0.05, 0.95, 0.5))
  conf_low = totalBootCI[1]
  conf_top = totalBootCI[2]
  #stat = mean(totalData[totalData$V2 == "sQTL", 'totalData']) - mean(totalData[totalData$V2 == "non_sQTL", 'totalData'])
  stat = totalBootCI[3]
  print(c(pval, stat, conf_low, conf_top))
  return(c(pval, stat, conf_low, conf_top))
}

s_ks = c()
c_l_ks = c()
c_t_ks = c()
s_f = c()
c_l_f = c()
c_t_f = c()
## READ AND SELECT COLOC DATA
#for (i in 1:2){
  #set.seed(i)
  # data_full = read.csv("Data/output_data_coloc.csv")
  # filt = readLines("Data/true_coloc_new_new.txt")
  # filt
  # data_full = data_full[data_full$phenotype_id %in% filt, ]
  # data = data_full[,24:ncol(data_full)]
  # data = subset(data, select = -c(Exon_coord) )
  #

  data = read.csv("Data/output_data_coloc.csv")
  ## READ AND SELECT NON_COLOC DATA
  data2 = read.csv("Data/output_data_non_coloc.csv")
  data2 = data2[,24:ncol(data2)]
  # new_data = read.table("Data/top_sQTLs_top_coloc.tsv", header = T, sep="\t")
  # datafull = data_full[order(phenotype_id),]
  #write.csv(data_full, "Data/output_data_coloc_top_new.csv", row.names = F)
  # new_data = new_data[new_data$top_pid %in% data_full$phenotype_id, ]
  # datafull = data_full[!duplicated(data_full$phenotype_id), ]
  # datafull = datafull[order(datafull$phenotype_id),]
  # newdata = new_data[order(datafull$phenotype_id),]
  # to_draw = cbind(newdata$mean_01_psi, datafull$ALIGNED.SEQ)
  # to_draw = data.frame(to_draw)
  # colnames(to_draw) = c("mean_01_psi", "ALIGNED.SEQ")
  # to_draw$mean_01_psi = as.numeric(to_draw$mean_01_psi)
  # to_draw$FOUND = (to_draw$ALIGNED.SEQ != "EXON NOT FOUND")  
  # table(to_draw$FOUND)
  # View(newdata)
  # ggplot() + geom_boxplot(data=to_draw, aes(mean_01_psi, fill=FOUND))
  # #ggsave("box_mean_psi_found_notfound.png", height = 3.11, width = 7.92,path = "Data/visuals/", device='png', dpi=700)
  # ks.test(to_draw[to_draw$FOUND == T,]$mean_01_psi, to_draw[to_draw$FOUND != T,]$mean_01_psi)
  # ggplot() + geom_density(data=to_draw, aes(mean_01_psi, fill=FOUND), alpha=0.2)
  #ggsave("mean_psi_found_notfound.png", height = 3.11, width = 7.92,path = "Data/visuals/", device='png', dpi=700)
  
  
  # View(data_full)
  nrow(data)
  nrow(data2)
  
  
  ## SELECT GENES WITH ONLY ONE SIGNIFICANT EXON
  # filt = readLines("Data/genes_with_only_one_significant_exon.txt")
  # filt
  # data = data[data$GeneID %in% filt, ]
  
  ## REMOVING NONUNIQUE ROWS
  data = unique.data.frame(data)
  data2 = unique.data.frame(data2)
  
  nrow(data)
  nrow(data2)
  data = data[!flag_outliers(data$LENGTH) & data$LENGTH >= 15,]
  data2 = data2[!flag_outliers(data2$LENGTH) & data2$LENGTH >= 15,]
  # data3 = data3[!flag_outliers(data3$LENGTH),]
  nrow(data)
  nrow(data2)
  min(data$LENGTH)
  min(data2$LENGTH)

  # View(data)
  # View(data[!flag_outliers(data$LENGTH) & data$LENGTH!=0 ,])
  #dev.off()
  
  
  # asn_table = cbind(data$LENGTH, data$ASN..)
  # asn_table = data.frame(asn_table)
  # asn_table$GROUP = "coloc"
  # colnames(asn_table) = c("LENGTH", "#ASN", "GROUP")
  # write.csv(asn_table, "Data/raw_asn_data_coloc.csv", row.names = F)
  # asn_table2 = cbind(data2$LENGTH, data2$ASN..)
  # asn_table2 = data.frame(asn_table2)
  # asn_table2$GROUP = "non_coloc"
  # colnames(asn_table2) = c("LENGTH", "#ASN", "GROUP")
  # write.csv(asn_table2, "Data/raw_asn_data_non_coloc.csv", row.names = F)
  # asn_table3 = cbind(data3$LENGTH, data3$ASN..)
  # asn_table3 = data.frame(asn_table3)
  # asn_table3$GROUP = "non_sQTL"
  # colnames(asn_table3) = c("LENGTH", "#ASN", "GROUP")
  # asn_table = rbind(asn_table, asn_table2)
  # asn_table = rbind(asn_table, asn_table3)
  # #colnames(asn_table) = c("LENGTH", "#ASN", "GROUP")
  # #View(asn_table)
  # write.csv(asn_table3, "Data/raw_asn_data_non_sQTL.csv", row.names = F)
  # ggplot(data=asn_table, aes(x=`#ASN`/LENGTH, color=GROUP)) + geom_histogram()

  ggplot() + geom_density(data=data[data$LENGTH!=0,], aes(LENGTH, fill="coloc"), alpha = 0.2) + geom_density(data=data2[data2$LENGTH!=0,], aes(LENGTH, fill="non coloc"), alpha = 0.2) 
  mean(data[data$LENGTH!=0,]$LENGTH)
  mean(data2[data2$LENGTH!=0 & !is.na(data2$LENGTH),]$LENGTH)
  ggplot() + geom_density(data=data[!flag_outliers(data$HYDROPATHICITY),], aes(HYDROPATHICITY, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$HYDROPATHICITY),], aes(HYDROPATHICITY, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data[!flag_outliers(data$ASN../data$LENGTH *100),], aes(ASN../LENGTH *100, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$ASN../data2$LENGTH *100),], aes(ASN../LENGTH * 100, fill="non coloc"), alpha = 0.2) 
  ggplot() + geom_density(data=data[!flag_outliers(data$CYS../data$LENGTH *100),], aes(CYS../LENGTH *100, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$CYS../data2$LENGTH *100),], aes(CYS../LENGTH * 100, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data, aes(CYS../LENGTH *100, fill="coloc"), alpha = 0.2) + geom_density(data=data2, aes(CYS../LENGTH * 100, fill="non coloc"), alpha = 0.2) 
  
  # data = cbind(data, data$LENGTH %% 3)
  # data2 = cbind(data2, data2$LENGTH %% 3)
  # data3 = cbind(data3, data3$LENGTH %% 3)
  
  ggplot() + geom_density(data=data[!flag_outliers(data$LENGTH %% 3),], aes(LENGTH %% 3, fill="coloc"), alpha=0.3) + geom_density(data=data2[!flag_outliers(data2$LENGTH %% 3),], aes(LENGTH %% 3, fill="non coloc"), alpha=0.3)
  
  
  ggplot() + geom_density(data=data[!flag_outliers(data$MIN),], aes(MIN, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$MIN),], aes(MIN, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data[!flag_outliers(data$Q1),], aes(Q1, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q1),], aes(Q1, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data[!flag_outliers(data$Q2),], aes(Q2, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q2),], aes(Q2, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data[!flag_outliers(data$Q3),], aes(Q3, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q3),], aes(Q3, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data[!flag_outliers(data$MAX),], aes(MAX, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$MAX),], aes(MAX, fill="non coloc"), alpha = 0.2)
  
  
  ggplot() + geom_density(data=data[!flag_outliers(data$MIN_pLLDT),], aes(MIN_pLLDT, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$MIN_pLLDT),], aes(MIN_pLLDT, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data[!flag_outliers(data$Q1_pLLDT),], aes(Q1_pLLDT, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q1_pLLDT),], aes(Q1_pLLDT, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data[!flag_outliers(data$Q2_pLLDT),], aes(Q2_pLLDT, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q2_pLLDT),], aes(Q2_pLLDT, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data[!flag_outliers(data$Q3_pLLDT),], aes(Q3_pLLDT, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$Q3_pLLDT),], aes(Q3_pLLDT, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data[!flag_outliers(data$MAX_pLLDT),], aes(MAX_pLLDT, fill="coloc"), alpha = 0.2) + geom_density(data=data2[!flag_outliers(data2$MAX_pLLDT),], aes(MAX_pLLDT, fill="non coloc"), alpha = 0.2)
  
  
  
  
  ggplot() + geom_density(data=data, aes(HELIX, fill="coloc"), alpha = 0.2) + geom_density(data=data2, aes(HELIX, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data, aes(SHEET, fill="coloc"), alpha = 0.2) + geom_density(data=data2, aes(SHEET, fill="non coloc"), alpha = 0.2)
  ggplot() + geom_density(data=data, aes(TURN, fill="coloc"), alpha = 0.2) + geom_density(data=data2, aes(TURN, fill="non coloc"), alpha = 0.2)
  
  data_unstruc = data[data$Q2 > 40.0 & data$Q2_pLLDT < 50.0, ]
  data_unstruc = data_unstruc[!is.na(data_unstruc$GeneID),] # 84/571 = 14,7% 67 signal
  nrow(data_unstruc)
  data2_unstruc = data2[data2$Q2 > 40.0 & data2$Q2_pLLDT < 50.0, ]
  data2_unstruc = data2_unstruc[!is.na(data2_unstruc$gene_id),] # 106/686 = 15,5% 68 signal
  nrow(data2_unstruc)

  
  data_struc = data[data$Q2 < 25.0 & data$Q2_pLLDT >= 90.0, ]
  data_struc = data_struc[!is.na(data_struc$GeneID),] # 105/571 = 18,4% 69 signals
  nrow(data_struc)
  data2_struc = data2[data2$Q2 < 25.0 & data2$Q2_pLLDT >= 90.0, ]
  data2_struc = data2_struc[!is.na(data2_struc$gene_id),] # 108/686 = 15,7% 66 signals
  nrow(data2_struc)
  
  data_cons = data[data$Q2_pLLDT >= 90.0, ]
  data_cons = data_cons[!is.na(data_cons$GeneID),] # 230/571 = 40,1% 156 signals
  nrow(data_cons)
  data2_cons = data2[data2$Q2_pLLDT >= 90.0, ]
  data2_cons = data2_cons[!is.na(data2_cons$gene_id),] # 246/686 = 35,8% 162 signals
  nrow(data2_cons)
  
  data_loop = data[data$TURN == "True" & data$HELIX == "False" & data$SHEET == "False", ]
  data_loop = data_loop[!is.na(data_loop$GeneID),] # 59/571 = 10,3 % 37 signal
  nrow(data_loop)
  data2_loop = data2[data2$TURN == "True" & data2$HELIX == "False" & data2$SHEET == "False", ]
  data2_loop = data2_loop[!is.na(data2_loop$gene_id),] # 54/686 = 7,9 % 32 signal
  nrow(data2_loop)
 
  #dev.off()
  to_draw = data.frame()
  round(nrow(data_unstruc) / nrow(data) * 100)
  pa = replicate(round(nrow(data_unstruc) / nrow(data) * 100), "coloc")
  pa = c(pa, replicate(round(nrow(data2_unstruc) / nrow(data2) * 100), "non_coloc"))
  pa = c(pa, replicate(round(nrow(data_struc) / nrow(data) * 100), "coloc"))
  pa = c(pa, replicate(round(nrow(data2_struc) / nrow(data2) * 100), "non_coloc"))
  pa = c(pa, replicate(round(nrow(data_cons) / nrow(data) * 100), "coloc"))
  pa = c(pa, replicate(round(nrow(data2_cons) / nrow(data2) * 100), "non_coloc"))
  pa = c(pa, replicate(round(nrow(data_loop) / nrow(data) * 100), "coloc"))
  pa = c(pa, replicate(round(nrow(data2_loop) / nrow(data2) * 100), "non_coloc"))
  trait = replicate(round(nrow(data_unstruc) / nrow(data)*100) + round(nrow(data2_unstruc) / nrow(data2)*100) , "Unstructured")
  trait = c(trait, replicate(round(nrow(data_struc) / nrow(data)*100) + round(nrow(data2_struc) / nrow(data2)*100)), "Structured")
  trait = c(trait, replicate(round(nrow(data_cons) / nrow(data)*100) + round(nrow(data2_cons) / nrow(data2)*100)), "Conserved")
  trait = c(trait, replicate(round(nrow(data_loop) / nrow(data)*100) + round(nrow(data2_loop) / nrow(data2)*100)), "Loop")
  print(length(pa))
  print(length(trait))
  to_draw = cbind(trait, pa)
  to_draw
  colnames(to_draw) = c("trait", "pair")
  to_draw = data.frame(to_draw)
  to_draw
  ggplot(to_draw, aes(fill = as.factor(pa), x=trait)) + geom_bar(position = position_dodge())  +theme(axis.text.x = element_text( vjust = 0.5, hjust=1, size = 10), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8))  + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5))) + ylab("proportion")
  #ggsave("structure_bars.png", height = 3.11, width = 7.92,path = "Data/visuals/", device='png', dpi=700)
  
  pval = c()
  pair = c()
  trait = c()
  test_type = c()
  conf_low = c()
  conf_top = c()
  stat = c()
  # 
  # #HYDROPATHICITY TESTS
  # trait = c(trait, "Hydropathicity")
  # pair = c(pair, "coloc vs non_coloc")
  # test_type = c(test_type, "ks")
  # out = ks_stat(data[!flag_outliers(data$HYDROPATHICITY),]$HYDROPATHICITY, 
  #               data2[!flag_outliers(data2$HYDROPATHICITY),]$HYDROPATHICITY)
  # pval = c(pval, out[1])
  # conf_low = c(conf_low, out[3])
  # conf_top = c(conf_top, out[4])
  # stat = c(stat, out[2])
  # 
  
  # LENGTH TESTS
  
  trait = c(trait, "Length")
  pair = c(pair, "coloc vs non_coloc")
  test_type = c(test_type, "ks")
  out = ks_stat(data$LENGTH, data2$LENGTH)
  pval = c(pval, out[1])
  conf_low = c(conf_low, out[3])
  conf_top = c(conf_top, out[4])
  stat = c(stat, out[2])
  
  # ASN... TESTS
  
  trait = c(trait, "% ASN")
  pair = c(pair, "coloc vs non_coloc")
  test_type = c(test_type, "ks")
  out = ks_stat(data[!flag_outliers(data$ASN../data$LENGTH *100),]$ASN../data$LENGTH*100, 
                data2[!flag_outliers(data2$ASN../data2$LENGTH *100),]$ASN../data2$LENGTH*100)
  pval = c(pval, out[1])
  conf_low = c(conf_low, out[3])
  conf_top = c(conf_top, out[4])
  stat = c(stat, out[2])

  
  # CYS.. TESTS
  
  trait = c(trait, "% CYS")
  pair = c(pair, "coloc vs non_coloc")
  test_type = c(test_type, "ks")
  out = ks_stat(data[!flag_outliers(data$CYS../data$LENGTH *100),]$CYS../data$LENGTH*100, 
                data2[!flag_outliers(data2$CYS../data2$LENGTH *100),]$CYS../data2$LENGTH*100)
  pval = c(pval, out[1])
  conf_low = c(conf_low, out[3])
  conf_top = c(conf_top, out[4])
  stat = c(stat, out[2])

  
  # SYM TESTS
  
  trait = c(trait, "Symm_dist")
  pair = c(pair, "coloc vs non_coloc")
  test_type = c(test_type, "ks")
  out = ks_stat(data[!flag_outliers(data$LENGTH %% 3),]$LENGTH %% 3, 
                data2[!flag_outliers(data2$LENGTH %% 3),]$LENGTH %% 3.)
  pval = c(pval, out[1])
  conf_low = c(conf_low, out[3])
  conf_top = c(conf_top, out[4])
  stat = c(stat, out[2])

  
  # MEAN.ACC TESTS
  # ks.test(data[!flag_outliers(data$MEAN.ACC),]$MEAN.ACC, data2[!flag_outliers(data2$MEAN.ACC),]$MEAN.ACC)
  # ks.test(data[!flag_outliers(data$MEAN.ACC),]$MEAN.ACC, data3[!flag_outliers(data3$MEAN.ACC),]$MEAN.ACC)
  # ks.test(data3[!flag_outliers(data3$MEAN.ACC),]$MEAN.ACC, data2[!flag_outliers(data2$MEAN.ACC),]$MEAN.ACC)
  
  # #MIN TESTS
  # 
  # trait = c(trait, "MIN_RSA")
  # pair = c(pair, "coloc vs non_coloc")
  # test_type = c(test_type, "ks")
  # out = ks_stat(data[!flag_outliers(data$MIN),]$MIN, data2[!flag_outliers(data2$MIN),]$MIN)
  # pval = c(pval, out[1])
  # conf_low = c(conf_low, out[3])
  # conf_top = c(conf_top, out[4])
  # stat = c(stat, out[2])

  #Q1 TESTS

  trait = c(trait, "Q1_RSA")
  pair = c(pair, "coloc vs non_coloc")
  test_type = c(test_type, "ks")
  out = ks_stat(data[!flag_outliers(data$Q1),]$Q1, data2[!flag_outliers(data2$Q1),]$Q1)
  pval = c(pval, out[1])
  conf_low = c(conf_low, out[3])
  conf_top = c(conf_top, out[4])
  stat = c(stat, out[2])
 
  
  # #Q2 TESTS
  # 
  # trait = c(trait, "Q2_RSA")
  # pair = c(pair, "coloc vs non_coloc")
  # test_type = c(test_type, "ks")
  # out = ks_stat(data[!flag_outliers(data$Q2),]$Q2, data2[!flag_outliers(data2$Q2),]$Q2)
  # pval = c(pval, out[1])
  # conf_low = c(conf_low, out[3])
  # conf_top = c(conf_top, out[4])
  # stat = c(stat, out[2])
  
  
  # #Q3 TESTS
  # 
  # trait = c(trait, "Q3_RSA")
  # pair = c(pair, "coloc vs non_coloc")
  # test_type = c(test_type, "ks")
  # out = ks_stat(data[!flag_outliers(data$Q3),]$Q3, data2[!flag_outliers(data2$Q3),]$Q3)
  # pval = c(pval, out[1])
  # conf_low = c(conf_low, out[3])
  # conf_top = c(conf_top, out[4])
  # stat = c(stat, out[2])
  # 
  # 
  # #MAX TESTS
  # 
  # trait = c(trait, "MAX_RSA")
  # pair = c(pair, "coloc vs non_coloc")
  # test_type = c(test_type, "ks")
  # out = ks_stat(data[!flag_outliers(data$MAX),]$MAX, data2[!flag_outliers(data2$MAX),]$MAX)
  # pval = c(pval, out[1])
  # conf_low = c(conf_low, out[3])
  # conf_top = c(conf_top, out[4])
  # stat = c(stat, out[2])
  # 
  # #MIN_pLDDT TESTS
  # 
  # trait = c(trait, "MIN_pLDDT")
  # pair = c(pair, "coloc vs non_coloc")
  # test_type = c(test_type, "ks")
  # out = ks_stat(data[!flag_outliers(data$MIN_pLLDT),]$MIN_pLLDT, 
  #               data2[!flag_outliers(data2$MIN_pLLDT),]$MIN_pLLDT)
  # pval = c(pval, out[1])
  # conf_low = c(conf_low, out[3])
  # conf_top = c(conf_top, out[4])
  # stat = c(stat, out[2])
  # 
  # #Q1_pLDDT TESTS
  # 
  # trait = c(trait, "Q1_pLDDT")
  # pair = c(pair, "coloc vs non_coloc")
  # test_type = c(test_type, "ks")
  # out = ks_stat(data[!flag_outliers(data$Q1_pLLDT),]$Q1_pLLDT, 
  #               data2[!flag_outliers(data2$Q1_pLLDT),]$Q1_pLLDT)
  # pval = c(pval, out[1])
  # conf_low = c(conf_low, out[3])
  # conf_top = c(conf_top, out[4])
  # stat = c(stat, out[2])
  
  
  # #Q2_pLDDT TESTS
  # 
  # trait = c(trait, "Q2_pLDDT")
  # pair = c(pair, "coloc vs non_coloc")
  # test_type = c(test_type, "ks")
  # out = ks_stat(data[!flag_outliers(data$Q2_pLLDT),]$Q2_pLLDT, data2[!flag_outliers(data2$Q2_pLLDT),]$Q2_pLLDT)
  # pval = c(pval, out[1])
  # conf_low = c(conf_low, out[3])
  # conf_top = c(conf_top, out[4])
  # stat = c(stat, out[2])
  
  #Q3_pLDDT TESTS

  trait = c(trait, "Q3_pLDDT")
  pair = c(pair, "coloc vs non_coloc")
  test_type = c(test_type, "ks")
  out = ks_stat(data[!flag_outliers(data$Q3_pLLDT),]$Q3_pLLDT, data2[!flag_outliers(data2$Q3_pLLDT),]$Q3_pLLDT)
  pval = c(pval, out[1])
  conf_low = c(conf_low, out[3])
  conf_top = c(conf_top, out[4])
  stat = c(stat, out[2])


  # #MAX_pLDDT TESTS
  # 
  # trait = c(trait, "MAX_pLDDT")
  # pair = c(pair, "coloc vs non_coloc")
  # test_type = c(test_type, "ks")
  # out = ks_stat(data[!flag_outliers(data$MAX_pLLDT),]$MAX_pLLDT, data2[!flag_outliers(data2$MAX_pLLDT),]$MAX_pLLDT)
  # pval = c(pval, out[1])
  # conf_low = c(conf_low, out[3])
  # conf_top = c(conf_top, out[4])
  # stat = c(stat, out[2])
  # 
  # 
  # # UNSTRUC TESTS
  # unstruc = data.frame("coloc" = c(nrow(data_unstruc), nrow(data) - nrow(data_unstruc)),
  #                    "non_coloc" = c(nrow(data2_unstruc), nrow(data2) - nrow(data2_unstruc)),
  #                    #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
  #                    row.names = c("unstruc", "other"),
  #                    stringsAsFactors = FALSE)
  # 
  # unstruc
  # test = fisher.test(unstruc)
  # pval = c(pval, test$p.value)
  # trait = c(trait, "Unstructured")
  # pair = c(pair, "coloc vs non_coloc")
  # conf_low = c(conf_low, test$conf.int[1])
  # conf_top = c(conf_top, test$conf.int[2])
  # stat = c(stat, test$estimate)
  # test_type = c(test_type, "fisher")
  # 
  # # 
  # # STRUC TESTS
  # struc = data.frame("coloc" = c(nrow(data_struc), nrow(data)-nrow(data_struc)),
  #                      "non_coloc" = c(nrow(data2_struc), nrow(data2)-nrow(data2_struc)),
  #                      #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
  #                      row.names = c("struc", "other"),
  #                      stringsAsFactors = FALSE)
  # 
  # struc
  # test = fisher.test(struc)
  # pval = c(pval, test$p.value)
  # trait = c(trait, "Structured")
  # pair = c(pair, "coloc vs non_coloc")
  # conf_low = c(conf_low, test$conf.int[1])
  # conf_top = c(conf_top, test$conf.int[2])
  # stat = c(stat, test$estimate)
  # test_type = c(test_type, "fisher")
  # 
  # 
  # # CONSERVED TESTS
  # cons = data.frame("coloc" = c(nrow(data_cons), nrow(data) - nrow(data_cons)),
  #                    "non_coloc" = c(nrow(data2_cons), nrow(data2) - nrow(data2_cons)),
  #                    #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
  #                    row.names = c("cons", "other"),
  #                    stringsAsFactors = FALSE)
  # 
  # cons
  # test = fisher.test(cons)
  # pval = c(pval, test$p.value)
  # trait = c(trait, "Conserved")
  # pair = c(pair, "coloc vs non_coloc")
  # conf_low = c(conf_low, test$conf.int[1])
  # conf_top = c(conf_top, test$conf.int[2])
  # stat = c(stat, test$estimate)
  # test_type = c(test_type, "fisher")
  
  
  # SYM TESTS
  symet = data.frame("coloc" = c(length(data[data$LENGTH %% 3 == 0,]$LENGTH), length(data[data$LENGTH %% 3 != 0,]$LENGTH)),
                     "non_coloc" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
                     #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                     row.names = c("sym", "non_sym"),
                     stringsAsFactors = FALSE)
  
  symet
  test = fisher.test(symet)
  pval = c(pval, test$p.value)
  trait = c(trait, "Symmetric")
  pair = c(pair, "coloc vs non_coloc")
  conf_low = c(conf_low, test$conf.int[1])
  conf_top = c(conf_top, test$conf.int[2])
  stat = c(stat, test$estimate)
  test_type = c(test_type, "fisher")
  
  
  pa = c("coloc", "non_coloc",  "coloc", "non_coloc",)
  tra = c("Symmetric", "Symmetric", "Non symmetric", "Non symmetric")
  count = c(length(data[data$LENGTH %% 3 == 0,]$LENGTH)/length(data$LENGTH),length(data2[data2$LENGTH %% 3 == 0,]$LENGTH)/length(data2$LENGTH), length(data[data$LENGTH %% 3 != 0,]$LENGTH)/length(data$LENGTH),  length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)/length(data2$LENGTH))
  to_draw = data.frame()
  to_draw = cbind(tra, count)
  to_draw = cbind(to_draw, pa)
  to_draw
  to_draw = data.frame(to_draw)
  colnames(to_draw) = c("trait", "proportion", "Group")
  to_draw$proportion = as.numeric(to_draw$proportion)
  to_draw 
  ggplot(to_draw, aes(fill = as.factor(Group), x=trait, y=proportion )) + geom_col(position = position_dodge())  +  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("proportion") + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
  #ggsave("symmetric_bars.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
  
  
  # SYM NOT FOUND TESTS
  symet_not_found = data.frame("coloc" = c(length(data[data$LENGTH %% 3 == 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data[data$LENGTH %% 3 != 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                     "non_coloc" = c(length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data2[data2$LENGTH %% 3 != 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                     #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                     row.names = c("sym", "non_sym"),
                     stringsAsFactors = FALSE)
  symet_not_found
  test = fisher.test(symet_not_found)
  pval = c(pval, test$p.value)
  trait = c(trait, "Symm_not_found")
  pair = c(pair, "coloc vs non_coloc")
  conf_low = c(conf_low, test$conf.int[1])
  conf_top = c(conf_top, test$conf.int[2])
  stat = c(stat, test$estimate)
  test_type = c(test_type, "fisher")
  
  
  pa = c("coloc", "non_coloc", "coloc", "non_coloc")
  tra = c("Symmetric", "Symmetric",  "Non symmetric", "Non symmetric")
  count = c(length(data[data$LENGTH %% 3 == 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data[data$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data2[data2$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH), length(data[data$LENGTH %% 3 != 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data[data$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),  length(data2[data2$LENGTH %% 3 !=  0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data2[data2$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH))
  to_draw = data.frame()
  to_draw = cbind(tra, count)
  to_draw = cbind(to_draw, pa)
  to_draw
  to_draw = data.frame(to_draw)
  colnames(to_draw) = c("trait", "proportion", "Group")
  to_draw$proportion = as.numeric(to_draw$proportion)
  to_draw 
  ggplot(to_draw, aes(fill = as.factor(Group), x=trait, y=proportion )) + geom_col(position = position_dodge())  +  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("proportion") + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
  #ggsave("symmetric_not_found_bars.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
  
  
  
  # HELIX TESTS
  hel = data.frame("coloc" = c(length(data[data$HELIX == "True",]$HELIX), length(data[data$HELIX != "True",]$HELIX)),
                   "non_coloc" = c(length(data2[data2$HELIX == "True",]$HELIX), length(data2[data2$HELIX != "True",]$HELIX)),
                   #"non_sQTL" = c(length(data3[data3$HELIX ="True",]$HELIX), length(data3[data3$HELIX !"True",]$HELIX)),
                     row.names = c("helix", "no_helix"),
                     stringsAsFactors = FALSE)
  
  hel
  test = fisher.test(hel)
  pval = c(pval, test$p.value)
  trait = c(trait, "Helix")
  pair = c(pair, "coloc vs non_coloc")
  conf_low = c(conf_low, test$conf.int[1])
  conf_top = c(conf_top, test$conf.int[2])
  stat = c(stat, test$estimate)
  test_type = c(test_type, "fisher")
  
  
  # SHEETS TESTS
  she = data.frame("coloc" = c(length(data[data$SHEET == "True",]$SHEET), length(data[data$SHEET != "True",]$SHEET)),
                   "non_coloc" = c(length(data2[data2$SHEET == "True",]$SHEET), length(data2[data2$SHEET != "True",]$SHEET)),
                   #"non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                   row.names = c("sheet", "no_sheet"),
                   stringsAsFactors = FALSE)
  
  she
  test = fisher.test(she)
  pval = c(pval, test$p.value)
  trait = c(trait, "Sheet")
  pair = c(pair, "coloc vs non_coloc")
  conf_low = c(conf_low, test$conf.int[1])
  conf_top = c(conf_top, test$conf.int[2])
  stat = c(stat, test$estimate)
  test_type = c(test_type, "fisher")
  
  # TURNS TESTS
  turn = data.frame("coloc" = c(length(data[data$TURN == "True",]$TURN), length(data[data$TURN != "True",]$TURN)),
                   "non_coloc" = c(length(data2[data2$TURN == "True",]$TURN), length(data2[data2$TURN != "True",]$TURN)),
                   #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                   row.names = c("turn", "no_turn"),
                   stringsAsFactors = FALSE)
  
  turn
  test = fisher.test(turn)
  pval = c(pval, test$p.value)
  trait = c(trait, "Turn")
  pair = c(pair, "coloc vs non_coloc")
  conf_low = c(conf_low, test$conf.int[1])
  conf_top = c(conf_top, test$conf.int[2])
  stat = c(stat, test$estimate)
  test_type = c(test_type, "fisher")
  
  #TRANSMEMBRANE TEST
  
  trans= data.frame("coloc" = c(length(data[data$TRANSMEMBRANE == "True",]$LENGTH), length(data[data$TRANSMEMBRANE != "True",]$LENGTH)),
                               "non_coloc" = c(length(data2[data2$TRANSMEMBRANE == "True",]$LENGTH), length(data2[data2$TRANSMEMBRANE != "True",]$LENGTH)),
                               #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                               row.names = c("trans", "non_trans"),
                               stringsAsFactors = FALSE)
  trans
  test = fisher.test(trans)
  pval = c(pval, test$p.value)
  trait = c(trait, "Transmembrane")
  pair = c(pair, "coloc vs non_coloc")
  conf_low = c(conf_low, test$conf.int[1])
  conf_top = c(conf_top, test$conf.int[2])
  stat = c(stat, test$estimate)
  test_type = c(test_type, "fisher")
  
  # #DOMAINS ANALYSIS
  # ggplot() + geom_bar(data=data[data$DOMAIN != "['Disordered']" & data$DOMAIN != "[] ", ], aes(DOMAIN, fill="coloc"), stat = "count", alpha=0.3) + geom_bar(data=data2[data2$DOMAIN != "['Disordered']" & data2$DOMAIN != "[]", ], aes(DOMAIN, fill="non coloc"), stat = "count", alpha=0.3) + geom_bar(data=data3[data3$DOMAIN != "['Disordered']" & data3$DOMAIN != "Actin" & data3$DOMAIN != "[]", ], aes(DOMAIN, fill="non sQTL"), stat = "count", alpha=0.3)
  # table(data[data$DOMAIN != "['Disordered']" & data$DOMAIN != "[]", ]$DOMAIN)[table(data$DOMAIN) >= 2]
  # table(data2[data2$DOMAIN != "Disordered  " & data2$DOMAIN != "", ]$DOMAIN)[table(data2$DOMAIN) >= 2]
  # table(data3[data3$DOMAIN != "Disordered  " & data3$DOMAIN != "", ]$DOMAIN)[table(data3$DOMAIN) >= 2]
  
    
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
  to_draw
  colnames(to_draw) = c("trait", "pval", "pair", "conf_low", "conf_top", "statistics", "test_type")
  to_draw = data.frame(to_draw)
  to_draw$pval = -as.numeric(to_draw$pval)
  to_draw
  line = replicate(length(pval), 1.5)
  ggplot(to_draw, aes(fill = as.factor(pair), x=trait, y=pval )) + geom_col(position = position_dodge())  + geom_hline(yintercept = line) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("-log10(p_value)") + labs(fill='Pairs') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
  file = paste(c("statistical_summary_without_correction_cross_val_", as.character(i),".png"), collapse="")
  ggsave(file, path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
  
  to_draw$trait = factor(to_draw$trait, levels = unique(trait))
  
  to_draw$conf_low = as.numeric(to_draw$conf_low)
  to_draw$conf_top = as.numeric(to_draw$conf_top)
  
  to_draw$alpha = ifelse((to_draw$test_type == "ks"  & ((to_draw$conf_low <= 0 & to_draw$conf_top >= 0) | 
                                                          (to_draw$conf_low >= 0 & to_draw$conf_top <= 0))) | 
                           (to_draw$test_type == "fisher"  & ((log(to_draw$conf_low) <= 0 & log(to_draw$conf_top >= 0)) | 
                                                (log(to_draw$conf_low) >= 0 & log(to_draw$conf_top) <= 0))), yes = 0.6, 1.0)
  to_draw[to_draw$test_type == "fisher",]$alpha = ifelse((log(to_draw[to_draw$test_type == "fisher",]$conf_low) <= 0 & log(to_draw[to_draw$test_type == "fisher",]$conf_top) >= 0) | (log(to_draw[to_draw$test_type == "fisher",]$conf_low) >= 0 & log(to_draw[to_draw$test_type == "fisher",]$conf_top) <= 0), yes = 0.2, 1.0)
  
  
  
  unique(trait)
  to_draw$trait = factor(to_draw$trait, levels = unique(trait))
  line1 = replicate(length(pval), 0.0)
  View(to_draw)
  
  ggplot(to_draw[to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)),  alpha = alpha)) +
    geom_point(size = 4, position=position_dodge2(width=0.9)) +
    geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low))), 
                  position=position_dodge2(width=0.9)) +  
    ylab("Enrichment in non coloc                                                                 Enrichment in coloc\n log(odd_ratio)") + 
    coord_flip() + geom_hline(yintercept = line1, color="red") + 
    theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "both")),
          axis.title.x = element_text(angle = 0)) + gtex_v8_figure_theme() + guides(alpha = FALSE)
  file = paste(c("statistical_summary_fisher_cross_val_", as.character(i),".png"), collapse="")
  ggsave('statistical_summary_fisher_coloc.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
  
  line1 = replicate(length(pval), 0.0)
  ggplot(to_draw[to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics),  alpha = alpha)) +
    geom_point(size = 4, position=position_dodge2(width=0.9)) +
    geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low)), position=position_dodge2(width=0.9)) + 
    ylab("Enrichment in non coloc                                                                  Enrichment in coloc\n m1 - m2") + 
    coord_flip() + geom_hline(yintercept = line1, color="red") + 
    theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "both")),
          axis.title.x = element_text(angle = 0)) + gtex_v8_figure_theme()  + guides(alpha = FALSE)
  file = paste(c("statistical_summary_ks_cross_val_", as.character(i),".png"), collapse="")
  ggsave('statistical_summary_ks_coloc.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
  
  ggplot() + geom_density(data=data[data$DIP_P <= 0.05, ], aes(DIP, fill="coloc"), alpha = 0.2) + geom_density(data=data2[data2$DIP_P <= 0.05, ], aes(DIP, fill="non_coloc"), alpha = 0.2) + gtex_v8_figure_theme()
  ggsave('DIP_stat_coloc.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
  
  ggplot() + geom_histogram(data=data, aes(DIP_P, fill="coloc"), alpha = 0.2) + gtex_v8_figure_theme()
  ggsave('DIP_stat_coloc.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
  
  to_plot = с(data$DIP, data2$DIP)
  to_plot = data.frame(to_plot)
  to_plot = to_plot[to_plot$to_plot < 0.5, ]
  ggplot() + stat_qq(data=data, aes(sample=DIP, colour="coloc")) + 
    stat_qq(data=data2[data2$DIP < 0.5, ], aes(sample=DIP, colour="non_coloc")) + 
    stat_qq_line(data=data, aes(sample=DIP, colour="coloc")) + 
    stat_qq_line(data=data2[data2$DIP < 0.5, ], aes(sample=DIP, colour="non_coloc")) +
    gtex_v8_figure_theme()
  ggsave('DIP_stat_coloc_qq_plot.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
  
  #ggsave('DIP_p_coloc.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
  
  # s_ks = cbind(s_ks, to_draw[to_draw$test_type == "ks", ]$statistics)
  # c_l_ks = data.frame(to_draw[to_draw$test_type == "ks", ]$conf_low)
  # c_t_ks = data.frame(to_draw[to_draw$test_type == "ks", ]$conf_top)
  # s_f = data.frame(to_draw[to_draw$test_type == "fisher", ]$statistics)
  # c_l_f = data.frame(to_draw[to_draw$test_type == "fisher", ]$conf_low)
  # c_t_f = data.frame(to_draw[to_draw$test_type == "fisher", ]$conf_top)
#}





# View(to_draw)
# 
# 
# View(to_draw[to_draw$test_type == "ks", ])
# write.csv(to_draw[to_draw$test_type == "ks", ], "enrichent_analysis_ks_stat.csv", quote=F, row.names = F)

  
  
  
  
# 
# 
# 
# pre_data_1 = data[, 19:28]
# pre_data_1$LENGTH = data$LENGTH
# #pre_data_1$LENGTH.ALIGN = data$LENGTH.ALIGN
# pre_data_1$STRUC = as.numeric(data$GeneID %in% data_struc$GeneID)
# pre_data_1$UNSTRUC = as.numeric(data$GeneID %in% data_unstruc$GeneID)
# pre_data_1$TRANSMEMBRANE = as.numeric(data$TRANSMEMBRANE == "True")
# pre_data_1$HELIX = as.numeric(data$HELIX == "True")
# pre_data_1$SHEET = as.numeric(data$SHEET == "True")
# pre_data_1$TURN = as.numeric(data$TURN == "True")
# pre_data_1$ASN = data$ASN../data$LENGTH *100
# pre_data_1$CYS = data$CYS../data$LENGTH *100
# pre_data_1$SIGNAL = as.numeric(data$SIGNAL != "[]" & data$SIGNAL != "['']" & data$SIGNAL != "['', '']")
# pre_data_1$GROUP = replicate(nrow(pre_data_1), "coloc")
# 
# 
# pre_data_2 = data2[, 24:33]
# pre_data_2$LENGTH = data2$LENGTH
# #pre_data_2$LENGTH.ALIGN = data2$LENGTH.ALIGN
# pre_data_2$STRUC = as.numeric(data2$gene_id %in% data2_struc$gene_id)
# pre_data_2$UNSTRUC = as.numeric(data2$gene_id %in% data2_unstruc$gene_id)
# pre_data_2$TRANSMEMBRANE = as.numeric(data2$TRANSMEMBRANE == "True")
# pre_data_2$HELIX = as.numeric(data2$HELIX == "True")
# pre_data_2$SHEET = as.numeric(data2$SHEET == "True")
# pre_data_2$TURN = as.numeric(data2$TURN == "True")
# pre_data_2$ASN = data2$ASN../data2$LENGTH *100
# pre_data_2$CYS = data2$CYS../data2$LENGTH *100
# pre_data_2$SIGNAL = as.numeric(data2$SIGNAL != "[]" & data2$SIGNAL != "['']" & data2$SIGNAL != "['', '']")
# pre_data_2$GROUP = replicate(nrow(pre_data_2), "non_coloc")
# pre_data_2 = sample_n(pre_data_2, 200)
# 
# pre_data_3 = data3[, 21:30]
# pre_data_3$LENGTH = data3$LENGTH
# #pre_data_3$LENGTH.ALIGN = data3$LENGTH.ALIGN
# pre_data_3$STRUC = as.numeric(data3$GENE.NAME %in% data3_struc$GENE.NAME)
# pre_data_3$UNSTRUC = as.numeric(data3$GENE.NAME %in% data3_unstruc$GENE.NAME)
# pre_data_3$TRANSMEMBRANE = as.numeric(data3$TRANSMEMBRANE == "True")
# pre_data_3$HELIX = as.numeric(data3$HELIX == "True")
# pre_data_3$SHEET = as.numeric(data3$SHEET == "True")
# pre_data_3$TURN = as.numeric(data3$TURN == "True")
# pre_data_3$ASN = data3$ASN../data3$LENGTH *100
# pre_data_3$CYS = data3$CYS../data3$LENGTH *100
# pre_data_3$SIGNAL = as.numeric(data3$SIGNAL != "[]" & data3$SIGNAL != "['']" & data3$SIGNAL != "['', '']")
# pre_data_3$GROUP = replicate(nrow(pre_data_3), "non_sQTL")
# pre_data_3 = sample_n(pre_data_3, 200)
# 
# data_for_pca = rbind(pre_data_1, pre_data_2, pre_data_3)
# data_for_pca = na.omit(data_for_pca)
# data_for_pca.pca = prcomp(data_for_pca[, 1:20], center = T, scale. = T)
# data_for_pca.pca
# View(data_for_pca)
# 
# library(ggfortify)
# autoplot(data_for_pca.pca, data = data_for_pca, colour = 'STRUC')
# autoplot(data_for_pca.pca, data = data_for_pca, colour = 'UNSTRUC')
# autoplot(data_for_pca.pca, data = data_for_pca, colour = 'GROUP')
# autoplot(data_for_pca.pca, data = data_for_pca, colour = 'TRANSMEMBRANE')
# autoplot(data_for_pca.pca, data = data_for_pca, colour = 'SIGNAL')
# 
# # wss <- (nrow(data_for_pca[, 1:10])-1)*sum(apply(data_for_pca[, 1:10],2,var))
# # for (i in 2:15) wss[i] <- sum(kmeans(data_for_pca[, 1:10],
# #                                      centers=i)$withinss)
# # plot(1:15, wss, type="b", xlab="Number of Clusters",
# #      ylab="Within groups sum of squares")
# comp <- data.frame(data_for_pca.pca$x[,1:4])
# k <- kmeans(comp, 2, nstart=25, iter.max=1000)
# library(RColorBrewer)
# library(scales)
# palette(alpha(brewer.pal(9,'Set1'), 0.5))
# plot(comp, col=k$clust, pch=16)
# comp$STRUC = data_for_pca$STRUC
# comp$UNSTRUC = data_for_pca$UNSTRUC
# comp$GROUP = data_for_pca$GROUP
# ggplot(comp, aes(x=PC1, y=PC2)) +
#         geom_point(aes(color=as.factor(k$cluster), shape=as.factor(GROUP)), size=2.5) + scale_shape_manual(values=c(3, 25, 16))
# #ggsave("PCA_try.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
# 
# 
# ggplot(comp, aes(x=PC1, y=PC2)) +
#         geom_point(aes(size=as.factor(UNSTRUC), color=as.factor(k$cluster), shape=as.factor(UNSTRUC))) + scale_shape_manual(values=c(3, 13, 17))
# sort(table(k$clust))
# clust <- names(sort(table(k$clust)))
# # First cluster
# table(data_for_pca[k$clust==clust[1], ]$GROUP)
# # Second Cluster
# table(data_for_pca[k$clust==clust[2], ]$GROUP)
# table(data_for_pca[k$clust==clust[3], ]$GROUP)
# # table(data_for_pca[k$clust==clust[4], ]$STRUC)
# table(data_for_pca[k$clust==clust[1], ]$UNSTRUC)
# table(data_for_pca[k$clust==clust[2], ]$UNSTRUC)
# table(data_for_pca[k$clust==clust[3], ]$UNSTRUC)
# 
# table(data_for_pca[k$clust==clust[1], ]$STRUC)
# table(data_for_pca[k$clust==clust[2], ]$STRUC)
# table(data_for_pca[k$clust==clust[3], ]$STRUC)
# # table(data_for_pca[k$clust==clust[4], ]$UNSTRUC)
# 
# mean1 = numeric(0)
# max1 = numeric(0)
# for (i in 1:(ncol(pre_data_1)-1)){
#         mean1 = c(mean1, mean(pre_data_1[,i], na.rm = T))
#         max1 = c(max1, max(pre_data_1[,i], na.rm = T))   
# }
# mean1
# 
# mean2 = numeric(0)
# max2 = numeric(0)
# for (i in 1:(ncol(pre_data_2)-1)){
#         mean2 = c(mean2, mean(pre_data_2[,i], na.rm = T))
#         max2 = c(max2, max(pre_data_2[,i], na.rm = T))   
# }
# 
# 
# mean3 = numeric(0)
# max3 = numeric(0)
# for (i in 1:(ncol(pre_data_3)-1)){
#         mean3 = c(mean3, mean(pre_data_3[,i], na.rm = T))
#         max3 = c(max3, max(pre_data_3[,i], na.rm = T))   
# }
# mean1 = mean1/max1
# mean2 = mean2/max2
# mean3 = mean3/max3
# data_mean = rbind(mean1, mean2, mean3)
# data_mean
# colnames(data_mean) = colnames(pre_data_1)[1:(ncol(pre_data_1)-1)]
# data_mean
# ggplot(data_mean, aes(x=replicate(3, "MIN"), y=MIN)) + geom_col(fill=c(1,2,3), position = position_dodge())  +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("-log10(p_value)") + labs(fill='Pairs') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
# #ggsave("statistical_summary_with_correction.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
# 
