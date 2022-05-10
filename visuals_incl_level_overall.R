setwd("~/splicing_project/")
###___________________________________________________CHANGED______________________________________________________
data_full = read.csv("Data/combined_sQTL_data.csv")
data_full_1 = read.csv("Data/output_non_sQTL_full.csv")
data_full_1 %>% select(intersect(colnames(data_full), colnames(data_full_1))) -> data_full_1
View(data_full_1)
ncol(data_full)
ncol(data_full_1)
data_full = data_full[!is.na(data_full$anc_allele_freq), ]
data_full %>% select(intersect(colnames(data_full), colnames(data_full_1))) -> data_full
data_full = rbind(data_full, data_full_1)
data2 = data_full[data_full$mean_01_psi < 0.5, ]

data = data_full[data_full$mean_01_psi > 0.5, ]
nrow(data)
nrow(data2)


View(data)

data = unique.data.frame(data)
data2 = unique.data.frame(data2)
# View(data)
nrow(data)
nrow(data2)
set.seed(157)
library(ggplot2)
library(stats)
library(dplyr)
library(entropy)
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

inv_norm <- function(x) qnorm((rank(x, na.last='keep', ties.method = "random") - 0.5)/sum(!is.na(x)))

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
        d = x[x[,1] %in% sample(x[,1], replace = TRUE),]
        colnames(d) = c("totalData", "V2")
        return(d)
}


ks_stat = function(x, y){
        #x = data[grepl("['Disordered']", data$DOMAIN, fixed = TRUE),]$MIN_pLLDT
        #x = x[!is.na(x)]
        #y = data[!grepl("['Disordered']", data$DOMAIN, fixed = TRUE),]$MIN_pLLDT
        #y = y[!is.na(y)]
        pval =  ks.test(x, y)$p.value
        totalData = c(x, y)
        totalData = cbind(totalData, c(replicate(length(x), "sQTL"), replicate(length(y), "non_sQTL")))
        #View(totalData)
        totalData[,1] = inv_norm(as.numeric(totalData[,1]))
        totalData = data.frame(totalData)
        totalData = totalData[!is.na(totalData[,1]), ]
        totalData[,1] = as.numeric(totalData[,1])
        print(median(totalData[totalData$V2 == "sQTL", 'totalData']) - median(totalData[totalData$V2 == "non_sQTL", 'totalData']))
        
        totalBoot = my_bootstrap(totalData, R=10000, medianDiff, select_samples)
        #totalBoot = boot(totalData, statistic = medianDiff, R = 10000)
        totalBootCI = quantile(totalBoot, probs=c(0.05, 0.95, 0.0, 1.0, 0.5))
        print(totalBootCI)
        conf_low = totalBootCI[1]
        conf_top = totalBootCI[2]
        mi = totalBootCI[3]
        ma = totalBootCI[4]
        stat = totalBootCI[5]
        #stat = median(totalData[totalData$V2 == "sQTL", 'totalData'], na.rm = T) - median(totalData[totalData$V2 == "non_sQTL", 'totalData'], na.rm = T)
        print(c(pval, stat, conf_low, conf_top, mi, ma))
        return(c(pval, stat, conf_low, conf_top, mi, ma))
}


data = data[!flag_outliers(data$LENGTH) & data$LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$LENGTH) & data2$LENGTH >= 15,]
nrow(data)
nrow(data2)
# View(data)
# View(data[!flag_outliers(data$LENGTH) & data$LENGTH!=0 ,])
#dev.off()






data_unstruc = data[data$Q2 > 40.0 & data$Q2_pLLDT < 50.0, ]
data_unstruc = data_unstruc[!is.na(data_unstruc$GENE.NAME),] # 84/571 = 14,7% 67 signal
nrow(data_unstruc)
data2_unstruc = data2[data2$Q2 > 40.0 & data2$Q2_pLLDT < 50.0, ]
data2_unstruc = data2_unstruc[!is.na(data2_unstruc$GENE.NAME),] # 106/686 = 15,5% 68 signal
nrow(data2_unstruc)





# View(data_unstruc)
# View(data2_unstruc)
# View(data3_unstruc)
#write.csv(data_unstruc, "Data/output_data_higher_in_unstruc.csv", row.names = F)


data_struc = data[data$Q2 < 25.0 & data$Q2_pLLDT >= 90.0, ]
data_struc = data_struc[!is.na(data_struc$GENE.NAME),] # 105/571 = 18,4% 69 signals
nrow(data_struc)
# View(data_struc)
data2_struc = data2[data2$Q2 < 25.0 & data2$Q2_pLLDT >= 90.0, ]
data2_struc = data2_struc[!is.na(data2_struc$GENE.NAME),] # 108/686 = 15,7% 66 signals
nrow(data2_struc)
# View(data2_struc)
# View(data_struc)
# View(data2_struc)
# View(data3_struc)
#write.csv(data_struc, "Data/output_data_higher_in_struc.csv", row.names = F)

data_cons = data[data$Q2_pLLDT >= 90.0, ]
data_cons = data_cons[!is.na(data_cons$GENE.NAME),] # 230/571 = 40,1% 156 signals
nrow(data_cons)
data2_cons = data2[data2$Q2_pLLDT >= 90.0, ]
data2_cons = data2_cons[!is.na(data2_cons$GENE.NAME),] # 246/686 = 35,8% 162 signals
nrow(data2_cons)

# View(data_cons)
# View(data2_cons)
# View(data3_cons)
#write.csv(data_cons, "Data/output_data_higher_in_potent_cons.csv", row.names = F)


data_loop = data[data$TURN == "True" & data$HELIX == "False" & data$SHEET == "False", ]
data_loop = data_loop[!is.na(data_loop$GENE.NAME),] # 59/571 = 10,3 % 37 signal
nrow(data_loop)
data2_loop = data2[data2$TURN == "True" & data2$HELIX == "False" & data2$SHEET == "False", ]
data2_loop = data2_loop[!is.na(data2_loop$GENE.NAME),] # 54/686 = 7,9 % 32 signal
nrow(data2_loop)

# View(data_loop)
# View(data2_loop)
# View(data3_loop)
#write.csv(data_loop, "Data/output_data_higher_in_loop_cut.csv", row.names = F)

dev.off()
to_draw = data.frame()
round(nrow(data_unstruc) / nrow(data) * 100)
pa = replicate(round(nrow(data_unstruc) / nrow(data) * 100), "higher_in")
pa = c(pa, replicate(round(nrow(data2_unstruc) / nrow(data2) * 100), "lower_in"))
pa = c(pa, replicate(round(nrow(data_struc) / nrow(data) * 100), "higher_in"))
pa = c(pa, replicate(round(nrow(data2_struc) / nrow(data2) * 100), "lower_in"))
pa = c(pa, replicate(round(nrow(data_cons) / nrow(data) * 100), "higher_in"))
pa = c(pa, replicate(round(nrow(data2_cons) / nrow(data2) * 100), "lower_in"))
pa = c(pa, replicate(round(nrow(data_loop) / nrow(data) * 100), "higher_in"))
pa = c(pa, replicate(round(nrow(data2_loop) / nrow(data2) * 100), "lower_in"))
trait = replicate(round(nrow(data_unstruc) / nrow(data)*100) + round(nrow(data2_unstruc) / nrow(data2)*100), "Unstructured")
trait = c(trait, replicate(round(nrow(data_struc) / nrow(data)*100) + round(nrow(data2_struc) / nrow(data2)*100), "Structured"))
trait = c(trait, replicate(round(nrow(data_cons) / nrow(data)*100) + round(nrow(data2_cons) / nrow(data2)*100), "Conserved"))
trait = c(trait, replicate(round(nrow(data_loop) / nrow(data)*100) + round(nrow(data2_loop) / nrow(data2)*100), "Loop"))
print(length(pa))
print(length(trait))
to_draw = cbind(trait, pa)
to_draw
colnames(to_draw) = c("trait", "pair")
to_draw = data.frame(to_draw)
to_draw
ggplot(to_draw, aes(fill = as.factor(pa), x=trait)) + geom_bar(position = position_dodge())  +theme(axis.text.x = element_text( vjust = 0.5, hjust=1, size = 10), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8))  + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5))) + ylab("proportion")
#ggsave("structure_bars_der_all.png", height = 3.11, width = 7.92,path = "Data/visuals/", device='png', dpi=700)
to_draw
table(to_draw)[1, ]
chisq.test(table(to_draw)[1, ])
chisq.test(table(to_draw)[2, ])
chisq.test(table(to_draw)[3, ])
chisq.test(table(to_draw)[4, ])
pval = c()
pair = c()
trait = c()
conf_low = c()
conf_top = c()
stat = c()
test_type = c()
mi = c()
ma = c()
# 
# #HYDROPATHICITY TESTS
# trait = c(trait, "Hydropathicity")
# pair = c(pair, "higher_in vs lower_in")
# test_type = c(test_type, "ks")
# out = ks_stat(data[!flag_outliers(data$HYDROPATHICITY),]$HYDROPATHICITY,
#               data2[!flag_outliers(data2$HYDROPATHICITY),]$HYDROPATHICITY)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# mi = c(mi, out[5])
# ma = c(ma, out[6])

#LENGTH TESTS

trait = c(trait, "Length")
pair = c(pair, "higher_in vs lower_in")
test_type = c(test_type, "ks")
out = ks_stat(data$LENGTH, data2$LENGTH)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])
mi = c(mi, out[5])
ma = c(ma, out[6])

#ASN... TESTS

trait = c(trait, "% ASN")
pair = c(pair, "higher_in vs lower_in")
test_type = c(test_type, "ks")
out = ks_stat(data[!flag_outliers(data$ASN../data$LENGTH *100),]$ASN../data$LENGTH*100, 
              data2[!flag_outliers(data2$ASN../data2$LENGTH *100),]$ASN../data2$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])
mi = c(mi, out[5])
ma = c(ma, out[6])

#CYS.. TESTS

trait = c(trait, "% CYS")
pair = c(pair, "higher_in vs lower_in")
test_type = c(test_type, "ks")
out = ks_stat(data[!flag_outliers(data$CYS../data$LENGTH *100),]$CYS../data$LENGTH*100, 
              data2[!flag_outliers(data2$CYS../data2$LENGTH *100),]$CYS../data2$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])
mi = c(mi, out[5])
ma = c(ma, out[6])

# SYM TESTS

trait = c(trait, "Symm_dist")
pair = c(pair, "higher_in vs lower_in")
test_type = c(test_type, "ks")
out = ks_stat(data[!flag_outliers(data$LENGTH %% 3),]$LENGTH %% 3, 
              data2[!flag_outliers(data2$LENGTH %% 3),]$LENGTH %% 3.)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])
mi = c(mi, out[5])
ma = c(ma, out[6])
# 
# #MIN TESTS
# 
# trait = c(trait, "MIN_RSA")
# pair = c(pair, "higher_in vs lower_in")
# test_type = c(test_type, "ks")
# out = ks_stat(data[!flag_outliers(data$MIN),]$MIN,
#               data2[!flag_outliers(data2$MIN),]$MIN)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# mi = c(mi, out[5])
# ma = c(ma, out[6])

#Q1 TESTS

trait = c(trait, "Q1_RSA")
pair = c(pair, "higher_in vs lower_in")
test_type = c(test_type, "ks")
out = ks_stat(data[!flag_outliers(data$Q1),]$Q1,
              data2[!flag_outliers(data2$Q1),]$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])
mi = c(mi, out[5])
ma = c(ma, out[6])

# #Q2 TESTS
# 
# trait = c(trait, "Q2_RSA")
# pair = c(pair, "higher_in vs lower_in")
# test_type = c(test_type, "ks")
# out = ks_stat(data[!flag_outliers(data$Q2),]$Q2, 
#               data2[!flag_outliers(data2$Q2),]$Q2)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# mi = c(mi, out[5])
# ma = c(ma, out[6])
# 
# #Q3 TESTS
# 
# trait = c(trait, "Q3_RSA")
# pair = c(pair, "higher_in vs lower_in")
# test_type = c(test_type, "ks")
# out = ks_stat(data[!flag_outliers(data$Q3),]$Q3, 
#               data2[!flag_outliers(data2$Q3),]$Q3)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# mi = c(mi, out[5])
# ma = c(ma, out[6])
# 
# #MAX TESTS
# 
# trait = c(trait, "MAX_RSA")
# pair = c(pair, "higher_in vs lower_in")
# test_type = c(test_type, "ks")
# out = ks_stat(data[!flag_outliers(data$MAX),]$MAX, 
#               data2[!flag_outliers(data2$MAX),]$MAX)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# mi = c(mi, out[5])
# ma = c(ma, out[6])
# 
# #MIN_pLDDT TESTS
# 
# trait = c(trait, "MIN_pLDDT")
# pair = c(pair, "higher_in vs lower_in")
# test_type = c(test_type, "ks")
# out = ks_stat(data[!flag_outliers(data$MIN_pLLDT),]$MIN_pLLDT, 
#               data2[!flag_outliers(data2$MIN_pLLDT),]$MIN_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# mi = c(mi, out[5])
# ma = c(ma, out[6])
# 
# #Q1_pLDDT TESTS
# 
# trait = c(trait, "Q1_pLDDT")
# pair = c(pair, "higher_in vs lower_in")
# test_type = c(test_type, "ks")
# out = ks_stat(data[!flag_outliers(data$Q1_pLLDT),]$Q1_pLLDT, 
#               data2[!flag_outliers(data2$Q1_pLLDT),]$Q1_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# mi = c(mi, out[5])
# ma = c(ma, out[6])

# #Q2_pLDDT TESTS
# 
# trait = c(trait, "Q2_pLDDT")
# pair = c(pair, "higher_in vs lower_in")
# test_type = c(test_type, "ks")
# out = ks_stat(data[!flag_outliers(data$Q2_pLLDT),]$Q2_pLLDT, data2[!flag_outliers(data2$Q2_pLLDT),]$Q2_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# mi = c(mi, out[5])
# ma = c(ma, out[6])

#Q3_pLDDT TESTS

trait = c(trait, "Q3_pLDDT")
pair = c(pair, "higher_in vs lower_in")
test_type = c(test_type, "ks")
out = ks_stat(data[!flag_outliers(data$Q3_pLLDT),]$Q3_pLLDT, data2[!flag_outliers(data2$Q3_pLLDT),]$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])
mi = c(mi, out[5])
ma = c(ma, out[6])

# #MAX_pLDDT TESTS
# 
# trait = c(trait, "MAX_pLDDT")
# pair = c(pair, "higher_in vs lower_in")
# test_type = c(test_type, "ks")
# out = ks_stat(data[!flag_outliers(data$MAX_pLLDT),]$MAX_pLLDT, data2[!flag_outliers(data2$MAX_pLLDT),]$MAX_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# mi = c(mi, out[5])
# ma = c(ma, out[6])
# 
# #UNSTRUC TESTS
# unstruc = data.frame("higher_in" = c(nrow(data_unstruc), nrow(data) - nrow(data_unstruc)),
#                      "lower_in" = c(nrow(data2_unstruc), nrow(data2) - nrow(data2_unstruc)),
#                      #"lower_in" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                      row.names = c("unstruc", "other"),
#                      stringsAsFactors = FALSE)
# 
# unstruc
# test = fisher.test(unstruc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Unstructured")
# pair = c(pair, "higher_in vs lower_in")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# mi = c(mi, 0.0)
# ma = c(ma, 0.0)
# 
# # STRUC TESTS
# struc = data.frame("higher_in" = c(nrow(data_struc), nrow(data)-nrow(data_struc)),
#                    "lower_in" = c(nrow(data2_struc), nrow(data2)-nrow(data2_struc)),
#                    #"lower_in" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                    row.names = c("struc", "other"),
#                    stringsAsFactors = FALSE)
# 
# struc
# test = fisher.test(struc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Structured")
# pair = c(pair, "higher_in vs lower_in")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# mi = c(mi, 0.0)
# ma = c(ma, 0.0)
# 
# # CONSERVED TESTS
# cons = data.frame("higher_in" = c(nrow(data_cons), nrow(data) - nrow(data_cons)),
#                   "lower_in" = c(nrow(data2_cons), nrow(data2) - nrow(data2_cons)),
#                   #"lower_in" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                   row.names = c("cons", "other"),
#                   stringsAsFactors = FALSE)
# 
# cons
# test = fisher.test(cons)
# pval = c(pval, test$p.value)
# trait = c(trait, "Conserved")
# pair = c(pair, "higher_in vs lower_in")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# mi = c(mi, 0.0)
# ma = c(ma, 0.0)

# SYM TESTS
symet = data.frame("higher_in" = c(length(data[data$LENGTH %% 3 == 0,]$LENGTH), length(data[data$LENGTH %% 3 != 0,]$LENGTH)),
                   "lower_in" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
                   #"lower_in" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
pval = c(pval, test$p.value)
trait = c(trait, "Symmetric")
pair = c(pair, "higher_in vs lower_in")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")
mi = c(mi, 0.0)
ma = c(ma, 0.0)

pa = c("higher_in", "lower_in", "higher_in", "lower_in")
tra = c("Symmetric", "Symmetric", "Non symmetric", "Non symmetric")
count = c(length(data[data$LENGTH %% 3 == 0,]$LENGTH)/length(data$LENGTH),length(data2[data2$LENGTH %% 3 == 0,]$LENGTH)/length(data2$LENGTH), length(data[data$LENGTH %% 3 != 0,]$LENGTH)/length(data$LENGTH),  length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)/length(data2$LENGTH))
to_draw = cbind(tra, count)
to_draw = cbind(to_draw, pa)
to_draw
to_draw = data.frame(to_draw)
colnames(to_draw) = c("trait", "proportion", "Group")
to_draw$proportion = as.numeric(to_draw$proportion)
to_draw 
ggplot(to_draw, aes(fill = as.factor(Group), x=trait, y=proportion )) + geom_col(position = position_dodge())  +  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("proportion") + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
#ggsave("symmetric_bars_der_all.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)


# SYM NOT FOUND TESTS
symet_not_found = data.frame("higher_in" = c(length(data[data$LENGTH %% 3 == 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data[data$LENGTH %% 3 != 0 & data$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             "lower_in" = c(length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data2[data2$LENGTH %% 3 != 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             #"lower_in" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             row.names = c("sym", "non_sym"),
                             stringsAsFactors = FALSE)
symet_not_found
test = fisher.test(symet_not_found)
pval = c(pval, test$p.value)
trait = c(trait, "Symm_not_found")
pair = c(pair, "higher_in vs lower_in")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")
mi = c(mi, 0.0)
ma = c(ma, 0.0)

pa = c("higher_in", "lower_in", "higher_in", "lower_in")
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
ggplot(to_draw, aes(fill = as.factor(Group), x=trait, y=proportion )) + geom_col(position = position_dodge())  +  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("proportion") + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
#ggsave("symmetric_not_found_bars_higher_in.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)



# HELIX TESTS
hel = data.frame("higher_in" = c(length(data[data$HELIX == "True",]$HELIX), length(data[data$HELIX != "True",]$HELIX)),
                 "lower_in" = c(length(data2[data2$HELIX == "True",]$HELIX), length(data2[data2$HELIX != "True",]$HELIX)),
                 #"lower_in" = c(length(data3[data3$HELIX ="True",]$HELIX), length(data3[data3$HELIX !"True",]$HELIX)),
                 row.names = c("helix", "no_helix"),
                 stringsAsFactors = FALSE)

hel
test = fisher.test(hel)
pval = c(pval, test$p.value)
trait = c(trait, "Helix")
pair = c(pair, "higher_in vs lower_in")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")
mi = c(mi, 0.0)
ma = c(ma, 0.0)

# SHEETS TESTS
she = data.frame("higher_in" = c(length(data[data$SHEET == "True",]$SHEET), length(data[data$SHEET != "True",]$SHEET)),
                 "lower_in" = c(length(data2[data2$SHEET == "True",]$SHEET), length(data2[data2$SHEET != "True",]$SHEET)),
                 #"lower_in" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 row.names = c("sheet", "no_sheet"),
                 stringsAsFactors = FALSE)

she
test = fisher.test(she)
pval = c(pval, test$p.value)
trait = c(trait, "Sheet")
pair = c(pair, "higher_in vs lower_in")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")
mi = c(mi, 0.0)
ma = c(ma, 0.0)

# TURNS TESTS
turn = data.frame("higher_in" = c(length(data[data$TURN == "True",]$TURN), length(data[data$TURN != "True",]$TURN)),
                  "lower_in" = c(length(data2[data2$TURN == "True",]$TURN), length(data2[data2$TURN != "True",]$TURN)),
                  #"lower_in" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("turn", "no_turn"),
                  stringsAsFactors = FALSE)

turn
test = fisher.test(turn)
pval = c(pval, test$p.value)
trait = c(trait, "Turn")
pair = c(pair, "higher_in vs lower_in")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")
mi = c(mi, 0.0)
ma = c(ma, 0.0)

#TRANSMEMBRANE TEST

trans= data.frame("higher_in" = c(length(data[data$TRANSMEMBRANE == "True",]$LENGTH), length(data[data$TRANSMEMBRANE != "True",]$LENGTH)),
                  "lower_in" = c(length(data2[data2$TRANSMEMBRANE == "True",]$LENGTH), length(data2[data2$TRANSMEMBRANE != "True",]$LENGTH)),
                  #"lower_in" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                  row.names = c("trans", "non_trans"),
                  stringsAsFactors = FALSE)
trans
test = fisher.test(trans)
pval = c(pval, test$p.value)
trait = c(trait, "Transmembrane")
pair = c(pair, "higher_in vs lower_in")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")
mi = c(mi, 0.0)
ma = c(ma, 0.0)

features = c("DOMAIN", "TRANSMEMBRANE", "MOTIF", "TOPO_DOM", "ACT_SITE", "MOD_RES",
             "REGION", "REPEAT", "TRANSMEM", "NP_BIND", "DNA_BIND", "CROSSLNK", 
             "ZN_FING", "METAL", "SITE", "INTRAMEM", "LIPID")  
data[data == ""] = NA
data2[data2 == ""] = NA
is_domain_data = data.frame(is.na(data[,features]))
is_domain_data2 = data.frame(is.na(data2[, features]))
dim(is_domain_data)
is_domain_data$SUM = rowSums(is_domain_data)
is_domain_data2$SUM = rowSums(is_domain_data2)
is_domain_data$SUM = ifelse(is_domain_data$SUM == length(features), F, T)
is_domain_data2$SUM = ifelse(is_domain_data2$SUM == length(features), F, T)
dim(is_domain_data)

station = data.frame("higher_in" = c(sum(is_domain_data$SUM), nrow(is_domain_data) - sum(is_domain_data$SUM)),
                     "lower_in" = c(sum(is_domain_data2$SUM), nrow(is_domain_data2) - sum(is_domain_data2$SUM)),
                     row.names = c("stat", "non_stat"),
                     stringsAsFactors = FALSE)
print(station)
test = fisher.test(station)
est = fisher.test(station)
pval = c(pval, test$p.value)
trait = c(trait, "Domain")
pair = c(pair, "higher_in vs lower_in")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

mi = c(mi, 0.0)
ma = c(ma, 0.0)

# #DOMAINS ANALYSIS
# ggplot() + geom_bar(data=data[data$DOMAIN != "['Disordered']" & data$DOMAIN != "[] ", ], aes(DOMAIN, fill="higher_in"), stat = "count", alpha=0.3) + geom_bar(data=data2[data2$DOMAIN != "['Disordered']" & data2$DOMAIN != "[]", ], aes(DOMAIN, fill="non higher_in"), stat = "count", alpha=0.3)
# table(data[data$DOMAIN != "['Disordered']" & data$DOMAIN != "[]", ]$DOMAIN)[table(data$DOMAIN) >= 2]
# table(data2[data2$DOMAIN != "Disordered  " & data2$DOMAIN != "", ]$DOMAIN)[table(data2$DOMAIN) >= 2]



# for (i in seq(1, length(pval), by=3)){
#           print(pval[i:(i+2)])
#           pval[i:(i+2)] = p.adjust(pval[i:(i+2)], method="hochberg")
# }

# pval = p.adjust(pval, method="hochberg")

to_draw = data.frame()
to_draw = cbind(trait, log10(pval))
to_draw = cbind(to_draw, pair)
to_draw = cbind(to_draw, conf_low)
to_draw = cbind(to_draw, conf_top)
to_draw = cbind(to_draw, mi)
to_draw = cbind(to_draw, ma)
to_draw = cbind(to_draw, stat)
to_draw = cbind(to_draw, test_type)
to_draw
colnames(to_draw) = c("trait", "pval", "pair", "conf_low", "conf_top", "min", "max", "statistics", "test_type")
to_draw = data.frame(to_draw)
to_draw$pval = -as.numeric(to_draw$pval)
to_draw
line = replicate(length(pval), 1.5)
ggplot(to_draw, aes(fill = as.factor(pair), x=trait, y=pval )) + geom_col(position = position_dodge(), color="light blue", fill="light blue")  + 
        geom_hline(yintercept = line) +  
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), 
              legend.text = element_text(size=7), axis.title.y = element_text(size = 8), 
              axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("-log10(p_value)") + 
        labs(fill='Pairs') + guides(shape = guide_legend(override.aes = list(size = 0.5)), 
                                    color = guide_legend(override.aes = list(size = 0.5))) +
        gtex_v8_figure_theme()
ggsave("statistical_summary_without_correction_inclusion_groups.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)

to_draw$conf_low = as.numeric(to_draw$conf_low)
to_draw$conf_top = as.numeric(to_draw$conf_top)
to_draw$min= as.numeric(to_draw$min)
to_draw$max = as.numeric(to_draw$max)
to_draw$trait = factor(to_draw$trait, levels = trait)
to_draw$alpha = ifelse((to_draw$test_type == "ks"  & ((to_draw$conf_low <= 0 & to_draw$conf_top >= 0) | 
                                                              (to_draw$conf_low >= 0 & to_draw$conf_top <= 0))) | 
                               (to_draw$test_type == "fisher"  & ((log(to_draw$conf_low) <= 0 & log(to_draw$conf_top >= 0)) | 
                                                                          (log(to_draw$conf_low) >= 0 & log(to_draw$conf_top) <= 0))), yes = 0.4, 1.0)
to_draw[to_draw$test_type == "fisher",]$alpha = ifelse((log(to_draw[to_draw$test_type == "fisher",]$conf_low) <= 0 & log(to_draw[to_draw$test_type == "fisher",]$conf_top) >= 0) | (log(to_draw[to_draw$test_type == "fisher",]$conf_low) >= 0 & log(to_draw[to_draw$test_type == "fisher",]$conf_top) <= 0), yes = 0.4, 1.0)

line1 = replicate(length(pval), 0.0)

ggplot(to_draw[to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha=alpha)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low))))  + 
        ylab("Enrichment in lower inclusion exons                                                                  Enrichment in higher inclusion exons\n log(odd_ratio)") + 
        coord_flip() + geom_hline(yintercept = line1, color="red") + 
        theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "both")),
              axis.title.x = element_text(angle = 0)) + gtex_v8_figure_theme() + guides(alpha = FALSE)
ggsave("statistical_summary_fisher_inclusion_groups.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)

line1 = replicate(length(pval), 0.0)
#col = c('3', '4', '5', '6', 'light blue', 'orange', 'orange', 'orange', 'orange', 'orange', '2', '2', '2', '2', '2')
ggplot(to_draw[to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha=alpha)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) +  
        ylab("Enrichment in lower inclusion exons                                                                  Enrichment in higher inclusion exons\n m1-m2") + 
        coord_flip() + geom_hline(yintercept = line1, color="red") + 
        theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "both")),
              axis.title.x = element_text(angle = 0)) + gtex_v8_figure_theme() + guides(alpha = FALSE)
#ggsave("statistical_summary_with_correction_der_all.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggsave("statistical_summary_ks_inclusion_groups.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)



write.csv(to_draw, "Data/to_draw_inclusion_groups_overall.csv", quote=F, row.names = F)





ggplot() + geom_density(data=data[!flag_outliers(data$LENGTH %% 3),], aes(LENGTH %% 3, fill="higher_in"), alpha=0.3) + geom_density(data=data2[!flag_outliers(data2$LENGTH %% 3),], aes(LENGTH %% 3, fill="non higher_in"), alpha=0.3)


ggplot() + geom_density(data=data[!flag_outliers(data$MIN),], aes(MIN, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$MIN),], aes(MIN, fill="non higher_in"), alpha = 0.2) 

ggplot() + geom_density(data=data[!flag_outliers(data$Q1),], aes(Q1, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$Q1),], aes(Q1, fill="non higher_in"), alpha = 0.2) +
        annotate(geom="text", x=80, y=0.03, label=10^(-to_draw[to_draw$trait == "Q1_RSA", ]$pval))


ggplot() + geom_density(data=data[!flag_outliers(data$Q2),], aes(Q2, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$Q2),], aes(Q2, fill="non higher_in"), alpha = 0.2)

ggplot() + geom_density(data=data[!flag_outliers(data$Q3),], aes(Q3, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$Q3),], aes(Q3, fill="non higher_in"), alpha = 0.2) 

ggplot() + geom_density(data=data[!flag_outliers(data$MAX),], aes(MAX, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$MAX),], aes(MAX, fill="non higher_in"), alpha = 0.2)


ggplot() + geom_density(data=data[!flag_outliers(data$MIN_pLLDT),], aes(MIN_pLLDT, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$MIN_pLLDT),], aes(MIN_pLLDT, fill="non higher_in"), alpha = 0.2)

ggplot() + geom_density(data=data[!flag_outliers(data$Q1_pLLDT),], aes(Q1_pLLDT, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$Q1_pLLDT),], aes(Q1_pLLDT, fill="non higher_in"), alpha = 0.2)

ggplot() + geom_density(data=data[!flag_outliers(data$Q2_pLLDT),], aes(Q2_pLLDT, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$Q2_pLLDT),], aes(Q2_pLLDT, fill="non higher_in"), alpha = 0.2)

ggplot() + geom_density(data=data[!flag_outliers(data$Q3_pLLDT),], aes(Q3_pLLDT, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$Q3_pLLDT),], aes(Q3_pLLDT, fill="non higher_in"), alpha = 0.2) + 
        annotate(geom="text", x=20, y=0.15, label=10^(-to_draw[to_draw$trait == "Q3_pLDDT", ]$pval))

ggplot() + geom_density(data=data[!flag_outliers(data$MAX_pLLDT),], aes(MAX_pLLDT, fill="higher_in"), alpha = 0.2) + 
        geom_density(data=data2[!flag_outliers(data2$MAX_pLLDT),], aes(MAX_pLLDT, fill="non higher_in"), alpha = 0.2)









