setwd("~/splicing_project/")
data_full = read.csv("Data/combined_sQTL_data.csv")
View(data_full)
nrow(data_full)
data_full = data_full[!is.na(data_full$anc_allele_freq), ]
nrow(data_full)
anc_psi = c(data_full[data_full$anc_allele == data_full$ref_allele, ]$mean_00_psi, data_full[data_full$anc_allele == data_full$alt_allele, ]$mean_11_psi)
der_psi = c(data_full[data_full$anc_allele == data_full$ref_allele, ]$mean_11_psi, data_full[data_full$anc_allele == data_full$alt_allele, ]$mean_00_psi)
length(anc_psi)
length(der_psi)
data_full$anc_psi = anc_psi
data_full$der_psi = der_psi

data1 = data_full[data_full$mean_01_psi < 0.5 & data_full$anc_psi < data_full$mean_01_psi & data_full$mean_01_psi < data_full$der_psi, ]
nrow(data1)
data2 = data_full[data_full$mean_01_psi < 0.5 & data_full$anc_psi > data_full$mean_01_psi & data_full$mean_01_psi > data_full$der_psi, ]
nrow(data2)
data3 = data_full[data_full$mean_01_psi > 0.5 & data_full$anc_psi < data_full$mean_01_psi & data_full$mean_01_psi < data_full$der_psi, ]
nrow(data3)
data4 = data_full[data_full$mean_01_psi > 0.5 & data_full$anc_psi > data_full$mean_01_psi & data_full$mean_01_psi > data_full$der_psi, ]
nrow(data4)

library(ggplot2)
library(gridExtra)
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
        print(mean(totalData[totalData$V2 == "sQTL", 'totalData']) - mean(totalData[totalData$V2 == "non_sQTL", 'totalData']))
        
        totalBoot = my_bootstrap(totalData, R=1000, meanDiff, select_samples)
        #totalBoot = boot(totalData, statistic = meanDiff, R = 10000)
        totalBootCI = quantile(totalBoot, probs=c(0.05, 0.95))
        conf_low = totalBootCI[1]
        conf_top = totalBootCI[2]
        stat = mean(totalData[totalData$V2 == "sQTL", 'totalData']) - mean(totalData[totalData$V2 == "non_sQTL", 'totalData'])
        print(c(pval, stat, conf_low, conf_top))
        return(c(pval, stat, conf_low, conf_top))
}

# 
# data = data[!flag_outliers(data$LENGTH),]
# data2 = data2[!flag_outliers(data2$LENGTH),]
# nrow(data)
# nrow(data2)
# # View(data)
# # View(data[!flag_outliers(data$LENGTH) & data$LENGTH!=0 ,])
# dev.off()


data = data[!flag_outliers(data$LENGTH) & data$LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$LENGTH) & data2$LENGTH >= 15,]
data3 = data3[!flag_outliers(data3$LENGTH) & data3$LENGTH >= 15,]
data4 = data4[!flag_outliers(data4$LENGTH) & data4$LENGTH >= 15,]

data1$GROUP = 1
data2$GROUP = 2
data3$GROUP = 3
data4$GROUP = 4

data_full = rbind(data1, data2, data3, data4)



View(data_full)

plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(LENGTH, group=GROUP)) + coord_flip() + ggtitle("LENGTH")
plot2 = ggplot(data = data_full, aes(LENGTH, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("LENGTH")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_length_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)




plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(ASN../LENGTH *100, group=GROUP)) + coord_flip() + ggtitle("ASN../LENGTH *100")
plot2 = ggplot(data = data_full, aes(ASN../LENGTH *100, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("ASN../LENGTH *100")

plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_ASN_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)



plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(CYS../LENGTH *100, group=GROUP)) + coord_flip() + ggtitle("CYS../LENGTH *100")
plot2 = ggplot(data = data_full, aes(CYS../LENGTH *100, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("CYS../LENGTH *100")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_CYS_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)


plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(MIN, group=GROUP)) + coord_flip() + ggtitle("MIN")
plot2 = ggplot(data = data_full, aes(MIN, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("MIN")

plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_MIN_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)



plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(Q1, group=GROUP)) + coord_flip() + ggtitle("Q1")
plot2 = ggplot(data = data_full, aes(Q1, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("Q1")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_Q1_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)


plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(Q2, group=GROUP)) + coord_flip() + ggtitle("Q2")
plot2 = ggplot(data = data_full, aes(Q2, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("Q2")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_Q2_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)


plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(Q3, group=GROUP)) + coord_flip() + ggtitle("Q3")
plot2 = ggplot(data = data_full, aes(Q3, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("Q3")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_Q3_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)


plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(MAX, group=GROUP)) + coord_flip() + ggtitle("MAX")
plot2 = ggplot(data = data_full, aes(MAX, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("MAX")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_MAX_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)


plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(MIN_pLLDT, group=GROUP)) + coord_flip() + ggtitle("MIN_pLLDT")
plot2 = ggplot(data = data_full, aes(MIN_pLLDT, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("MIN_pLLDT")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_MIN_pLLDT_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)


plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(Q1_pLLDT, group=GROUP)) + coord_flip() + ggtitle("Q1_pLLDT")
plot2 = ggplot(data = data_full, aes(Q1_pLLDT, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("Q1_pLLDT")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_Q1_pLLDT_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)


plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(Q2_pLLDT, group=GROUP)) + coord_flip() + ggtitle("Q2_pLLDT")
plot2 = ggplot(data = data_full, aes(Q2_pLLDT, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("Q2_pLLDT")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_Q2_pLLDT_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)


plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(Q3_pLLDT, group=GROUP)) + coord_flip() + ggtitle("Q3_pLLDT")
plot2 = ggplot(data = data_full, aes(Q3_pLLDT, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("Q3_pLLDT")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_Q3_pLLDT_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)


plot1 = ggplot() + gtex_v8_figure_theme() + geom_boxplot(data=data_full, aes(MAX_pLLDT, group=GROUP)) + coord_flip() + ggtitle("MAX_pLLDT")
plot2 = ggplot(data = data_full, aes(MAX_pLLDT, color=as.factor(GROUP))) + gtex_v8_figure_theme() + geom_density(alpha = 0.2) + ggtitle("MAX_pLLDT")
plot = arrangeGrob(plot1, plot2, widths = c(1, 1),
                   layout_matrix = rbind(c(1, 2)))

ggsave("descriptive_analysis_MAX_pLLDT_psi.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)

data_full %>% 
        count(HELIX, GROUP) %>% 
        mutate(perc = n) -> data_full_2
data_full_2 = data_full_2[1:8, ]
data_full_2$perc[1] = data_full_2$perc[1]/(data_full_2$n[1] + data_full_2$n[5])
data_full_2$perc[2] = data_full_2$perc[2]/(data_full_2$n[2] + data_full_2$n[6])
data_full_2$perc[3] = data_full_2$perc[3]/(data_full_2$n[3] + data_full_2$n[7])
data_full_2$perc[4] = data_full_2$perc[4]/(data_full_2$n[4] + data_full_2$n[8])
data_full_2$perc[5] = data_full_2$perc[5]/(data_full_2$n[1] + data_full_2$n[5])
data_full_2$perc[6] = data_full_2$perc[6]/(data_full_2$n[2] + data_full_2$n[6])
data_full_2$perc[7] = data_full_2$perc[7]/(data_full_2$n[3] + data_full_2$n[7])
data_full_2$perc[8] = data_full_2$perc[8]/(data_full_2$n[4] + data_full_2$n[8])

ggplot(data_full_2, aes(x = HELIX, y=perc, fill=as.factor(GROUP))) + geom_bar(stat="identity", position=position_dodge()) + gtex_v8_figure_theme() + ggtitle("HELIX")
ggsave("descriptive_analysis_HELIX_psi.png", path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)

data_full %>% 
        count(SHEET, GROUP) %>% 
        mutate(perc = n ) -> data_full_2

data_full_2 = data_full_2[1:8, ]
data_full_2$perc[1] = data_full_2$perc[1]/(data_full_2$n[1] + data_full_2$n[5])
data_full_2$perc[2] = data_full_2$perc[2]/(data_full_2$n[2] + data_full_2$n[6])
data_full_2$perc[3] = data_full_2$perc[3]/(data_full_2$n[3] + data_full_2$n[7])
data_full_2$perc[4] = data_full_2$perc[4]/(data_full_2$n[4] + data_full_2$n[8])
data_full_2$perc[5] = data_full_2$perc[5]/(data_full_2$n[1] + data_full_2$n[5])
data_full_2$perc[6] = data_full_2$perc[6]/(data_full_2$n[2] + data_full_2$n[6])
data_full_2$perc[7] = data_full_2$perc[7]/(data_full_2$n[3] + data_full_2$n[7])
data_full_2$perc[8] = data_full_2$perc[8]/(data_full_2$n[4] + data_full_2$n[8])

ggplot(data_full_2, aes(x = SHEET, y = perc, fill=as.factor(GROUP))) + geom_bar(stat = "identity", position=position_dodge()) + gtex_v8_figure_theme() + ggtitle("SHEET")
ggsave("descriptive_analysis_SHEET_psi.png", path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)

data_full %>% 
        count(TURN, GROUP) %>% 
        mutate(perc = n ) -> data_full_2


data_full_2 = data_full_2[1:8, ]
data_full_2$perc[1] = data_full_2$perc[1]/(data_full_2$n[1] + data_full_2$n[5])
data_full_2$perc[2] = data_full_2$perc[2]/(data_full_2$n[2] + data_full_2$n[6])
data_full_2$perc[3] = data_full_2$perc[3]/(data_full_2$n[3] + data_full_2$n[7])
data_full_2$perc[4] = data_full_2$perc[4]/(data_full_2$n[4] + data_full_2$n[8])
data_full_2$perc[5] = data_full_2$perc[5]/(data_full_2$n[1] + data_full_2$n[5])
data_full_2$perc[6] = data_full_2$perc[6]/(data_full_2$n[2] + data_full_2$n[6])
data_full_2$perc[7] = data_full_2$perc[7]/(data_full_2$n[3] + data_full_2$n[7])
data_full_2$perc[8] = data_full_2$perc[8]/(data_full_2$n[4] + data_full_2$n[8])

ggplot(data_full_2, aes(x = TURN, y = perc, fill=as.factor(GROUP))) + geom_bar(stat = "identity", position=position_dodge()) + gtex_v8_figure_theme() + ggtitle("TURN")
ggsave("descriptive_analysis_TURN_psi.png", path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)



data1_unstruc = data1[data1$Q2 > 40.0 & data1$Q2_pLLDT < 50.0, ]
data1_unstruc = data1_unstruc[!is.na(data1_unstruc$GENE.NAME),] # 84/571 = 14,7% 67 signal
nrow(data1_unstruc)
data2_unstruc = data2[data2$Q2 > 40.0 & data2$Q2_pLLDT < 50.0, ]
data2_unstruc = data2_unstruc[!is.na(data2_unstruc$GENE.NAME),] # 106/686 = 15,5% 68 signal
nrow(data2_unstruc)
data3_unstruc = data3[data3$Q2 > 40.0 & data3$Q2_pLLDT < 50.0, ]
data3_unstruc = data3_unstruc[!is.na(data3_unstruc$GENE.NAME),] # 84/571 = 14,7% 67 signal
nrow(data3_unstruc)
data4_unstruc = data4[data4$Q2 > 40.0 & data4$Q2_pLLDT < 50.0, ]
data4_unstruc = data4_unstruc[!is.na(data4_unstruc$GENE.NAME),] # 106/686 = 15,5% 68 signal
nrow(data4_unstruc)




# View(data_unstruc)
# View(data2_unstruc)
# View(data3_unstruc)
#write.csv(data_unstruc, "Data/output_data_sQTL_unstruc.csv", row.names = F)


data1_struc = data1[data1$Q2 < 25.0 & data1$Q2_pLLDT >= 90.0, ]
data1_struc = data1_struc[!is.na(data1_struc$GENE.NAME),] # 105/571 = 18,4% 69 signals
nrow(data1_struc)

data2_struc = data2[data2$Q2 < 25.0 & data2$Q2_pLLDT >= 90.0, ]
data2_struc = data2_struc[!is.na(data2_struc$GENE.NAME),] # 108/686 = 15,7% 66 signals
nrow(data2_struc)

data3_struc = data3[data3$Q2 < 25.0 & data3$Q2_pLLDT >= 90.0, ]
data3_struc = data3_struc[!is.na(data3_struc$GENE.NAME),] # 105/571 = 18,4% 69 signals
nrow(data3_struc)

data4_struc = data4[data4$Q2 < 25.0 & data4$Q2_pLLDT >= 90.0, ]
data4_struc = data4_struc[!is.na(data4_struc$GENE.NAME),] # 108/686 = 15,7% 66 signals
nrow(data4_struc)
# View(data2_struc)
# View(data_struc)
# View(data2_struc)
# View(data3_struc)
#write.csv(data_struc, "Data/output_data_sQTL_struc.csv", row.names = F)

data1_cons = data1[data1$Q2_pLLDT >= 90.0, ]
data1_cons = data1_cons[!is.na(data1_cons$GENE.NAME),] # 230/571 = 40,1% 156 signals
nrow(data1_cons)
data2_cons = data2[data2$Q2_pLLDT >= 90.0, ]
data2_cons = data2_cons[!is.na(data2_cons$GENE.NAME),] # 246/686 = 35,8% 162 signals
nrow(data2_cons)
data3_cons = data3[data3$Q2_pLLDT >= 90.0, ]
data3_cons = data3_cons[!is.na(data3_cons$GENE.NAME),] # 230/571 = 40,1% 156 signals
nrow(data3_cons)
data4_cons = data4[data4$Q2_pLLDT >= 90.0, ]
data4_cons = data4_cons[!is.na(data4_cons$GENE.NAME),] # 246/686 = 35,8% 162 signals
nrow(data4_cons)

# View(data_cons)
# View(data2_cons)
# View(data3_cons)
#write.csv(data_cons, "Data/output_data_sQTL_potent_cons.csv", row.names = F)


data1_loop = data1[data1$TURN == "True" & data1$HELIX == "False" & data1$SHEET == "False", ]
data1_loop = data1_loop[!is.na(data1_loop$GENE.NAME),] # 59/571 = 10,3 % 37 signal
nrow(data1_loop)
data2_loop = data2[data2$TURN == "True" & data2$HELIX == "False" & data2$SHEET == "False", ]
data2_loop = data2_loop[!is.na(data2_loop$GENE.NAME),] # 54/686 = 7,9 % 32 signal
nrow(data2_loop)
data3_loop = data3[data3$TURN == "True" & data3$HELIX == "False" & data3$SHEET == "False", ]
data3_loop = data3_loop[!is.na(data3_loop$GENE.NAME),] # 59/571 = 10,3 % 37 signal
nrow(data3_loop)
data4_loop = data4[data4$TURN == "True" & data4$HELIX == "False" & data4$SHEET == "False", ]
data4_loop = data4_loop[!is.na(data4_loop$GENE.NAME),] # 54/686 = 7,9 % 32 signal
nrow(data4_loop)
# View(data_loop)
# View(data2_loop)
# View(data3_loop)
#write.csv(data_loop, "Data/output_data_sQTL_loop_cut.csv", row.names = F)

dev.off()
to_draw = data.frame()

pa = replicate(round(nrow(data1_unstruc) / nrow(data1) * 100), "1")
pa = c(pa, replicate(round(nrow(data2_unstruc) / nrow(data2) * 100), "2"))
pa = c(pa, replicate(round(nrow(data3_unstruc) / nrow(data3) * 100), "3"))
pa = c(pa, replicate(round(nrow(data4_unstruc) / nrow(data4) * 100), "4"))
length(pa)
pa = c(pa, replicate(round(nrow(data1_struc) / nrow(data1) * 100), "1"))
pa = c(pa, replicate(round(nrow(data2_struc) / nrow(data2) * 100), "2"))
pa = c(pa, replicate(round(nrow(data3_struc) / nrow(data3) * 100), "3"))
pa = c(pa, replicate(round(nrow(data4_struc) / nrow(data4) * 100), "4"))

pa = c(pa, replicate(round(nrow(data1_cons) / nrow(data1) * 100), "1"))
pa = c(pa, replicate(round(nrow(data2_cons) / nrow(data2) * 100), "2"))
pa = c(pa, replicate(round(nrow(data3_cons) / nrow(data3) * 100), "3"))
pa = c(pa, replicate(round(nrow(data4_cons) / nrow(data4) * 100), "4"))

pa = c(pa, replicate(round(nrow(data1_loop) / nrow(data1) * 100), "1"))
pa = c(pa, replicate(round(nrow(data2_loop) / nrow(data2) * 100), "2"))
pa = c(pa, replicate(round(nrow(data3_loop) / nrow(data3) * 100), "3"))
pa = c(pa, replicate(round(nrow(data4_loop) / nrow(data4) * 100), "4"))

trait = replicate(round(nrow(data1_unstruc) / nrow(data1)*100) + round(nrow(data2_unstruc) / nrow(data2)*100) + round(nrow(data3_unstruc) / nrow(data3)*100) + round(nrow(data4_unstruc) / nrow(data4)*100), "Unstructured")
trait = c(trait, replicate(round(nrow(data1_struc) / nrow(data1)*100) + round(nrow(data2_struc) / nrow(data2)*100) + round(nrow(data3_struc) / nrow(data3)*100) + round(nrow(data4_struc) / nrow(data4)*100), "Structured"))
trait = c(trait, replicate(round(nrow(data1_cons) / nrow(data1)*100) + round(nrow(data2_cons) / nrow(data2)*100) + round(nrow(data3_cons) / nrow(data3)*100) + round(nrow(data4_cons) / nrow(data4)*100), "Conserved"))
trait = c(trait, replicate(round(nrow(data1_loop) / nrow(data1)*100) + round(nrow(data2_loop) / nrow(data2)*100) + round(nrow(data3_loop) / nrow(data3)*100) + round(nrow(data4_loop) / nrow(data4)*100), "Loop"))

print(length(pa))
print(length(trait))
to_draw = cbind(trait, pa)
to_draw
colnames(to_draw) = c("trait", "pair")
to_draw = data.frame(to_draw)
to_draw
ggplot(to_draw, aes(fill = as.factor(pa), x=trait)) + gtex_v8_figure_theme() + geom_bar(position = position_dodge())  +theme(axis.text.x = element_text( vjust = 0.5, hjust=1, size = 10), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8))  + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5))) + ylab("proportion")
ggsave("structure_bars_among_groups.png", height = 3.11, width = 7.92,path = "Data/visuals/", device='png', dpi=700)
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
# 
# 
# #HYDROPATHICITY TESTS
# trait = c(trait, "Hydropathicity")
# pair = c(pair, "1 vs 2")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$HYDROPATHICITY),]$HYDROPATHICITY,
#               data2[!flag_outliers(data2$HYDROPATHICITY),]$HYDROPATHICITY)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Hydropathicity")
# pair = c(pair, "1 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$HYDROPATHICITY),]$HYDROPATHICITY,
#               data3[!flag_outliers(data3$HYDROPATHICITY),]$HYDROPATHICITY)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "Hydropathicity")
# pair = c(pair, "1 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$HYDROPATHICITY),]$HYDROPATHICITY,
#               data4[!flag_outliers(data4$HYDROPATHICITY),]$HYDROPATHICITY)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Hydropathicity")
# pair = c(pair, "2 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$HYDROPATHICITY),]$HYDROPATHICITY,
#               data3[!flag_outliers(data3$HYDROPATHICITY),]$HYDROPATHICITY)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Hydropathicity")
# pair = c(pair, "2 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$HYDROPATHICITY),]$HYDROPATHICITY,
#               data4[!flag_outliers(data4$HYDROPATHICITY),]$HYDROPATHICITY)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "Hydropathicity")
# pair = c(pair, "3 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data3[!flag_outliers(data3$HYDROPATHICITY),]$HYDROPATHICITY,
#               data4[!flag_outliers(data4$HYDROPATHICITY),]$HYDROPATHICITY)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 

#LENGTH TESTS

trait = c(trait, "Length")
pair = c(pair, "1 vs 2")
test_type = c(test_type, "ks")
out = ks_stat(data1$LENGTH, data2$LENGTH)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Length")
pair = c(pair, "1 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data1$LENGTH, data3$LENGTH)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Length")
pair = c(pair, "1 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data1$LENGTH, data4$LENGTH)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Length")
pair = c(pair, "2 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data2$LENGTH, data3$LENGTH)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Length")
pair = c(pair, "2 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data2$LENGTH, data4$LENGTH)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Length")
pair = c(pair, "3 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data3$LENGTH, data4$LENGTH)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


#ASN... TESTS

trait = c(trait, "% ASN")
pair = c(pair, "1 vs 2")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$ASN../data1$LENGTH *100),]$ASN../data1$LENGTH*100, 
              data2[!flag_outliers(data2$ASN../data2$LENGTH *100),]$ASN../data2$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "% ASN")
pair = c(pair, "1 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$ASN../data1$LENGTH *100),]$ASN../data1$LENGTH*100, 
              data3[!flag_outliers(data3$ASN../data3$LENGTH *100),]$ASN../data3$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "% ASN")
pair = c(pair, "1 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$ASN../data1$LENGTH *100),]$ASN../data1$LENGTH*100, 
              data4[!flag_outliers(data4$ASN../data4$LENGTH *100),]$ASN../data4$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "% ASN")
pair = c(pair, "2 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$ASN../data2$LENGTH *100),]$ASN../data2$LENGTH*100, 
              data3[!flag_outliers(data3$ASN../data3$LENGTH *100),]$ASN../data3$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "% ASN")
pair = c(pair, "2 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$ASN../data2$LENGTH *100),]$ASN../data2$LENGTH*100, 
              data4[!flag_outliers(data4$ASN../data4$LENGTH *100),]$ASN../data4$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "% ASN")
pair = c(pair, "3 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data3[!flag_outliers(data3$ASN../data3$LENGTH *100),]$ASN../data3$LENGTH*100, 
              data4[!flag_outliers(data4$ASN../data4$LENGTH *100),]$ASN../data4$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])





#CYS.. TESTS

trait = c(trait, "% CYS")
pair = c(pair, "1 vs 2")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$CYS../data1$LENGTH *100),]$CYS../data1$LENGTH*100, 
              data2[!flag_outliers(data2$CYS../data2$LENGTH *100),]$CYS../data2$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

trait = c(trait, "% CYS")
pair = c(pair, "1 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$CYS../data1$LENGTH *100),]$CYS../data1$LENGTH*100, 
              data3[!flag_outliers(data3$CYS../data3$LENGTH *100),]$CYS../data3$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

trait = c(trait, "% CYS")
pair = c(pair, "1 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$CYS../data1$LENGTH *100),]$CYS../data1$LENGTH*100, 
              data4[!flag_outliers(data4$CYS../data4$LENGTH *100),]$CYS../data4$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

trait = c(trait, "% CYS")
pair = c(pair, "2 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$CYS../data2$LENGTH *100),]$CYS../data2$LENGTH*100, 
              data3[!flag_outliers(data3$CYS../data3$LENGTH *100),]$CYS../data3$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "% CYS")
pair = c(pair, "2 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$CYS../data2$LENGTH *100),]$CYS../data2$LENGTH*100, 
              data4[!flag_outliers(data4$CYS../data4$LENGTH *100),]$CYS../data3$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

trait = c(trait, "% CYS")
pair = c(pair, "3 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data3[!flag_outliers(data3$CYS../data3$LENGTH *100),]$CYS../data3$LENGTH*100, 
              data4[!flag_outliers(data4$CYS../data4$LENGTH *100),]$CYS../data4$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])



# SYM TESTS

trait = c(trait, "Symm_dist")
pair = c(pair, "1 vs 2")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$LENGTH %% 3),]$LENGTH %% 3, 
              data2[!flag_outliers(data2$LENGTH %% 3),]$LENGTH %% 3.)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

trait = c(trait, "Symm_dist")
pair = c(pair, "1 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$LENGTH %% 3),]$LENGTH %% 3, 
              data3[!flag_outliers(data3$LENGTH %% 3),]$LENGTH %% 3.)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Symm_dist")
pair = c(pair, "1 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$LENGTH %% 3),]$LENGTH %% 3, 
              data4[!flag_outliers(data4$LENGTH %% 3),]$LENGTH %% 3.)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Symm_dist")
pair = c(pair, "2 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$LENGTH %% 3),]$LENGTH %% 3, 
              data3[!flag_outliers(data3$LENGTH %% 3),]$LENGTH %% 3.)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

trait = c(trait, "Symm_dist")
pair = c(pair, "2 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$LENGTH %% 3),]$LENGTH %% 3, 
              data4[!flag_outliers(data4$LENGTH %% 3),]$LENGTH %% 3.)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Symm_dist")
pair = c(pair, "3 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data3[!flag_outliers(data3$LENGTH %% 3),]$LENGTH %% 3, 
              data4[!flag_outliers(data4$LENGTH %% 3),]$LENGTH %% 3.)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

# 
# 
# #MIN TESTS
# 
# trait = c(trait, "MIN_RSA")
# pair = c(pair, "1 vs 2")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MIN),]$MIN, data2[!flag_outliers(data2$MIN),]$MIN)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "MIN_RSA")
# pair = c(pair, "1 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MIN),]$MIN, data3[!flag_outliers(data3$MIN),]$MIN)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "MIN_RSA")
# pair = c(pair, "1 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MIN),]$MIN, data4[!flag_outliers(data4$MIN),]$MIN)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MIN_RSA")
# pair = c(pair, "2 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$MIN),]$MIN, data3[!flag_outliers(data3$MIN),]$MIN)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "MIN_RSA")
# pair = c(pair, "2 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$MIN),]$MIN, data4[!flag_outliers(data4$MIN),]$MIN)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MIN_RSA")
# pair = c(pair, "3 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data3[!flag_outliers(data3$MIN),]$MIN, data4[!flag_outliers(data4$MIN),]$MIN)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# 

#Q1 TESTS

trait = c(trait, "Q1_RSA")
pair = c(pair, "1 vs 2")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$Q1),]$Q1, data2[!flag_outliers(data2$Q1),]$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Q1_RSA")
pair = c(pair, "1 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$Q1),]$Q1, data3[!flag_outliers(data3$Q1),]$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Q1_RSA")
pair = c(pair, "1 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$Q1),]$Q1, data4[!flag_outliers(data4$Q1),]$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Q1_RSA")
pair = c(pair, "2 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$Q1),]$Q1, data3[!flag_outliers(data3$Q1),]$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

trait = c(trait, "Q1_RSA")
pair = c(pair, "2 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$Q1),]$Q1, data4[!flag_outliers(data4$Q1),]$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Q1_RSA")
pair = c(pair, "3 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data3[!flag_outliers(data3$Q1),]$Q1, data4[!flag_outliers(data4$Q1),]$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

# 
# 
# #Q2 TESTS
# 
# trait = c(trait, "Q2_RSA")
# pair = c(pair, "1 vs 2")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q2),]$Q2, data2[!flag_outliers(data2$Q2),]$Q2)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "Q2_RSA")
# pair = c(pair, "1 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q2),]$Q2, data3[!flag_outliers(data3$Q2),]$Q2)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q2_RSA")
# pair = c(pair, "1 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q2),]$Q2, data4[!flag_outliers(data4$Q2),]$Q2)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q2_RSA")
# pair = c(pair, "2 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$Q2),]$Q2, data3[!flag_outliers(data3$Q2),]$Q2)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q2_RSA")
# pair = c(pair, "2 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$Q2),]$Q2, data4[!flag_outliers(data4$Q2),]$Q2)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q2_RSA")
# pair = c(pair, "3 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data3[!flag_outliers(data3$Q2),]$Q2, data4[!flag_outliers(data4$Q2),]$Q2)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# #Q3 TESTS
# 
# trait = c(trait, "Q3_RSA")
# pair = c(pair, "1 vs 2")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q3),]$Q3, data2[!flag_outliers(data2$Q3),]$Q3)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q3_RSA")
# pair = c(pair, "1 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q3),]$Q3, data3[!flag_outliers(data3$Q3),]$Q3)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "Q3_RSA")
# pair = c(pair, "1 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q3),]$Q3, data4[!flag_outliers(data4$Q3),]$Q3)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "Q3_RSA")
# pair = c(pair, "2 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$Q3),]$Q3, data3[!flag_outliers(data3$Q3),]$Q3)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q3_RSA")
# pair = c(pair, "2 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$Q3),]$Q3, data4[!flag_outliers(data4$Q3),]$Q3)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "Q3_RSA")
# pair = c(pair, "3 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data3[!flag_outliers(data3$Q3),]$Q3, data4[!flag_outliers(data4$Q3),]$Q3)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# 
# 
# 
# #MAX TESTS
# 
# trait = c(trait, "MAX_RSA")
# pair = c(pair, "1 vs 2")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MAX),]$MAX, data2[!flag_outliers(data2$MAX),]$MAX)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MAX_RSA")
# pair = c(pair, "1 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MAX),]$MAX, data3[!flag_outliers(data3$MAX),]$MAX)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MAX_RSA")
# pair = c(pair, "1 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MAX),]$MAX, data4[!flag_outliers(data4$MAX),]$MAX)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MAX_RSA")
# pair = c(pair, "2 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$MAX),]$MAX, data3[!flag_outliers(data3$MAX),]$MAX)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "MAX_RSA")
# pair = c(pair, "2 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$MAX),]$MAX, data4[!flag_outliers(data4$MAX),]$MAX)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "MAX_RSA")
# pair = c(pair, "3 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data3[!flag_outliers(data3$MAX),]$MAX, data4[!flag_outliers(data4$MAX),]$MAX)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# 
# #MIN_pLDDT TESTS
# 
# trait = c(trait, "MIN_pLDDT")
# pair = c(pair, "1 vs 2")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MIN_pLLDT),]$MIN_pLLDT, 
#               data2[!flag_outliers(data2$MIN_pLLDT),]$MIN_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MIN_pLDDT")
# pair = c(pair, "1 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MIN_pLLDT),]$MIN_pLLDT, 
#               data3[!flag_outliers(data3$MIN_pLLDT),]$MIN_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MIN_pLDDT")
# pair = c(pair, "1 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MIN_pLLDT),]$MIN_pLLDT, 
#               data4[!flag_outliers(data4$MIN_pLLDT),]$MIN_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MIN_pLDDT")
# pair = c(pair, "2 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$MIN_pLLDT),]$MIN_pLLDT, 
#               data3[!flag_outliers(data3$MIN_pLLDT),]$MIN_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MIN_pLDDT")
# pair = c(pair, "2 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$MIN_pLLDT),]$MIN_pLLDT, 
#               data4[!flag_outliers(data4$MIN_pLLDT),]$MIN_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MIN_pLDDT")
# pair = c(pair, "3 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data3[!flag_outliers(data3$MIN_pLLDT),]$MIN_pLLDT, 
#               data4[!flag_outliers(data4$MIN_pLLDT),]$MIN_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# #Q1_pLDDT TESTS
# 
# trait = c(trait, "Q1_pLDDT")
# pair = c(pair, "1 vs 2")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q1_pLLDT),]$Q1_pLLDT, 
#               data2[!flag_outliers(data2$Q1_pLLDT),]$Q1_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q1_pLDDT")
# pair = c(pair, "1 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q1_pLLDT),]$Q1_pLLDT, 
#               data3[!flag_outliers(data3$Q1_pLLDT),]$Q1_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q1_pLDDT")
# pair = c(pair, "1 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q1_pLLDT),]$Q1_pLLDT, 
#               data4[!flag_outliers(data4$Q1_pLLDT),]$Q1_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q1_pLDDT")
# pair = c(pair, "2 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$Q1_pLLDT),]$Q1_pLLDT, 
#               data3[!flag_outliers(data3$Q1_pLLDT),]$Q1_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q1_pLDDT")
# pair = c(pair, "2 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$Q1_pLLDT),]$Q1_pLLDT, 
#               data4[!flag_outliers(data4$Q1_pLLDT),]$Q1_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q1_pLDDT")
# pair = c(pair, "3 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data3[!flag_outliers(data3$Q1_pLLDT),]$Q1_pLLDT, 
#               data4[!flag_outliers(data4$Q1_pLLDT),]$Q1_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 

# 
# #Q2_pLDDT TESTS
# 
# trait = c(trait, "Q2_pLDDT")
# pair = c(pair, "1 vs 2")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q2_pLLDT),]$Q2_pLLDT, 
#               data2[!flag_outliers(data2$Q2_pLLDT),]$Q2_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q2_pLDDT")
# pair = c(pair, "1 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q2_pLLDT),]$Q2_pLLDT, 
#               data3[!flag_outliers(data3$Q2_pLLDT),]$Q2_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "Q2_pLDDT")
# pair = c(pair, "1 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$Q2_pLLDT),]$Q2_pLLDT, 
#               data4[!flag_outliers(data4$Q2_pLLDT),]$Q2_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q2_pLDDT")
# pair = c(pair, "2 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$Q2_pLLDT),]$Q2_pLLDT, 
#               data3[!flag_outliers(data3$Q2_pLLDT),]$Q2_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q2_pLDDT")
# pair = c(pair, "2 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$Q2_pLLDT),]$Q2_pLLDT, 
#               data4[!flag_outliers(data4$Q2_pLLDT),]$Q2_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "Q2_pLDDT")
# pair = c(pair, "3 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data3[!flag_outliers(data3$Q2_pLLDT),]$Q2_pLLDT, 
#               data4[!flag_outliers(data4$Q2_pLLDT),]$Q2_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])


#Q3_pLDDT TESTS

trait = c(trait, "Q3_pLDDT")
pair = c(pair, "1 vs 2")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$Q3_pLLDT),]$Q3_pLLDT,
              data2[!flag_outliers(data2$Q3_pLLDT),]$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Q3_pLDDT")
pair = c(pair, "1 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$Q3_pLLDT),]$Q3_pLLDT,
              data3[!flag_outliers(data3$Q3_pLLDT),]$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Q3_pLDDT")
pair = c(pair, "1 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$Q3_pLLDT),]$Q3_pLLDT,
              data4[!flag_outliers(data4$Q3_pLLDT),]$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Q3_pLDDT")
pair = c(pair, "2 vs 3")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$Q3_pLLDT),]$Q3_pLLDT,
              data3[!flag_outliers(data3$Q3_pLLDT),]$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

trait = c(trait, "Q3_pLDDT")
pair = c(pair, "2 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data2[!flag_outliers(data2$Q3_pLLDT),]$Q3_pLLDT,
              data4[!flag_outliers(data4$Q3_pLLDT),]$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


trait = c(trait, "Q3_pLDDT")
pair = c(pair, "3 vs 4")
test_type = c(test_type, "ks")
out = ks_stat(data3[!flag_outliers(data3$Q3_pLLDT),]$Q3_pLLDT,
              data4[!flag_outliers(data4$Q3_pLLDT),]$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])



# #MAX_pLDDT TESTS
# 
# trait = c(trait, "MAX_pLDDT")
# pair = c(pair, "1 vs 2")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MAX_pLLDT),]$MAX_pLLDT, 
#               data2[!flag_outliers(data2$MAX_pLLDT),]$MAX_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MAX_pLDDT")
# pair = c(pair, "1 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MAX_pLLDT),]$MAX_pLLDT, 
#               data3[!flag_outliers(data3$MAX_pLLDT),]$MAX_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MAX_pLDDT")
# pair = c(pair, "1 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data1[!flag_outliers(data1$MAX_pLLDT),]$MAX_pLLDT, 
#               data4[!flag_outliers(data4$MAX_pLLDT),]$MAX_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MAX_pLDDT")
# pair = c(pair, "2 vs 3")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(dat2$MAX_pLLDT),]$MAX_pLLDT, 
#               data3[!flag_outliers(data3$MAX_pLLDT),]$MAX_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# trait = c(trait, "MAX_pLDDT")
# pair = c(pair, "2 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data2[!flag_outliers(data2$MAX_pLLDT),]$MAX_pLLDT, 
#               data4[!flag_outliers(data4$MAX_pLLDT),]$MAX_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# trait = c(trait, "MAX_pLDDT")
# pair = c(pair, "3 vs 4")
# test_type = c(test_type, "ks")
# out = ks_stat(data3[!flag_outliers(data3$MAX_pLLDT),]$MAX_pLLDT, 
#               data4[!flag_outliers(data4$MAX_pLLDT),]$MAX_pLLDT)
# pval = c(pval, out[1])
# conf_low = c(conf_low, out[3])
# conf_top = c(conf_top, out[4])
# stat = c(stat, out[2])
# 
# 
# #UNSTRUC TESTS
# unstruc = data.frame("sQTL" = c(nrow(data1_unstruc), nrow(data1) - nrow(data1_unstruc)),
#                      "non_sQTL" = c(nrow(data2_unstruc), nrow(data2) - nrow(data2_unstruc)),
#                      #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                      row.names = c("unstruc", "other"),
#                      stringsAsFactors = FALSE)
# 
# unstruc
# test = fisher.test(unstruc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Unstructured")
# pair = c(pair, "1 vs 2")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# unstruc = data.frame("sQTL" = c(nrow(data1_unstruc), nrow(data1) - nrow(data1_unstruc)),
#                      "non_sQTL" = c(nrow(data3_unstruc), nrow(data3) - nrow(data3_unstruc)),
#                      #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                      row.names = c("unstruc", "other"),
#                      stringsAsFactors = FALSE)
# 
# unstruc
# test = fisher.test(unstruc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Unstructured")
# pair = c(pair, "1 vs 3")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# unstruc = data.frame("sQTL" = c(nrow(data1_unstruc), nrow(data1) - nrow(data1_unstruc)),
#                      "non_sQTL" = c(nrow(data4_unstruc), nrow(data4) - nrow(data4_unstruc)),
#                      #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                      row.names = c("unstruc", "other"),
#                      stringsAsFactors = FALSE)
# 
# unstruc
# test = fisher.test(unstruc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Unstructured")
# pair = c(pair, "1 vs 4")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# unstruc = data.frame("sQTL" = c(nrow(data2_unstruc), nrow(data2) - nrow(data2_unstruc)),
#                      "non_sQTL" = c(nrow(data3_unstruc), nrow(data3) - nrow(data3_unstruc)),
#                      #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                      row.names = c("unstruc", "other"),
#                      stringsAsFactors = FALSE)
# 
# unstruc
# test = fisher.test(unstruc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Unstructured")
# pair = c(pair, "2 vs 3")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# unstruc = data.frame("sQTL" = c(nrow(data2_unstruc), nrow(data2) - nrow(data2_unstruc)),
#                      "non_sQTL" = c(nrow(data4_unstruc), nrow(data4) - nrow(data4_unstruc)),
#                      #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                      row.names = c("unstruc", "other"),
#                      stringsAsFactors = FALSE)
# 
# unstruc
# test = fisher.test(unstruc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Unstructured")
# pair = c(pair, "2 vs 4")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# unstruc = data.frame("sQTL" = c(nrow(data3_unstruc), nrow(data3) - nrow(data3_unstruc)),
#                      "non_sQTL" = c(nrow(data4_unstruc), nrow(data4) - nrow(data4_unstruc)),
#                      #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                      row.names = c("unstruc", "other"),
#                      stringsAsFactors = FALSE)
# 
# unstruc
# test = fisher.test(unstruc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Unstructured")
# pair = c(pair, "3 vs 4")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# 
# 
# 
# # STRUC TESTS
# struc = data.frame("sQTL" = c(nrow(data1_struc), nrow(data1)-nrow(data1_struc)),
#                    "non_sQTL" = c(nrow(data2_struc), nrow(data2)-nrow(data2_struc)),
#                    #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                    row.names = c("struc", "other"),
#                    stringsAsFactors = FALSE)
# 
# struc
# test = fisher.test(struc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Structured")
# pair = c(pair, "1 vs 2")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# struc = data.frame("sQTL" = c(nrow(data1_struc), nrow(data1)-nrow(data1_struc)),
#                    "non_sQTL" = c(nrow(data3_struc), nrow(data3)-nrow(data3_struc)),
#                    #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                    row.names = c("struc", "other"),
#                    stringsAsFactors = FALSE)
# 
# struc
# test = fisher.test(struc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Structured")
# pair = c(pair, "1 vs 3")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# struc = data.frame("sQTL" = c(nrow(data1_struc), nrow(data1)-nrow(data1_struc)),
#                    "non_sQTL" = c(nrow(data4_struc), nrow(data4)-nrow(data4_struc)),
#                    #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                    row.names = c("struc", "other"),
#                    stringsAsFactors = FALSE)
# 
# struc
# test = fisher.test(struc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Structured")
# pair = c(pair, "1 vs 4")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# struc = data.frame("sQTL" = c(nrow(data2_struc), nrow(data2)-nrow(data2_struc)),
#                    "non_sQTL" = c(nrow(data3_struc), nrow(data3)-nrow(data3_struc)),
#                    #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                    row.names = c("struc", "other"),
#                    stringsAsFactors = FALSE)
# 
# struc
# test = fisher.test(struc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Structured")
# pair = c(pair, "2 vs 3")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# struc = data.frame("sQTL" = c(nrow(data2_struc), nrow(data2)-nrow(data2_struc)),
#                    "non_sQTL" = c(nrow(data4_struc), nrow(data4)-nrow(data4_struc)),
#                    #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                    row.names = c("struc", "other"),
#                    stringsAsFactors = FALSE)
# 
# struc
# test = fisher.test(struc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Structured")
# pair = c(pair, "2 vs 4")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# struc = data.frame("sQTL" = c(nrow(data3_struc), nrow(data3)-nrow(data3_struc)),
#                    "non_sQTL" = c(nrow(data4_struc), nrow(data4)-nrow(data4_struc)),
#                    #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                    row.names = c("struc", "other"),
#                    stringsAsFactors = FALSE)
# 
# struc
# test = fisher.test(struc)
# pval = c(pval, test$p.value)
# trait = c(trait, "Structured")
# pair = c(pair, "3 vs 4")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# 
# 
# # CONSERVED TESTS
# cons = data.frame("sQTL" = c(nrow(data1_cons), nrow(data1) - nrow(data1_cons)),
#                   "non_sQTL" = c(nrow(data2_cons), nrow(data2) - nrow(data2_cons)),
#                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                   row.names = c("cons", "other"),
#                   stringsAsFactors = FALSE)
# 
# cons
# test = fisher.test(cons)
# pval = c(pval, test$p.value)
# trait = c(trait, "Conserved")
# pair = c(pair, "1 vs 2")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# cons = data.frame("sQTL" = c(nrow(data1_cons), nrow(data1) - nrow(data1_cons)),
#                   "non_sQTL" = c(nrow(data3_cons), nrow(data3) - nrow(data3_cons)),
#                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                   row.names = c("cons", "other"),
#                   stringsAsFactors = FALSE)
# 
# cons
# test = fisher.test(cons)
# pval = c(pval, test$p.value)
# trait = c(trait, "Conserved")
# pair = c(pair, "1 vs 3")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# cons = data.frame("sQTL" = c(nrow(data1_cons), nrow(data1) - nrow(data1_cons)),
#                   "non_sQTL" = c(nrow(data4_cons), nrow(data4) - nrow(data4_cons)),
#                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                   row.names = c("cons", "other"),
#                   stringsAsFactors = FALSE)
# 
# cons
# test = fisher.test(cons)
# pval = c(pval, test$p.value)
# trait = c(trait, "Conserved")
# pair = c(pair, "1 vs 4")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# cons = data.frame("sQTL" = c(nrow(data2_cons), nrow(data2) - nrow(data2_cons)),
#                   "non_sQTL" = c(nrow(data3_cons), nrow(data3) - nrow(data3_cons)),
#                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                   row.names = c("cons", "other"),
#                   stringsAsFactors = FALSE)
# 
# cons
# test = fisher.test(cons)
# pval = c(pval, test$p.value)
# trait = c(trait, "Conserved")
# pair = c(pair, "2 vs 3")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# cons = data.frame("sQTL" = c(nrow(data2_cons), nrow(data2) - nrow(data2_cons)),
#                   "non_sQTL" = c(nrow(data4_cons), nrow(data4) - nrow(data4_cons)),
#                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                   row.names = c("cons", "other"),
#                   stringsAsFactors = FALSE)
# 
# cons
# test = fisher.test(cons)
# pval = c(pval, test$p.value)
# trait = c(trait, "Conserved")
# pair = c(pair, "2 vs 4")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")
# 
# cons = data.frame("sQTL" = c(nrow(data3_cons), nrow(data3) - nrow(data3_cons)),
#                   "non_sQTL" = c(nrow(data4_cons), nrow(data4) - nrow(data4_cons)),
#                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
#                   row.names = c("cons", "other"),
#                   stringsAsFactors = FALSE)
# 
# cons
# test = fisher.test(cons)
# pval = c(pval, test$p.value)
# trait = c(trait, "Conserved")
# pair = c(pair, "3 vs 4")
# conf_low = c(conf_low, test$conf.int[1])
# conf_top = c(conf_top, test$conf.int[2])
# stat = c(stat, test$estimate)
# test_type = c(test_type, "fisher")



# SYM TESTS
symet = data.frame("sQTL" = c(length(data1[data1$LENGTH %% 3 == 0,]$LENGTH), length(data1[data1$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_sQTL" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
pval = c(pval, test$p.value)
trait = c(trait, "Symmetric")
pair = c(pair, "1 vs 2")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet = data.frame("sQTL" = c(length(data1[data1$LENGTH %% 3 == 0,]$LENGTH), length(data1[data1$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
pval = c(pval, test$p.value)
trait = c(trait, "Symmetric")
pair = c(pair, "1 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet = data.frame("sQTL" = c(length(data1[data1$LENGTH %% 3 == 0,]$LENGTH), length(data1[data1$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_sQTL" = c(length(data4[data4$LENGTH %% 3 == 0,]$LENGTH), length(data4[data4$LENGTH %% 3 != 0,]$LENGTH)),
                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
pval = c(pval, test$p.value)
trait = c(trait, "Symmetric")
pair = c(pair, "1 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet = data.frame("sQTL" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
pval = c(pval, test$p.value)
trait = c(trait, "Symmetric")
pair = c(pair, "2 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet = data.frame("sQTL" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_sQTL" = c(length(data4[data4$LENGTH %% 3 == 0,]$LENGTH), length(data4[data4$LENGTH %% 3 != 0,]$LENGTH)),
                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
pval = c(pval, test$p.value)
trait = c(trait, "Symmetric")
pair = c(pair, "2 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet = data.frame("sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_sQTL" = c(length(data4[data4$LENGTH %% 3 == 0,]$LENGTH), length(data4[data4$LENGTH %% 3 != 0,]$LENGTH)),
                   #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0,]$LENGTH), length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
pval = c(pval, test$p.value)
trait = c(trait, "Symmetric")
pair = c(pair, "3 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")



pa = c("1", "2", "3", "4", "1", "2", "3", "4")
tra = c("Symmetric", "Symmetric", "Symmetric", "Symmetric","Non symmetric", "Non symmetric", "Non symmetric", "Non symmetric")
count = c(length(data1[data1$LENGTH %% 3 == 0,]$LENGTH)/length(data1$LENGTH),
          length(data2[data2$LENGTH %% 3 == 0,]$LENGTH)/length(data2$LENGTH),
          length(data3[data3$LENGTH %% 3 == 0,]$LENGTH)/length(data3$LENGTH),
          length(data4[data4$LENGTH %% 3 == 0,]$LENGTH)/length(data4$LENGTH), 
          length(data1[data1$LENGTH %% 3 != 0,]$LENGTH)/length(data1$LENGTH),  
          length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)/length(data2$LENGTH),
          length(data3[data3$LENGTH %% 3 != 0,]$LENGTH)/length(data3$LENGTH),  
          length(data4[data4$LENGTH %% 3 != 0,]$LENGTH)/length(data4$LENGTH))
to_draw = cbind(tra, count)
to_draw = cbind(to_draw, pa)
to_draw
to_draw = data.frame(to_draw)
colnames(to_draw) = c("trait", "proportion", "Group")
to_draw$proportion = as.numeric(to_draw$proportion)
to_draw 
ggplot(to_draw, aes(fill = as.factor(Group), x=trait, y=proportion )) + gtex_v8_figure_theme() + geom_col(position = position_dodge())  +  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("proportion") + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
ggsave("symmetric_bars_among_groups.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)


# SYM NOT FOUND TESTS
symet_not_found = data.frame("sQTL" = c(length(data1[data1$LENGTH %% 3 == 0 & data1$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data1[data1$LENGTH %% 3 != 0 & data1$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             "non_sQTL" = c(length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data2[data2$LENGTH %% 3 != 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             row.names = c("sym", "non_sym"),
                             stringsAsFactors = FALSE)
symet_not_found
test = fisher.test(symet_not_found)
pval = c(pval, test$p.value)
trait = c(trait, "Symm_not_found")
pair = c(pair, "1 vs 2")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet_not_found = data.frame("sQTL" = c(length(data1[data1$LENGTH %% 3 == 0 & data1$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data1[data1$LENGTH %% 3 != 0 & data1$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             "non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             row.names = c("sym", "non_sym"),
                             stringsAsFactors = FALSE)
symet_not_found
test = fisher.test(symet_not_found)
pval = c(pval, test$p.value)
trait = c(trait, "Symm_not_found")
pair = c(pair, "1 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet_not_found = data.frame("sQTL" = c(length(data1[data1$LENGTH %% 3 == 0 & data1$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data1[data1$LENGTH %% 3 != 0 & data1$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             "non_sQTL" = c(length(data4[data4$LENGTH %% 3 == 0 & data4$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data4[data4$LENGTH %% 3 != 0 & data4$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             row.names = c("sym", "non_sym"),
                             stringsAsFactors = FALSE)
symet_not_found
test = fisher.test(symet_not_found)
pval = c(pval, test$p.value)
trait = c(trait, "Symm_not_found")
pair = c(pair, "1 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet_not_found = data.frame("sQTL" = c(length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data2[data2$LENGTH %% 3 != 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             "non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             row.names = c("sym", "non_sym"),
                             stringsAsFactors = FALSE)
symet_not_found
test = fisher.test(symet_not_found)
pval = c(pval, test$p.value)
trait = c(trait, "Symm_not_found")
pair = c(pair, "2 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet_not_found = data.frame("sQTL" = c(length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data2[data2$LENGTH %% 3 != 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             "non_sQTL" = c(length(data4[data4$LENGTH %% 3 == 0 & data4$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data4[data4$LENGTH %% 3 != 0 & data4$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             row.names = c("sym", "non_sym"),
                             stringsAsFactors = FALSE)
symet_not_found
test = fisher.test(symet_not_found)
pval = c(pval, test$p.value)
trait = c(trait, "Symm_not_found")
pair = c(pair, "2 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

symet_not_found = data.frame("sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             "non_sQTL" = c(length(data4[data4$LENGTH %% 3 == 0 & data4$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data4[data4$LENGTH %% 3 != 0 & data4$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                             row.names = c("sym", "non_sym"),
                             stringsAsFactors = FALSE)
symet_not_found
test = fisher.test(symet_not_found)
pval = c(pval, test$p.value)
trait = c(trait, "Symm_not_found")
pair = c(pair, "3 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

pa = c("1", "2", "3", "4", "1", "2", "3", "4")
tra = c("Symmetric", "Symmetric", "Symmetric", "Symmetric", "Non symmetric", "Non symmetric", "Non symmetric", "Non symmetric")
count = c(length(data1[data1$LENGTH %% 3 == 0 & data1$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data1[data1$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),
          length(data2[data2$LENGTH %% 3 == 0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data2[data2$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),
          length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data3[data3$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),
          length(data4[data4$LENGTH %% 3 == 0 & data4$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data4[data4$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH), 
          length(data1[data1$LENGTH %% 3 != 0 & data1$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data1[data1$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),  
          length(data2[data2$LENGTH %% 3 !=  0 & data2$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data2[data2$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),
          length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data3[data3$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH),  
          length(data4[data4$LENGTH %% 3 !=  0 & data4$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)/length(data4[data4$ALIGNED.SEQ == "EXON NOT FOUND", ]$LENGTH))

to_draw = data.frame()
to_draw = cbind(tra, count)
to_draw = cbind(to_draw, pa)
to_draw
to_draw = data.frame(to_draw)
colnames(to_draw) = c("trait", "proportion", "Group")
to_draw$proportion = as.numeric(to_draw$proportion)
to_draw 
ggplot(to_draw, aes(fill = as.factor(Group), x=trait, y=proportion )) + gtex_v8_figure_theme() + geom_col(position = position_dodge())  +  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("proportion") + labs(fill='Group') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
ggsave("symmetric_not_found_bars_among_groups.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)



# HELIX TESTS
hel = data.frame("sQTL" = c(length(data1[data1$HELIX == "True",]$HELIX), length(data1[data1$HELIX != "True",]$HELIX)),
                 "non_sQTL" = c(length(data2[data2$HELIX == "True",]$HELIX), length(data2[data2$HELIX != "True",]$HELIX)),
                 #"non_sQTL" = c(length(data3[data3$HELIX ="True",]$HELIX), length(data3[data3$HELIX !"True",]$HELIX)),
                 row.names = c("helix", "no_helix"),
                 stringsAsFactors = FALSE)

hel
test = fisher.test(hel)
pval = c(pval, test$p.value)
trait = c(trait, "Helix")
pair = c(pair, "1 vs 2")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

hel = data.frame("sQTL" = c(length(data1[data1$HELIX == "True",]$HELIX), length(data1[data1$HELIX != "True",]$HELIX)),
                 "non_sQTL" = c(length(data3[data3$HELIX == "True",]$HELIX), length(data3[data3$HELIX != "True",]$HELIX)),
                 #"non_sQTL" = c(length(data3[data3$HELIX ="True",]$HELIX), length(data3[data3$HELIX !"True",]$HELIX)),
                 row.names = c("helix", "no_helix"),
                 stringsAsFactors = FALSE)

hel
test = fisher.test(hel)
pval = c(pval, test$p.value)
trait = c(trait, "Helix")
pair = c(pair, "1 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

hel = data.frame("sQTL" = c(length(data1[data1$HELIX == "True",]$HELIX), length(data1[data1$HELIX != "True",]$HELIX)),
                 "non_sQTL" = c(length(data4[data4$HELIX == "True",]$HELIX), length(data4[data4$HELIX != "True",]$HELIX)),
                 #"non_sQTL" = c(length(data3[data3$HELIX ="True",]$HELIX), length(data3[data3$HELIX !"True",]$HELIX)),
                 row.names = c("helix", "no_helix"),
                 stringsAsFactors = FALSE)

hel
test = fisher.test(hel)
pval = c(pval, test$p.value)
trait = c(trait, "Helix")
pair = c(pair, "1 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

hel = data.frame("sQTL" = c(length(data2[data2$HELIX == "True",]$HELIX), length(data2[data2$HELIX != "True",]$HELIX)),
                 "non_sQTL" = c(length(data3[data3$HELIX == "True",]$HELIX), length(data3[data3$HELIX != "True",]$HELIX)),
                 #"non_sQTL" = c(length(data3[data3$HELIX ="True",]$HELIX), length(data3[data3$HELIX !"True",]$HELIX)),
                 row.names = c("helix", "no_helix"),
                 stringsAsFactors = FALSE)

hel
test = fisher.test(hel)
pval = c(pval, test$p.value)
trait = c(trait, "Helix")
pair = c(pair, "2 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

hel = data.frame("sQTL" = c(length(data2[data2$HELIX == "True",]$HELIX), length(data2[data2$HELIX != "True",]$HELIX)),
                 "non_sQTL" = c(length(data4[data4$HELIX == "True",]$HELIX), length(data4[data4$HELIX != "True",]$HELIX)),
                 #"non_sQTL" = c(length(data3[data3$HELIX ="True",]$HELIX), length(data3[data3$HELIX !"True",]$HELIX)),
                 row.names = c("helix", "no_helix"),
                 stringsAsFactors = FALSE)

hel
test = fisher.test(hel)
pval = c(pval, test$p.value)
trait = c(trait, "Helix")
pair = c(pair, "2 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

hel = data.frame("sQTL" = c(length(data3[data3$HELIX == "True",]$HELIX), length(data3[data3$HELIX != "True",]$HELIX)),
                 "non_sQTL" = c(length(data4[data4$HELIX == "True",]$HELIX), length(data4[data4$HELIX != "True",]$HELIX)),
                 #"non_sQTL" = c(length(data3[data3$HELIX ="True",]$HELIX), length(data3[data3$HELIX !"True",]$HELIX)),
                 row.names = c("helix", "no_helix"),
                 stringsAsFactors = FALSE)

hel
test = fisher.test(hel)
pval = c(pval, test$p.value)
trait = c(trait, "Helix")
pair = c(pair, "3 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")





# SHEETS TESTS
she = data.frame("sQTL" = c(length(data1[data1$SHEET == "True",]$SHEET), length(data1[data1$SHEET != "True",]$SHEET)),
                 "non_sQTL" = c(length(data2[data2$SHEET == "True",]$SHEET), length(data2[data2$SHEET != "True",]$SHEET)),
                 #"non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 row.names = c("sheet", "no_sheet"),
                 stringsAsFactors = FALSE)

she
test = fisher.test(she)
pval = c(pval, test$p.value)
trait = c(trait, "Sheet")
pair = c(pair, "1 vs 2")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

she = data.frame("sQTL" = c(length(data1[data1$SHEET == "True",]$SHEET), length(data1[data1$SHEET != "True",]$SHEET)),
                 "non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 #"non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 row.names = c("sheet", "no_sheet"),
                 stringsAsFactors = FALSE)

she
test = fisher.test(she)
pval = c(pval, test$p.value)
trait = c(trait, "Sheet")
pair = c(pair, "1 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

she = data.frame("sQTL" = c(length(data1[data1$SHEET == "True",]$SHEET), length(data1[data1$SHEET != "True",]$SHEET)),
                 "non_sQTL" = c(length(data4[data4$SHEET == "True",]$SHEET), length(data4[data4$SHEET != "True",]$SHEET)),
                 #"non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 row.names = c("sheet", "no_sheet"),
                 stringsAsFactors = FALSE)

she
test = fisher.test(she)
pval = c(pval, test$p.value)
trait = c(trait, "Sheet")
pair = c(pair, "1 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

she = data.frame("sQTL" = c(length(data2[data2$SHEET == "True",]$SHEET), length(data2[data2$SHEET != "True",]$SHEET)),
                 "non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 #"non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 row.names = c("sheet", "no_sheet"),
                 stringsAsFactors = FALSE)

she
test = fisher.test(she)
pval = c(pval, test$p.value)
trait = c(trait, "Sheet")
pair = c(pair, "2 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

she = data.frame("sQTL" = c(length(data2[data2$SHEET == "True",]$SHEET), length(data2[data2$SHEET != "True",]$SHEET)),
                 "non_sQTL" = c(length(data4[data4$SHEET == "True",]$SHEET), length(data4[data4$SHEET != "True",]$SHEET)),
                 #"non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 row.names = c("sheet", "no_sheet"),
                 stringsAsFactors = FALSE)

she
test = fisher.test(she)
pval = c(pval, test$p.value)
trait = c(trait, "Sheet")
pair = c(pair, "2 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

she = data.frame("sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 "non_sQTL" = c(length(data4[data4$SHEET == "True",]$SHEET), length(data4[data4$SHEET != "True",]$SHEET)),
                 #"non_sQTL" = c(length(data3[data3$SHEET == "True",]$SHEET), length(data3[data3$SHEET != "True",]$SHEET)),
                 row.names = c("sheet", "no_sheet"),
                 stringsAsFactors = FALSE)

she
test = fisher.test(she)
pval = c(pval, test$p.value)
trait = c(trait, "Sheet")
pair = c(pair, "3 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")



# TURNS TESTS
turn = data.frame("sQTL" = c(length(data1[data1$TURN == "True",]$TURN), length(data1[data1$TURN != "True",]$TURN)),
                  "non_sQTL" = c(length(data2[data2$TURN == "True",]$TURN), length(data2[data2$TURN != "True",]$TURN)),
                  #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("turn", "no_turn"),
                  stringsAsFactors = FALSE)

turn
test = fisher.test(turn)
pval = c(pval, test$p.value)
trait = c(trait, "Turn")
pair = c(pair, "1 vs 2")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

turn = data.frame("sQTL" = c(length(data1[data1$TURN == "True",]$TURN), length(data1[data1$TURN != "True",]$TURN)),
                  "non_sQTL" = c(length(data3[data3$TURN == "True",]$TURN), length(data3[data3$TURN != "True",]$TURN)),
                  #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("turn", "no_turn"),
                  stringsAsFactors = FALSE)

turn
test = fisher.test(turn)
pval = c(pval, test$p.value)
trait = c(trait, "Turn")
pair = c(pair, "1 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

turn = data.frame("sQTL" = c(length(data1[data1$TURN == "True",]$TURN), length(data1[data1$TURN != "True",]$TURN)),
                  "non_sQTL" = c(length(data4[data4$TURN == "True",]$TURN), length(data4[data4$TURN != "True",]$TURN)),
                  #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("turn", "no_turn"),
                  stringsAsFactors = FALSE)

turn
test = fisher.test(turn)
pval = c(pval, test$p.value)
trait = c(trait, "Turn")
pair = c(pair, "1 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

turn = data.frame("sQTL" = c(length(data2[data2$TURN == "True",]$TURN), length(data2[data2$TURN != "True",]$TURN)),
                  "non_sQTL" = c(length(data3[data3$TURN == "True",]$TURN), length(data3[data3$TURN != "True",]$TURN)),
                  #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("turn", "no_turn"),
                  stringsAsFactors = FALSE)

turn
test = fisher.test(turn)
pval = c(pval, test$p.value)
trait = c(trait, "Turn")
pair = c(pair, "2 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

turn = data.frame("sQTL" = c(length(data2[data2$TURN == "True",]$TURN), length(data2[data2$TURN != "True",]$TURN)),
                  "non_sQTL" = c(length(data4[data4$TURN == "True",]$TURN), length(data4[data4$TURN != "True",]$TURN)),
                  #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("turn", "no_turn"),
                  stringsAsFactors = FALSE)

turn
test = fisher.test(turn)
pval = c(pval, test$p.value)
trait = c(trait, "Turn")
pair = c(pair, "2 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

turn = data.frame("sQTL" = c(length(data3[data3$TURN == "True",]$TURN), length(data3[data3$TURN != "True",]$TURN)),
                  "non_sQTL" = c(length(data4[data4$TURN == "True",]$TURN), length(data4[data4$TURN != "True",]$TURN)),
                  #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("turn", "no_turn"),
                  stringsAsFactors = FALSE)

turn
test = fisher.test(turn)
pval = c(pval, test$p.value)
trait = c(trait, "Turn")
pair = c(pair, "3 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

#TRANSMEMBRANE TEST

trans= data.frame("sQTL" = c(length(data1[data1$TRANSMEMBRANE == "True",]$LENGTH), length(data1[data1$TRANSMEMBRANE != "True",]$LENGTH)),
                  "non_sQTL" = c(length(data2[data2$TRANSMEMBRANE == "True",]$LENGTH), length(data2[data2$TRANSMEMBRANE != "True",]$LENGTH)),
                  #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                  row.names = c("trans", "non_trans"),
                  stringsAsFactors = FALSE)
trans
test = fisher.test(trans)
pval = c(pval, test$p.value)
trait = c(trait, "Transmembrane")
pair = c(pair, "1 vs 2")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

trans= data.frame("sQTL" = c(length(data1[data1$TRANSMEMBRANE == "True",]$LENGTH), length(data1[data1$TRANSMEMBRANE != "True",]$LENGTH)),
                  "non_sQTL" = c(length(data3[data3$TRANSMEMBRANE == "True",]$LENGTH), length(data3[data3$TRANSMEMBRANE != "True",]$LENGTH)),
                  #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                  row.names = c("trans", "non_trans"),
                  stringsAsFactors = FALSE)
trans
test = fisher.test(trans)
pval = c(pval, test$p.value)
trait = c(trait, "Transmembrane")
pair = c(pair, "1 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

trans= data.frame("sQTL" = c(length(data1[data1$TRANSMEMBRANE == "True",]$LENGTH), length(data1[data1$TRANSMEMBRANE != "True",]$LENGTH)),
                  "non_sQTL" = c(length(data4[data4$TRANSMEMBRANE == "True",]$LENGTH), length(data4[data4$TRANSMEMBRANE != "True",]$LENGTH)),
                  #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                  row.names = c("trans", "non_trans"),
                  stringsAsFactors = FALSE)
trans
test = fisher.test(trans)
pval = c(pval, test$p.value)
trait = c(trait, "Transmembrane")
pair = c(pair, "1 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

trans= data.frame("sQTL" = c(length(data2[data2$TRANSMEMBRANE == "True",]$LENGTH), length(data2[data2$TRANSMEMBRANE != "True",]$LENGTH)),
                  "non_sQTL" = c(length(data3[data3$TRANSMEMBRANE == "True",]$LENGTH), length(data3[data3$TRANSMEMBRANE != "True",]$LENGTH)),
                  #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                  row.names = c("trans", "non_trans"),
                  stringsAsFactors = FALSE)
trans
test = fisher.test(trans)
pval = c(pval, test$p.value)
trait = c(trait, "Transmembrane")
pair = c(pair, "2 vs 3")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

trans= data.frame("sQTL" = c(length(data2[data2$TRANSMEMBRANE == "True",]$LENGTH), length(data2[data2$TRANSMEMBRANE != "True",]$LENGTH)),
                  "non_sQTL" = c(length(data4[data4$TRANSMEMBRANE == "True",]$LENGTH), length(data4[data4$TRANSMEMBRANE != "True",]$LENGTH)),
                  #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                  row.names = c("trans", "non_trans"),
                  stringsAsFactors = FALSE)
trans
test = fisher.test(trans)
pval = c(pval, test$p.value)
trait = c(trait, "Transmembrane")
pair = c(pair, "2 vs 4")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")

trans= data.frame("sQTL" = c(length(data3[data3$TRANSMEMBRANE == "True",]$LENGTH), length(data3[data3$TRANSMEMBRANE != "True",]$LENGTH)),
                  "non_sQTL" = c(length(data4[data4$TRANSMEMBRANE == "True",]$LENGTH), length(data4[data4$TRANSMEMBRANE != "True",]$LENGTH)),
                  #"non_sQTL" = c(length(data3[data3$LENGTH %% 3 == 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH), length(data3[data3$LENGTH %% 3 != 0 & data3$ALIGNED.SEQ == "EXON NOT FOUND",]$LENGTH)),
                  row.names = c("trans", "non_trans"),
                  stringsAsFactors = FALSE)
trans
test = fisher.test(trans)
pval = c(pval, test$p.value)
trait = c(trait, "Transmembrane")
pair = c(pair, "3 vs 4")
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


#pval = p.adjust(pval, method="hochberg")
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
to_draw$pair
groups = unlist(strsplit(to_draw$pair, ' vs '))
groups
to_draw$group1 = as.factor(groups[seq(from=1, to=length(groups), by=2)])
to_draw$group2 = as.factor(groups[seq(from=2, to=length(groups), by=2)])
View(to_draw)
to_draw$conf_low = as.numeric(to_draw$conf_low)
to_draw$conf_top = as.numeric(to_draw$conf_top)
to_draw$alpha = ifelse((to_draw$test_type == "ks"  & ((to_draw$conf_low <= 0 & to_draw$conf_top >= 0) | 
                                                              (to_draw$conf_low >= 0 & to_draw$conf_top <= 0))) | 
                               (to_draw$test_type == "fisher"  & ((log(to_draw$conf_low) <= 0 & log(to_draw$conf_top) >= 0) | 
                                                                          (log(to_draw$conf_low) >= 0 & log(to_draw$conf_top) <= 0))), yes = 0.8, 1.0)



line = replicate(length(pval), 1.5)
ggplot(to_draw, aes(fill = as.factor(pair), x=trait, y=pval )) + 
        geom_col(position = position_dodge())  + geom_hline(yintercept = line) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), 
              legend.text = element_text(size=7), axis.title.y = element_text(size = 8), 
              axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("-log10(p_value)") + 
        labs(fill='Pairs') + 
        guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5))) + gtex_v8_figure_theme()

line1 = replicate(length(pval), 0.0)

ggplot(to_draw[to_draw$pair == "1 vs 2" & to_draw$test_type == "fisher" , ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) +  ylab("Enrichment in group 2                                                                Enrichment in group 1\n log(odd_ratio)")+
        coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("1 vs 2") + gtex_v8_figure_theme() + guides(alpha = FALSE)
ggsave("statistical_summary_fisher_groups_1_2.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "1 vs 3" & to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) + ylab("Enrichment in group 3                                                                Enrichment in group 1\n log(odd_ratio)")+ coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("1 vs 3") + gtex_v8_figure_theme()
ggsave("statistical_summary_fisher_groups_1_3.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "1 vs 4" & to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) +  ylab("Enrichment in group 4                                                                Enrichment in group 1\n log(odd_ratio)")+coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("1 vs 4") + gtex_v8_figure_theme()
ggsave("statistical_summary_fisher_groups_1_4.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "2 vs 3" & to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) +  ylab("Enrichment in group 3                                                                Enrichment in group 2\n log(odd_ratio)")+coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("2 vs 3") + gtex_v8_figure_theme()
ggsave("statistical_summary_fisher_groups_2_3.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "2 vs 4" & to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) + ylab("Enrichment in group 4                                                                Enrichment in group 2\n log(odd_ratio)")+ coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("2 vs 4") + gtex_v8_figure_theme()
ggsave("statistical_summary_fisher_groups_2_4.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "3 vs 4" & to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) +  ylab("Enrichment in group 4                                                                Enrichment in group 3\n log(odd_ratio)")+ coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("3 vs 4") + gtex_v8_figure_theme()
ggsave("statistical_summary_fisher_groups_3_4.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)



ggplot(to_draw[to_draw$pair == "1 vs 2" & to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + ylab("Enrichment in group 2                                                                 Enrichment in group 1\n m1 - m2") +  coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("1 vs 2") + gtex_v8_figure_theme()
ggsave("statistical_summary_ks_groups_1_2.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "1 vs 3" & to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha)) +
        geom_point(size = 4)  + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + ylab("Enrichment in group 3                                                                 Enrichment in group 1\n m1 - m2") + coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("1 vs 3") + gtex_v8_figure_theme()
ggsave("statistical_summary_ks_groups_1_3.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "1 vs 4" & to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + ylab("Enrichment in group 4                                                                 Enrichment in group 1\n m1 - m2") + coord_flip()  +geom_hline(yintercept = line1, color="red") + ggtitle("1 vs 4") + gtex_v8_figure_theme()
ggsave("statistical_summary_ks_groups_1_4.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "2 vs 3" & to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + ylab("Enrichment in group 3                                                                 Enrichment in group 2\n m1 - m2") + coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("2 vs 3") + gtex_v8_figure_theme()
ggsave("statistical_summary_ks_groups_2_3.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "2 vs 4" & to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + ylab("Enrichment in group 4                                                                 Enrichment in group 2\n m1 - m2") + coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("2 vs 4") + gtex_v8_figure_theme()
ggsave("statistical_summary_ks_groups_2_4.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
ggplot(to_draw[to_draw$pair == "3 vs 4" & to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha)) +
        geom_point(size = 4) + guides(alpha = FALSE) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + ylab("Enrichment in group 4                                                                 Enrichment in group 3\n m1 - m2") + coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("3 vs 4") + gtex_v8_figure_theme()
ggsave("statistical_summary_ks_groups_3_4.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)

#ggsave("statistical_summary_with_correction_sQTL_among_groups.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)


ggplot(to_draw[to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), color=pair, alpha = alpha)) +
        facet_grid(group1 ~ group2) + 
        geom_point(size = 4, position = position_dodge(width=0.9)) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low), 
                          color=pair), position=position_dodge(width=0.9)) + 
        ylab("Enrichment in trait 2                                                                 Enrichment in trait 1\n m1 - m2")+ 
        coord_flip()  + geom_hline(yintercept = line1) + guides(alpha = FALSE) + gtex_v8_figure_theme()
ggsave("statistical_summary_ks_groups_all.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)




ggplot(to_draw[to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), color=pair, alpha = alpha)) +
        geom_point(size = 4, position=position_dodge2(width = 0.9)) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)), color=pair), 
                      position=position_dodge2(width = 0.9)) + 
        ylab("Enrichment in trait 2                                                                Enrichment in trait 1\n log(odd_ratio)")+
        coord_flip()  + geom_hline(yintercept = line1) + 
        facet_grid(group1 ~ group2) + guides(alpha = FALSE) + gtex_v8_figure_theme()
ggsave("statistical_summary_fisher_groups_all.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)

box1 = c()
box2 = c()
box3 = c()
box4 = c()
box1 = cbind(data1$anc_psi, "AA")
box1 = rbind(box1, cbind(data1$mean_01_psi, "AD"))
box1 = rbind(box1, cbind(data1$der_psi, "DD"))
box1 = data.frame(box1)

box2 = cbind(data2$anc_psi, "AA")
box2 = rbind(box2, cbind(data2$mean_01_psi, "AD"))
box2 = rbind(box2, cbind(data2$der_psi, "DD"))
box2 = data.frame(box2)

box3 = cbind(data3$anc_psi, "AA")
box3 = rbind(box3, cbind(data3$mean_01_psi, "AD"))
box3 = rbind(box3, cbind(data3$der_psi, "DD"))
box3 = data.frame(box3)

box4 = cbind(data4$anc_psi, "AA")
box4 = rbind(box4, cbind(data4$mean_01_psi, "AD"))
box4 = rbind(box4, cbind(data4$der_psi, "DD"))
box4 = data.frame(box4)
dev.off()
plot1 = ggplot(box1, aes(x=as.numeric(X1), group=as.factor(X2))) +
        geom_boxplot(fill="2")+ xlab("PSI") +  ylab("Genotype") + ggtitle("Derived Higher Inclusion <50% mean PSI") + gtex_v8_figure_theme() + coord_flip() + scale_x_discrete(breaks=c("-0.2","0.0","0.2"),
                                                                                                                                                                               labels=c("AA", "AD", "DD"))

plot2 = ggplot(box2, aes(x=as.numeric(X1), group=as.factor(X2)), fill='green') +
        geom_boxplot(fill="3")+ xlab("PSI") +  ylab("Genotype") + ggtitle("Derived Lower Inclusion <50% mean PSI") + gtex_v8_figure_theme() + coord_flip()

plot3 = ggplot(to_draw[to_draw$pair == "1 vs 2" & to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) +  ylab("Enrichment in group 2                   Enrichment in group 1\n log(odd_ratio)")+ coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("Derived Higher Inclusion <50% mean PSI \nvs Derived Lower Inclusion <50% mean PSI") + gtex_v8_figure_theme() + guides(alpha = FALSE)

plot4 = ggplot(box3, aes(x=as.numeric(X1), group=as.factor(X2)), fill='blue') +
        geom_boxplot(fill="4")+ xlab("PSI") +  ylab("Genotype")+ ggtitle("Derived Higher Inclusion >50% mean PSI") + gtex_v8_figure_theme() + coord_flip() + scale_x_discrete(breaks=c("-0.2","0.0","0.2"),
                                                                                                                                                                              labels=c("AA", "AD", "DD"))

plot5 = ggplot(box4, aes(x=as.numeric(X1), group=as.factor(X2)), fill='violet') +
        geom_boxplot(fill="5")+ xlab("PSI") +  ylab("Genotype")+ ggtitle("Derived Lower Inclusion >50% mean PSI") + gtex_v8_figure_theme() + coord_flip()


plot6 = ggplot(to_draw[to_draw$pair == "3 vs 4" & to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) + ylab("Enrichment in group 4                   Enrichment in group 3\n log(odd_ratio)")+ coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("Derived Higher Inclusion >50% mean PSI \nvs Derived Lower Inclusion >50% mean PSI") + gtex_v8_figure_theme() + guides(alpha = FALSE)

plot = arrangeGrob(plot1, plot2, plot3, plot4, plot5, plot6, widths = c(1, 1, 1),
             layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)))

ggsave("fancy_statistical_summary_fisher_groups_all.png", plot=plot, path = "Data/visuals/", height = 5.11, width = 12.0,device='png', dpi=700)



dev.off()



plot1 = ggplot(box1, aes(x=as.numeric(X1), group=as.factor(X2))) +
        geom_boxplot(fill="2") + xlab("PSI") + ylab("Genotype") + ggtitle("Derived Higher Inclusion <50% mean PSI") + gtex_v8_figure_theme() + coord_flip()

plot2 = ggplot(box2, aes(x=as.numeric(X1), group=as.factor(X2)), fill='green') +
        geom_boxplot(fill="3") + xlab("PSI") +  ylab("Genotype") + ggtitle("Derived Lower Inclusion <50% mean PSI") + gtex_v8_figure_theme() + coord_flip()

plot3 = ggplot(to_draw[to_draw$pair == "1 vs 2" & to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + ylab("Enrichment in group 2                   Enrichment in group 1\n m1 - m2") + coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("Derived Higher Inclusion <50% mean PSI \nvs Derived Lower Inclusion <50% mean PSI") + gtex_v8_figure_theme() + guides(alpha = FALSE)

plot4 = ggplot(box3, aes(x=as.numeric(X1), group=as.factor(X2)), fill='blue') +
        geom_boxplot(fill="4")+ xlab("PSI") +  ylab("Genotype")+ ggtitle("Derived Higher Inclusion >50% mean PSI") + gtex_v8_figure_theme() + coord_flip() + scale_x_discrete(breaks=c("-0.2","0.0","0.2"),
                                                                                                                                                                              labels=c("AA", "AD", "DD"))

plot5 = ggplot(box4, aes(x=as.numeric(X1), group=as.factor(X2)), fill='violet') +
        geom_boxplot(fill="5")+ xlab("PSI") +  ylab("Genotype")+ ggtitle("Derived Lower Inclusion >50% mean PSI") + gtex_v8_figure_theme() + coord_flip()


plot6 = ggplot(to_draw[to_draw$pair == "3 vs 4" & to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) + ylab("Enrichment in group 4                   Enrichment in group 3\n m1 - m2") +coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("Derived Higher Inclusion >50% mean PSI \nvs Derived Lower Inclusion >50% mean PSI") + gtex_v8_figure_theme() + guides(alpha = "none")

plot = arrangeGrob(plot1, plot2, plot3, plot4, plot5, plot6, widths = c(1, 1, 1),
                    layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)))

ggsave("fancy_statistical_summary_ks_groups_all.png", plot=plot, path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)

View(to_draw)




library(diptest)
hist(data1$LENGTH)
dip.test(data1$LENGTH)
dip.test(data1$MIN_pLLDT)
# 
# 
# pre_data_1 = data1
# pre_data_1$LENGTH = data1$LENGTH
# #pre_data_1$LENGTH.ALIGN = data$LENGTH.ALIGN
# pre_data_1$STRUC = as.numeric(data1$GENE.NAME %in% data1_struc$GENE.NAME)
# pre_data_1$UNSTRUC = as.numeric(data1$GENE.NAME %in% data1_unstruc$GENE.NAME)
# pre_data_1$TRANSMEMBRANE = as.numeric(data1$TRANSMEMBRANE == "True")
# pre_data_1$HELIX = as.numeric(data1$HELIX == "True")
# pre_data_1$SHEET = as.numeric(data1$SHEET == "True")
# pre_data_1$TURN = as.numeric(data1$TURN == "True")
# pre_data_1$ASN = data1$ASN../data1$LENGTH *100
# pre_data_1$CYS = data1$CYS../data1$LENGTH *100
# pre_data_1$SIGNAL = as.numeric(data1$SIGNAL != "[]" & data1$SIGNAL != "['']" & data1$SIGNAL != "['', '']")
# pre_data_1$GROUP = replicate(nrow(pre_data_1), "1")
# 
# 
# pre_data_2 = data2
# pre_data_2$LENGTH = data2$LENGTH
# #pre_data_2$LENGTH.ALIGN = data2$LENGTH.ALIGN
# pre_data_2$STRUC = as.numeric(data2$GENE.NAME %in% data2_struc$GENE.NAME)
# pre_data_2$UNSTRUC = as.numeric(data2$GENE.NAME %in% data2_unstruc$GENE.NAME)
# pre_data_2$TRANSMEMBRANE = as.numeric(data2$TRANSMEMBRANE == "True")
# pre_data_2$HELIX = as.numeric(data2$HELIX == "True")
# pre_data_2$SHEET = as.numeric(data2$SHEET == "True")
# pre_data_2$TURN = as.numeric(data2$TURN == "True")
# pre_data_2$ASN = data2$ASN../data2$LENGTH *100
# pre_data_2$CYS = data2$CYS../data2$LENGTH *100
# pre_data_2$SIGNAL = as.numeric(data2$SIGNAL != "[]" & data2$SIGNAL != "['']" & data2$SIGNAL != "['', '']")
# pre_data_2$GROUP = replicate(nrow(pre_data_2), "2")
# #pre_data_2 = sample_n(pre_data_2, 200)
# 
# 
# pre_data_3 = data3
# pre_data_3$LENGTH = data3$LENGTH
# #pre_data_2$LENGTH.ALIGN = data2$LENGTH.ALIGN
# pre_data_3$STRUC = as.numeric(data3$GENE.NAME %in% data3_struc$GENE.NAME)
# pre_data_3$UNSTRUC = as.numeric(data3$GENE.NAME %in% data3_unstruc$GENE.NAME)
# pre_data_3$TRANSMEMBRANE = as.numeric(data3$TRANSMEMBRANE == "True")
# pre_data_3$HELIX = as.numeric(data3$HELIX == "True")
# pre_data_3$SHEET = as.numeric(data3$SHEET == "True")
# pre_data_3$TURN = as.numeric(data3$TURN == "True")
# pre_data_3$ASN = data3$ASN../data3$LENGTH *100
# pre_data_3$CYS = data3$CYS../data3$LENGTH *100
# pre_data_3$SIGNAL = as.numeric(data3$SIGNAL != "[]" & data3$SIGNAL != "['']" & data3$SIGNAL != "['', '']")
# pre_data_3$GROUP = replicate(nrow(pre_data_3), "3")
# 
# 
# pre_data_4 = data4
# pre_data_4$LENGTH = data4$LENGTH
# #pre_data_4$LENGTH.ALIGN = data4$LENGTH.ALIGN
# pre_data_4$STRUC = as.numeric(data4$GENE.NAME %in% data4_struc$GENE.NAME)
# pre_data_4$UNSTRUC = as.numeric(data4$GENE.NAME %in% data4_unstruc$GENE.NAME)
# pre_data_4$TRANSMEMBRANE = as.numeric(data4$TRANSMEMBRANE == "True")
# pre_data_4$HELIX = as.numeric(data4$HELIX == "True")
# pre_data_4$SHEET = as.numeric(data4$SHEET == "True")
# pre_data_4$TURN = as.numeric(data4$TURN == "True")
# pre_data_4$ASN = data4$ASN../data4$LENGTH *100
# pre_data_4$CYS = data4$CYS../data4$LENGTH *100
# pre_data_4$SIGNAL = as.numeric(data4$SIGNAL != "[]" & data4$SIGNAL != "['']" & data4$SIGNAL != "['', '']")
# pre_data_4$GROUP = replicate(nrow(pre_data_4), "4")
# data_for_pca = rbind(pre_data_1, pre_data_2, pre_data_3, pre_data_4)
# #data_for_pca = na.omit(data_for_pca)
# View(data_for_pca)
# data_for_pca = subset(data_for_pca, select=-c(anc_psi, der_psi, TOPOLOGY, DOMAIN, ASN.., CYS.., ALIGN.COORDS, LENGTH.ALIGN))
# ncol(data_for_pca)
# data_for_pca = data_for_pca[!is.na(data_for_pca$MEAN.ACC), ]
# data_for_pca.pca = prcomp(data_for_pca[, 75:96], center = T, scale. = T)
# data_for_pca.pca
# View(data_for_pca)
# library(dplyr)
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
# k <- kmeans(comp, 4, nstart=25, iter.max=1000)
# library(RColorBrewer)
# library(scales)
# palette(alpha(brewer.pal(9,'Set1'), 0.5))
# plot(comp, col=k$clust, pch=16)
# comp$STRUC = data_for_pca$STRUC
# comp$UNSTRUC = data_for_pca$UNSTRUC
# comp$GROUP = data_for_pca$GROUP
# ggplot(comp, aes(x=PC1, y=PC2)) +
#         geom_point(aes(color=as.factor(k$cluster), shape=as.factor(GROUP)), size=2.5) + scale_shape_manual(values=c(3, 25, 16, 11))
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
# mean1 = mean1/max1
# mean2 = mean2/max2
# data_mean = rbind(mean1, mean2, mean3)
# data_mean
# colnames(data_mean) = colnames(pre_data_1)[1:(ncol(pre_data_1)-1)]
# data_mean
# ggplot(data_mean, aes(x=replicate(3, "MIN"), y=MIN)) + geom_col(fill=c(1,2,3), position = position_dodge())  +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("-log10(p_value)") + labs(fill='Pairs') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)))
# #ggsave("statistical_summary_with_correction.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
# 
