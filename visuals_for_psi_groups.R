library(ggplot2)
library(gridExtra)
library(stats)
library(dplyr)
library(boot)
library(mltools)
library(foreach)
library(doParallel)
library(doSNOW)
library(ggpubr)
library(ggExtra)


setwd("~/splicing_project/")

source("MedianBootstrap.R")
source("FlagOutlier.R")
source("NiceFigures.R")

data_full = read.csv("Data/combined_top_sQTL_with_stop.csv")


data_full = data_full[!is.na(data_full$anc_allele_freq), ]
anc_psi = c(data_full[data_full$anc_allele == data_full$ref_allele, ]$mean_00_psi, data_full[data_full$anc_allele == data_full$alt_allele, ]$mean_11_psi)
der_psi = c(data_full[data_full$anc_allele == data_full$ref_allele, ]$mean_11_psi, data_full[data_full$anc_allele == data_full$alt_allele, ]$mean_00_psi)
data_full$anc_psi = anc_psi
data_full$der_psi = der_psi
data1 = data_full[(data_full$anc_psi < data_full$mean_01_psi) & (data_full$mean_01_psi < data_full$der_psi), ]
data2 = data_full[(data_full$anc_psi > data_full$mean_01_psi) & (data_full$mean_01_psi > data_full$der_psi), ]

data1 = data[!flag_outliers(data1$LENGTH) & data1$LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$LENGTH) & data2$LENGTH >= 15,]

data1$GROUP = "Increasing"
data2$GROUP = "Decreasing"

data1 = data1[colnames(data1) %in% colnames(data2)]
data2 = data2[colnames(data2) %in% colnames(data1)]

to_draw = rbind(data1, data2)
to_draw$SYMM = to_draw$LENGTH %% 3

ggdensity(to_draw, x="SYMM", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")
ggdensity(to_draw, x="MIN", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")
ggdensity(to_draw, x="Q1", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")
ggdensity(to_draw, x="Q2", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")
ggdensity(to_draw, x="Q3", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")
ggdensity(to_draw, x="MAX", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")

ggdensity(to_draw, x="MIN_pLLDT", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")
ggdensity(to_draw, x="Q1_pLLDT", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")
ggdensity(to_draw, x="Q2_pLLDT", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")
ggdensity(to_draw, x="Q3_pLLDT", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")
ggdensity(to_draw, x="MAX_pLLDT", add = "median", rug = TRUE, color = "GROUP", fill = "GROUP")

ggboxplot(to_draw, x="GROUP", y="MIN", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()
ggboxplot(to_draw, x="GROUP", y="Q1", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()
ggboxplot(to_draw, x="GROUP", y="Q2", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()
ggboxplot(to_draw, x="GROUP", y="Q3", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()
ggboxplot(to_draw, x="GROUP", y="MAX", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()

ggboxplot(to_draw, x="GROUP", y="MIN_pLLDT", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()
ggboxplot(to_draw, x="GROUP", y="Q1_pLLDT", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()
ggboxplot(to_draw, x="GROUP", y="Q2_pLLDT", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()
ggboxplot(to_draw, x="GROUP", y="Q3_pLLDT", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()
ggboxplot(to_draw, x="GROUP", y="MAX_pLLDT", add = "jitter", rug = TRUE, color = "GROUP") +
        stat_compare_means()


pval = c()
trait = c()
conf_low = c()
conf_top = c()
stat = c()
test_type = c()

#LENGTH TESTS

trait = c(trait, "Length")
test_type = c(test_type, "ks")
out = ks_stat(data1$LENGTH, data2$LENGTH)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


#ASN... TESTS

trait = c(trait, "% ASN")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$ASN../data1$LENGTH *100),]$ASN../data1$LENGTH*100, 
              data2[!flag_outliers(data2$ASN../data2$LENGTH *100),]$ASN../data2$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])






#CYS.. TESTS

trait = c(trait, "% CYS")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$CYS../data1$LENGTH *100),]$CYS../data1$LENGTH*100, 
              data2[!flag_outliers(data2$CYS../data2$LENGTH *100),]$CYS../data2$LENGTH*100)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


#Q1 TESTS

trait = c(trait, "Q1_RSA")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$Q1),]$Q1, data2[!flag_outliers(data2$Q1),]$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


#Q3_pLDDT TESTS

trait = c(trait, "Q3_pLDDT")
test_type = c(test_type, "ks")
out = ks_stat(data1[!flag_outliers(data1$Q3_pLLDT),]$Q3_pLLDT,
              data2[!flag_outliers(data2$Q3_pLLDT),]$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


# SYM TESTS
symet = data.frame("sQTL" = c(length(data1[data1$LENGTH %% 3 == 0,]$LENGTH), length(data1[data1$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_sQTL" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
                   row.names = c("sym", "non_sym"),
                   stringsAsFactors = FALSE)

symet
test = fisher.test(symet)
pval = c(pval, test$p.value)
trait = c(trait, "Symmetric")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")


#TRANSMEMBRANE TEST

trans= data.frame("sQTL" = c(length(data1[data1$TRANSMEMBRANE == "True",]$LENGTH), length(data1[data1$TRANSMEMBRANE != "True",]$LENGTH)),
                  "non_sQTL" = c(length(data2[data2$TRANSMEMBRANE == "True",]$LENGTH), length(data2[data2$TRANSMEMBRANE != "True",]$LENGTH)),
                  row.names = c("trans", "non_trans"),
                  stringsAsFactors = FALSE)
trans
test = fisher.test(trans)
pval = c(pval, test$p.value)
trait = c(trait, "Transmembrane")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")


# Domain analysis
features = c("DOMAIN", "TRANSMEMBRANE", "MOTIF", "TOPO_DOM", "ACT_SITE", "MOD_RES",
             "REGION", "REPEAT", "TRANSMEM", "NP_BIND", "DNA_BIND", "CROSSLNK", 
             "ZN_FING", "METAL", "SITE", "INTRAMEM", "LIPID")  
data[data == ""] = NA
data2[data2 == ""] = NA
is_domain_data = is.na(data[,features])
is_domain_data2 = is.na(data2[, features])
dim(is_domain_data2)
View(is_domain_data)
is_domain_data$SUM = rowSums(is_domain_data)
is_domain_data2$SUM = rowSums(is_domain_data2)
is_domain_data$SUM = ifelse(is_domain_data$SUM == length(features), F, T)
is_domain_data2$SUM = ifelse(is_domain_data2$SUM == length(features), F, T)


station = data.frame("sQTL" = c(sum(is_domain_data$SUM), nrow(is_domain_data) - sum(is_domain_data$SUM)),
                     "non_sQTL" = c(sum(is_domain_data2$SUM), nrow(is_domain_data2) - sum(is_domain_data2$SUM)),
                     row.names = c("stat", "non_stat"),
                     stringsAsFactors = FALSE)
print(station)
test = fisher.test(station)
est = fisher.test(station)
pval = c(pval, test$p.value)
trait = c(trait, "Domain")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")


pval = p.adjust(pval, method="hochberg")
to_draw = data.frame()
to_draw = cbind(trait, log10(pval))
to_draw = cbind(to_draw, conf_low)
to_draw = cbind(to_draw, conf_top)
to_draw = cbind(to_draw, stat)
to_draw = cbind(to_draw, test_type)
colnames(to_draw) = c("trait", "pval", "conf_low", "conf_top", "statistics", "test_type")
to_draw = data.frame(to_draw)
to_draw$pval = -as.numeric(to_draw$pval)
to_draw$group1 = as.factor(groups[seq(from=1, to=length(groups), by=2)])
to_draw$group2 = as.factor(groups[seq(from=2, to=length(groups), by=2)])
to_draw$conf_low = as.numeric(to_draw$conf_low)
to_draw$conf_top = as.numeric(to_draw$conf_top)
to_draw$alpha = ifelse((to_draw$test_type == "ks"  & ((to_draw$conf_low <= 0 & to_draw$conf_top >= 0) | 
                                                              (to_draw$conf_low >= 0 & to_draw$conf_top <= 0))) | 
                               (to_draw$test_type == "fisher"  & ((log(to_draw$conf_low) <= 0 & log(to_draw$conf_top) >= 0) | 
                                                                          (log(to_draw$conf_low) >= 0 & log(to_draw$conf_top) <= 0))), yes = 0.4, 1.0)



line = replicate(length(pval), 1.5)



gghistogram(to_draw, x = "trait", y="pval", rug = F,  stat="identity", position = "dodge")  + 
        geom_hline(yintercept = line, linetype = "dashed") +
        theme(axis.text.x = element_text(angle = 90))


ggsave("statistical_summary_all_groups.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)




line1 = replicate(length(pval), 0.0)

ggplot(to_draw[to_draw$test_type == "fisher" , ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_pointrange(aes(ymin = log(conf_low), ymax = log(conf_top))) +
        ylab("Enrichment in group decreasing                                                                Enrichment in group increasing\n log(odd_ratio)")+
        coord_flip()  + geom_hline(yintercept = line1, color="red")  + gtex_v8_figure_theme() + guides(alpha = FALSE)
ggsave("statistical_summary_fisher_groups_inc_dec.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)


ggplot(to_draw[to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha)) +
        geom_pointrange(aes(ymin = log(conf_low), ymax = log(conf_top))) +
        ylab("Enrichment in group decreasing                                                                 Enrichment in group increasing\n m1 - m2") +  coord_flip()  + geom_hline(yintercept = line1, color="red") + ggtitle("1 vs 2") + gtex_v8_figure_theme()
ggsave("statistical_summary_ks_groups_inc_dec.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)

write.csv(to_draw, "Data/to_draw_psi_groups.csv", quote=F, row.names = F)
