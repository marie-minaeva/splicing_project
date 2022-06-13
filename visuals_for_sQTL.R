library(ggplot2)
library(stats)
library(dplyr)
library(mltools)
library(foreach)
library(doParallel)
library(doSNOW)

setwd("~/splicing_project/")

source("NiceFigures.R")
source("MedianBootstrap.R")
source("FlagOutlier.R")








data = read.csv("Data/combined_top_sQTL_with_stop.csv")
data2 = read.csv("Data/output_non_sQTL.csv")

data = unique.data.frame(data)
data2 = unique.data.frame(data2)


data = data[!flag_outliers(data$LENGTH) & data$LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$LENGTH) & data2$LENGTH >= 15,]
nrow(data)
nrow(data2)


data = data[colnames(data) %in% colnames(data2)]
data2 = data2[colnames(data2) %in% colnames(data)]
data$GROUP = "sQTL"
data2$GROUP = "non_sQTL"
to_draw = rbind(data, data2)
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
pair = c()
trait = c()
conf_low = c()
conf_top = c()
stat = c()
test_type = c()

#LENGTH TESTS
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

#Q3_pLDDT TESTS

trait = c(trait, "Q3_pLDDT")
pair = c(pair, "sQTL vs non_sQTL")
test_type = c(test_type, "ks")
out = ks_stat(data$Q3_pLLDT, data2$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])


# SYM TESTS
symet = data.frame("sQTL" = c(length(data[data$LENGTH %% 3 == 0,]$LENGTH), length(data[data$LENGTH %% 3 != 0,]$LENGTH)),
                   "non_sQTL" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
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


#TRANSMEMBRANE TEST

trans= data.frame("sQTL" = c(length(data[data$TRANSMEMBRANE == "True",]$LENGTH), length(data[data$TRANSMEMBRANE != "True",]$LENGTH)),
                  "non_sQTL" = c(length(data2[data2$TRANSMEMBRANE == "True",]$LENGTH), length(data2[data2$TRANSMEMBRANE != "True",]$LENGTH)),
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

station = data.frame("sQTL" = c(sum(is_domain_data$SUM), nrow(is_domain_data) - sum(is_domain_data$SUM)),
                     "non_sQTL" = c(sum(is_domain_data2$SUM), nrow(is_domain_data2) - sum(is_domain_data2$SUM)),
                     row.names = c("stat", "non_stat"),
                     stringsAsFactors = FALSE)
print(station)
test = fisher.test(station)
test
pval = c(pval, test$p.value)
trait = c(trait, "Domain")
pair = c(pair, "sQTL vs non_sQTL")
conf_low = c(conf_low, test$conf.int[1])
conf_top = c(conf_top, test$conf.int[2])
stat = c(stat, test$estimate)
test_type = c(test_type, "fisher")


pval = p.adjust(pval, method="hochberg")

to_draw = data.frame()
to_draw = cbind(trait, log10(pval))
to_draw = cbind(to_draw, pair)
to_draw = cbind(to_draw, conf_low)
to_draw = cbind(to_draw, conf_top)
to_draw = cbind(to_draw, stat)
to_draw = cbind(to_draw, test_type)
colnames(to_draw) = c("trait", "pval", "pair", "conf_low", "conf_top", "statistics", "test_type")
to_draw = data.frame(to_draw)
to_draw$pval = -as.numeric(to_draw$pval)
line = replicate(length(pval), 1.5)
ggplot(to_draw, aes(fill = as.factor(pair), x=trait, y=pval )) + geom_col(position = position_dodge())  + geom_hline(yintercept = line) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), legend.text = element_text(size=7), axis.title.y = element_text(size = 8), axis.title.x=element_blank(), legend.position = c(0.9, 0.8)) + ylab("-log10(p_value)") + labs(fill='Pairs') + guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5))) + gtex_v8_figure_theme()
ggsave('statistical_summary_without_correction_sQTL_non_sQTL.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
to_draw$conf_low = as.numeric(to_draw$conf_low)
to_draw$conf_top = as.numeric(to_draw$conf_top)

write.csv(to_draw, "to_draw.csv", row.names = F)
line1 = replicate(length(pval), 0.0)
to_draw$alpha = ifelse((to_draw$test_type == "ks"  & ((to_draw$conf_low <= 0 & to_draw$conf_top >= 0) | 
                                                        (to_draw$conf_low >= 0 & to_draw$conf_top <= 0))) | 
                         (to_draw$test_type == "fisher"  & ((log(to_draw$conf_low) <= 0 & log(to_draw$conf_top) >= 0) | 
                                                              (log(to_draw$conf_low) >= 0 & log(to_draw$conf_top) <= 0))), yes = 0.4, 1.0)


ggplot(to_draw[to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha)) +
        geom_pointrange(aes(x = trait, y = log(as.numeric(statistics)), ymax = log(as.numeric(conf_top)), ymin = log(as.numeric(conf_low)))) +
        ylab("Enrichment in non sQTLs                                                                  Enrichment in sQTLs\n log(odd_ratio)") + 
  geom_hline(yintercept = line1, color="red") + guides(alpha = FALSE) + gtex_v8_figure_theme()
ggsave('statistical_summary_fisher_sQTL_non_sQTL.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)

line1 = replicate(length(pval), 0.0)
ggplot(to_draw[to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics))) +
        geom_pointrange(aes(x = trait, y = as.numeric(statistics), ymax = as.numeric(conf_top), ymin = as.numeric(conf_low))) +
        ylab("Enrichment in non sQTLs                                                                  Enrichment in sQTLs\n m1-m2") + 
  geom_hline(yintercept = line1, color="red") + gtex_v8_figure_theme() + guides(alpha = FALSE)


ggsave('statistical_summary_ks_sQTL_non_sQTL.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)




write.csv(to_draw, "Data/to_draw_sQTL.csv", quote=F, row.names = F)

