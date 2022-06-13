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
library(ggpubr)

source("NiceFigures.R")
source("FlagOutlier.R")
source("MedianBootstrap.R")

s_ks = c()
c_l_ks = c()
c_t_ks = c()
s_f = c()
c_l_f = c()
c_t_f = c()
data = read.csv("Data/output_data_coloc.csv")
data2 = read.csv("Data/output_data_non_coloc.csv")

data = unique.data.frame(data)
data2 = unique.data.frame(data2)

data = data[!flag_outliers(data$LENGTH) & data$LENGTH >= 15,]
data2 = data2[!flag_outliers(data2$LENGTH) & data2$LENGTH >= 15,]

data = data[colnames(data) %in% colnames(data2)]
data2 = data2[colnames(data2) %in% colnames(data)]
data$GROUP = "coloc"
data2$GROUP = "non_coloc"
ncol(data)
ncol(data2)
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
test_type = c()
conf_low = c()
conf_top = c()
stat = c()

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


#Q1 TESTS

trait = c(trait, "Q1_RSA")
pair = c(pair, "coloc vs non_coloc")
test_type = c(test_type, "ks")
out = ks_stat(data[!flag_outliers(data$Q1),]$Q1, data2[!flag_outliers(data2$Q1),]$Q1)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

#Q3_pLDDT TESTS

trait = c(trait, "Q3_pLDDT")
pair = c(pair, "coloc vs non_coloc")
test_type = c(test_type, "ks")
out = ks_stat(data[!flag_outliers(data$Q3_pLLDT),]$Q3_pLLDT, data2[!flag_outliers(data2$Q3_pLLDT),]$Q3_pLLDT)
pval = c(pval, out[1])
conf_low = c(conf_low, out[3])
conf_top = c(conf_top, out[4])
stat = c(stat, out[2])

# SYM TESTS
symet = data.frame("coloc" = c(length(data[data$LENGTH %% 3 == 0,]$LENGTH), length(data[data$LENGTH %% 3 != 0,]$LENGTH)),
 "non_coloc" = c(length(data2[data2$LENGTH %% 3 == 0,]$LENGTH), length(data2[data2$LENGTH %% 3 != 0,]$LENGTH)),
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

#TRANSMEMBRANE TEST

trans= data.frame("coloc" = c(length(data[data$TRANSMEMBRANE == "True",]$LENGTH), length(data[data$TRANSMEMBRANE != "True",]$LENGTH)),
 "non_coloc" = c(length(data2[data2$TRANSMEMBRANE == "True",]$LENGTH), length(data2[data2$TRANSMEMBRANE != "True",]$LENGTH)),
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

features = c("DOMAIN", "TRANSMEMBRANE", "MOTIF", "TOPO_DOM", "ACT_SITE", "MOD_RES",
 "REGION", "REPEAT", "TRANSMEM", "NP_BIND", "DNA_BIND", "CROSSLNK", 
 "ZN_FING", "METAL", "SITE", "INTRAMEM", "LIPID")
is_domain_data = is.na(data[,features])
is_domain_data2 = is.na(data2[, features])
dim(is_domain_data)
View(is_domain_data)
is_domain_data$SUM = rowSums(is_domain_data)
is_domain_data2$SUM = rowSums(is_domain_data2)
is_domain_data$SUM = ifelse(is_domain_data$SUM == 21, F, T)
is_domain_data2$SUM = ifelse(is_domain_data2$SUM == 21, F, T)


station = data.frame("coloc" = c(sum(is_domain_data$SUM), nrow(is_domain_data) - sum(is_domain_data$SUM)),
 "non_coloc" = c(sum(is_domain_data2$SUM), nrow(is_domain_data2) - sum(is_domain_data2$SUM)),
 row.names = c("stat", "non_stat"),
 stringsAsFactors = FALSE)
print(station)
test = fisher.test(station)
est = fisher.test(station)
pval = c(pval, test$p.value)
trait = c(trait, "Domain")
pair = c(pair, "coloc vs non_coloc")
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
to_draw
colnames(to_draw) = c("trait", "pval", "pair", "conf_low", "conf_top", "statistics", "test_type")
to_draw = data.frame(to_draw)
to_draw$pval = -as.numeric(to_draw$pval)
to_draw



line = replicate(length(pval), 1.5)
gghistogram(to_draw, x = "trait", y="pval", rug = F, color = "pair", fill = "pair", stat="identity", position = "dodge")+ 
geom_hline(yintercept = line, linetype = "dashed") +
theme(axis.text.x = element_text(angle = 90))
file = paste(c("statistical_summary_with_correction.png"), collapse="")
ggsave(file, path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)

to_draw$trait = factor(to_draw$trait, levels = unique(trait))

to_draw$conf_low = as.numeric(to_draw$conf_low)
to_draw$conf_top = as.numeric(to_draw$conf_top)

to_draw$alpha = ifelse((to_draw$test_type == "ks"& ((to_draw$conf_low <= 0 & to_draw$conf_top >= 0) | 
(to_draw$conf_low >= 0 & to_draw$conf_top <= 0))) | 
 (to_draw$test_type == "fisher"& ((log(to_draw$conf_low) <= 0 & log(to_draw$conf_top >= 0)) | 
(log(to_draw$conf_low) >= 0 & log(to_draw$conf_top) <= 0))), yes = 0.6, 1.0)
to_draw[to_draw$test_type == "fisher",]$alpha = ifelse((log(to_draw[to_draw$test_type == "fisher",]$conf_low) <= 0 & log(to_draw[to_draw$test_type == "fisher",]$conf_top) >= 0) | (log(to_draw[to_draw$test_type == "fisher",]$conf_low) >= 0 & log(to_draw[to_draw$test_type == "fisher",]$conf_top) <= 0), yes = 0.2, 1.0)



to_draw$trait = factor(to_draw$trait, levels = unique(trait))
line1 = replicate(length(pval), 0.0)

ggplot(to_draw[to_draw$test_type == "fisher", ], aes(x = trait, y = log(as.numeric(statistics)),alpha = alpha)) +
        geom_pointrange(size = 4, position=position_dodge2(width=0.9)) +
        ylab("Enrichment in non coloc Enrichment in coloc\n log(odd_ratio)") + 
        coord_flip() + geom_hline(yintercept = line1, color="red") + 
        theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "both")),
        axis.title.x = element_text(angle = 0)) + gtex_v8_figure_theme() + guides(alpha = FALSE)
ggsave('statistical_summary_fisher_coloc.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)

line1 = replicate(length(pval), 0.0)
ggplot(to_draw[to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics),alpha = alpha)) +
        geom_pointrange(size = 4, position=position_dodge2(width=0.9)) +
        ylab("Enrichment in non colocEnrichment in coloc\n m1 - m2") + 
        coord_flip() + geom_hline(yintercept = line1, color="red") + 
        theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "both")),
        axis.title.x = element_text(angle = 0)) + gtex_v8_figure_theme()+ guides(alpha = FALSE)
ggsave('statistical_summary_ks_coloc.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
