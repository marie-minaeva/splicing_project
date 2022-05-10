setwd("~/splicing_project/")

library(ggplot2)
library(stats)
library(dplyr)
library(boot)
library(mltools)
library(foreach)
library(doParallel)
library(doSNOW)
library(ggpubr)

gtex_v8_figure_theme <- function() {
        return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_sQTL = read.csv("Data/to_draw_sQTL.csv")
data_inc = read.csv("Data/to_draw_inclusion_groups.csv")
data_inc_over = read.csv("Data/to_draw_inclusion_groups_overall.csv")
data_psi = read.csv("Data/to_draw_psi_groups.csv")
data_rand = read.csv("Data/to_draw_inclusion_groups_random_shuffle.csv")

colnames(data_sQTL)
colnames(data_psi)
colnames(data_inc)
colnames(data_inc_over)
cols = c("trait", "pval", "conf_low", "conf_top", "statistics", "test_type", "alpha") 
data_sQTL = data_sQTL[cols]
data_inc = data_inc[cols]
data_psi = data_psi[cols]
data_rand = data_rand[cols]
data_inc_over = data_inc_over[cols]

data_sQTL$group = "sQTLs vs. non_sQTLs"
data_psi$group = "Increasing vs. decreasing"
data_inc$group = "Highly sQTLs vs. lowly sQTLs"
data_inc_over$group = "Highly vs. lowly"
data_rand$group = "Random Shuffle"
to_draw = rbind(data_sQTL, data_inc, data_psi, data_rand, data_inc_over)
line1 = replicate(length(pval), 0.0)
to_draw = to_draw[to_draw$trait != "Helix" & to_draw$trait != "Turn" & to_draw$trait != "Sheet",]

ggplot(to_draw[to_draw$test_type == "fisher" , ], aes(x = trait, y = log(as.numeric(statistics)), alpha = alpha, color=group, group=group)) +
        geom_pointrange(aes(ymin = log(conf_low), ymax = log(conf_top)), position=position_dodge2(width = 0.3)) +
        ylab("Enrichment in group 1                                                              Enrichment in group 2\n log(odd_ratio)") +
        coord_flip() + gtex_v8_figure_theme() + guides(alpha = FALSE) + 
        theme(legend.position = "bottom") + scale_alpha(range = c(0.25, 1.0))
ggsave("statistical_summary_fisher.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)


ggplot(to_draw[to_draw$test_type == "ks", ], aes(x = trait, y = as.numeric(statistics), alpha = alpha, color=group, group=group)) +
        geom_pointrange(aes(ymin = conf_low, ymax = conf_top), position=position_dodge2(width = 0.3)) +
        ylab("Enrichment in group 1                                                                 Enrichment in group 2\n m1 - m2") +  
        coord_flip()  + gtex_v8_figure_theme() + 
        guides(alpha = FALSE) + theme(legend.position = "bottom") + scale_alpha(range = c(0.25, 1.0))
ggsave("statistical_summary_ks.png", path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
