library(dplyr)
library(ggplot2)
library(reshape2)

setwd("~/splicing_project/")
data1 = read.csv("Data/output_top_sQTL.csv")
data2 = read.csv("Data/output_non_sQTL_full.csv")
data_full = rbind(data1, data2)
ncol(data1)
ncol(data2)
data1 = data1[, colnames(data1) %in% colnames(data2)]
data2 = data2[, colnames(data2) %in% colnames(data1)]
ncol(data1)
ncol(data2)
data_full = rbind(data1, data2)

data_full %>% select(where(is.numeric)) -> data_for_cor
data_for_cor$CYS.. = data_for_cor$CYS.. / data_for_cor$LENGTH
data_for_cor$ASN.. = data_for_cor$ASN.. / data_for_cor$LENGTH
cor_mat = round(cor(data_for_cor, use='na.or.complete'), 2)

melted_cormat <- melt(cor_mat)
# png(file="~/splicing_project/Data/visuals/feature_correlation_heatmap.png", 1024, 1024, res=100, pointsize = 10)
# heatmap(cor_mat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile(color = "white") + 
        geom_text(aes(Var2, Var1, label = value), color = "black", size = 2.5) +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name="Pearson\nCorrelation") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size =10, hjust = 1), 
              axis.text.y = element_text(vjust = 1, 
                                         size =10, hjust = 1),
              axis.title.x = element_blank(), axis.title.y = element_blank()) + 
        coord_fixed()

ggsave('feature_correlation_heatmap.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=700)
