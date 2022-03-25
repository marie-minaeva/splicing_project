setwd("~/splicing_project/")


library(ggplot2)

data_non = read.csv("Data/output_non_sQTL_full.csv")
data_non$MOTIF = ifelse(data_non$MOTIF == '', NA, data_non$MOTIF)

data = read.csv("Data/output_top_sQTL.csv")
data$MOTIF = ifelse(data$MOTIF == '', NA, data$MOTIF)
# data_non$MOTIF = ifelse(length(data_non$MOTIF) > 20, data_non$MOTIF[1:20], data_non$MOTIF)
# ggplot(data_non[!is.na(data_non$MOTIF),], aes(x = "", fill=MOTIF)) + geom_bar(width = 1, position = position_fill()) +
#         xlab("") + theme(legend.position="bottom")
table(data_non$MOTIF)
tab_non = data.frame(table(data_non[!is.na(data_non$MOTIF),]$MOTIF))
tab_non$dataset = "non_sQTL"

tab = data.frame(table(data[!is.na(data$MOTIF),]$MOTIF))
tab$dataset = "sQTL"

tab = rbind(tab, tab_non)

tab = tab[tab$Freq > 1,]
tab[tab$dataset == "sQTL",]$Freq = tab[tab$dataset == "sQTL",]$Freq / sum(tab[tab$dataset == "sQTL",]$Freq)
tab[tab$dataset == "non_sQTL",]$Freq = tab[tab$dataset == "non_sQTL",]$Freq / sum(tab[tab$dataset == "non_sQTL",]$Freq)
p = ggplot(tab, aes(x = dataset, fill=Var1, y=Freq)) + 
        geom_bar(width = 1, stat = "identity") +
        xlab("") + theme(legend.position="bottom") + 
        labs(fill="Domains with \nmore than 1 evidence") 
        # geom_text(stat='identity', aes(label=Freq/sum(Freq)))
p
View(tab)

ggsave(p, filename = "Data/visuals/motif_comparison_sQTL_non_sQTL.png", width = 12, height = 4)
