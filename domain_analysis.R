setwd("~/splicing_project/")


library(ggplot2)
library(stringr)



data_non = read.csv("Data/output_non_sQTL_full.csv")
data = read.csv("Data/output_data_coloc.csv")
features = colnames(data_non)[50:ncol(data_non)]
features = features[colnames(data_non)[50:ncol(data_non)] %in% colnames(data)]
features
enrich = data.frame()
for (feature in features){
        data_non[,feature] = sapply(data_non[,feature],function(x){
                return(unlist(str_split(x, pattern = ";"))[1])
        })
        data[,feature] = sapply(data[,feature],function(x){
                return(unlist(str_split(x, pattern = ";"))[1])
        })
        data_non[,feature] = ifelse(data_non[,feature] == '', NA, data_non[,feature])
        data[,feature] = ifelse(data[,feature] == '', NA, data[,feature])
        tab_non = data.frame(table(data_non[!is.na(data_non[,feature]),feature]))
        tab_non$dataset = "non_coloc"
        
        tab = data.frame(table(data[!is.na(data[,feature]),feature]))
        try(tab$dataset <- "coloc")
        station = data.frame("coloc" = c(sum(tab$Freq), length(data$GENE.NAME) - sum(tab$Freq)),
                             "non_coloc" = c(sum(tab_non$Freq), length(data_non$GENE.NAME) - sum(tab_non$Freq)),
                             row.names = c("stat", "non_stat"),
                             stringsAsFactors = FALSE)
        print(feature)
        print(station)
        test = fisher.test(station)
        enrich = rbind(enrich, c(feature, test$p.value, test$estimate))
        tab = rbind(tab, tab_non)
        
        tab = tab[tab$Freq > 1,]
        tab[tab$dataset == "coloc",]$Freq = tab[tab$dataset == "coloc",]$Freq / sum(tab[tab$dataset == "coloc",]$Freq)
        tab[tab$dataset == "non_coloc",]$Freq = tab[tab$dataset == "non_coloc",]$Freq / sum(tab[tab$dataset == "non_coloc",]$Freq)
        p = ggplot(tab, aes(x = dataset, fill=Var1, y=Freq)) + 
                geom_bar(width = 1, stat = "identity", position = position_stack()) +
                xlab("") + theme(legend.position="bottom") + 
                labs(fill=paste(feature,"with \nmore than 1 evidence", sep=" "))
        if (feature == "REGION"){
                p = ggplot(tab[tab$Var1 != "Disordered", ], aes(x = dataset, fill=Var1, y=Freq)) + 
                        geom_bar(width = 1, stat = "identity", position = position_stack()) +
                        xlab("") + theme(legend.position="bottom") + 
                        labs(fill=paste(feature,"with \nmore than 1 evidence", sep=" ")) + 
                        theme(legend.key.size = unit(0.2, 'cm'), legend.text=element_text(size=3))
        } 
        ggsave(p, filename = paste0("Data/visuals/", feature, "_comparison_coloc_non_coloc.png"), width = 12, height = 4)
        
}
colnames(enrich) = c("feature", "p.value", "odds.ratio")
View(enrich)
