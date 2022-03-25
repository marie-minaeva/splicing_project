library(dplyr)
library(stringr)
library(multimode)
library(ggplot2)
library(tidyverse)

#one column - one exon
fun = function(x) {
        data = readLines(x)
        data = unlist(str_split(data, " "))
        data = as.numeric(data)
        data = as.numeric(c(data, rep(NA, 1862 - length(data))))
        # if (length(data) >= 80){
        #         data = as.numeric(c(data, rep(NA, 1862 - length(data))))}else{
        #                 data = rep(NA, 1862)
        #         }
}

gtex_v8_figure_theme <- function() {
        return(theme(plot.title = element_text(face="plain",size=14), text = element_text(size=10),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}

plot_pLDDT = function(data, files, path) {
        for (i in 1:ncol(data)){
                f = unlist(str_split(as.character(files[i]), pattern = '.txt'))[1]
                f = unlist(str_split(f, pattern = '//'))[2]
                ggplot(data = data, aes(x=1:nrow(data), y= data[, i])) + geom_col() + gtex_v8_figure_theme() + xlab("position") + ylab("pLDDT")
                ggsave(paste(f, ".png"), path = path, height = 5.11, width = 7.92,device='png', dpi=150)
 
        }
}


setwd("~/splicing_project/")



non_coloc_data = c(0)
files <- list.files(path="Data/pLDDTs_non_coloc/", pattern="*.txt", full.names=TRUE, recursive=FALSE)
non_coloc_data = lapply(files, fun)
non_coloc_data = data.frame(non_coloc_data)



# non_coloc_data %>% select(where(not_all_na)) -> non_coloc_data
non_coloc_data = data.frame(non_coloc_data)
non_coloc_data = t(non_coloc_data)
non_coloc_data = cbind(non_coloc_data, "non_coloc")
non_coloc_data = data.frame(non_coloc_data)
View(non_coloc_data)




data = rbind(non_coloc_data)



#data = coloc_data
View(data)
labels = data[,1863]
data = data[, -1863]
#data = data[, 1:80]
data = t(data)

data = data.frame(data)
data=lapply(data,as.numeric)
data = data.frame(data)
colnames(data) = c()
View(data)

pval_non_coloc = c()
nmod_non_coloc = c()
for (i in 1:ncol(data)){
        pval_non_coloc = c(pval_non_coloc, try(modetest(data[,i])$p.value))
        nmod_non_coloc = c(nmod_non_coloc, try(nmodes(data[,i], bw=10.0)))
}
length(pval_non_coloc)






coloc_data = c(0)
files <- list.files(path="Data/pLDDTs_coloc/", pattern="*.txt", full.names=TRUE, recursive=FALSE)
coloc_data = lapply(files, fun)
coloc_data = data.frame(coloc_data)


# not_all_na <- function(x) all(!is.na(x))
# coloc_data %>% select(where(not_all_na)) -> coloc_data
coloc_data = data.frame(coloc_data)
coloc_data = t(coloc_data)
coloc_data = cbind(coloc_data, "coloc")
coloc_data = data.frame(coloc_data)
View(coloc_data)

data = rbind(coloc_data)



#data = coloc_data
View(data)
labels = data[,1863]
data = data[, -1863]
#data = data[, 1:80]
data = t(data)

data = data.frame(data)
data=lapply(data,as.numeric)
data = data.frame(data)
colnames(data) = c()
View(data)


pval_coloc = c()
nmod_coloc = c()
for (i in 1:ncol(data)){
        pval_coloc = c(pval_coloc, try(modetest(data[,i])$p.value))
        nmod_coloc = c(nmod_coloc, try(nmodes(data[,i], bw=10.0)))
}
length(pval_coloc)


pval_coloc = data.frame(pval_coloc)
pval_non_coloc = data.frame(pval_non_coloc)
unim_p = data.frame("coloc" = c(length(pval_coloc[pval_coloc$pval_coloc <= 0.05,]), length(pval_coloc[pval_coloc$pval_coloc > 0.05,])),
                  "non_coloc" = c(length(pval_non_coloc[pval_non_coloc$pval_non_coloc <= 0.05,]), length(pval_non_coloc[pval_non_coloc$pval_non_coloc > 0.05,])),
                  #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("multimodal", "unimodal"),
                  stringsAsFactors = FALSE)

unim_p
test = fisher.test(unim_p)
test




nmod_coloc = data.frame(nmod_coloc)
nmod_non_coloc = data.frame(nmod_non_coloc)

unim = data.frame("coloc" = c(length(nmod_coloc[nmod_coloc$nmod_coloc == 2,]), length(nmod_coloc[nmod_coloc$nmod_coloc == 1,])),
                  "non_coloc" = c(length(nmod_non_coloc[nmod_non_coloc$nmod_non_coloc == 2,]), length(nmod_non_coloc[nmod_non_coloc$nmod_non_coloc == 1,])),
                  #"non_sQTL" = c(length(data3[data2$TURN == "True",]$TURN), length(data3[data2$TURN != "True",]$TURN)),
                  row.names = c("multimodal", "unimodal"),
                  stringsAsFactors = FALSE)

unim
test = fisher.test(unim)
test



line = replicate(length(pval), 0.05)
ggplot(pval_coloc, aes(round(pval_coloc, 1))) +
        geom_bar(position = position_dodge2()) +
        geom_vline(xintercept = line, colour="red") +
        geom_text(stat='count', aes(label=..count..), vjust=0)
ggsave('coloc_pval_multimod_mass_excess.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=150)



ggplot(pval_non_coloc, aes(round(pval_non_coloc, 1))) +
        geom_bar(position = position_dodge2()) +
        geom_vline(xintercept = line, colour="red") +
        geom_text(stat='count', aes(label=..count..), vjust=0)
ggsave('non_coloc_pval_multimod_mass_excess.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=150)


ggplot(data.frame(nmod_coloc), aes(nmod_coloc)) +
        geom_bar(position = position_dodge2()) +
        geom_text(stat='count', aes(label=..count..), vjust=0)
ggsave('coloc_nmods_multimod_mass_excess.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=150)


ggplot(data.frame(nmod_non_coloc), aes(nmod_non_coloc)) +
        geom_bar(position = position_dodge2()) +
        geom_text(stat='count', aes(label=..count..), vjust=0)
ggsave('non_coloc_nmods_multimod_mass_excess.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=150)


# plot_pLDDT(data, files, "Data/visuals/non_coloc_pLDDT_position/")

# data = data.frame(data)
# summary(data)
# data[! is.numeric(data)]
# View(data)
# data.pca = prcomp(data, center = F)
# library(ggfortify)
# autoplot(data.pca, data = data)
# 
# 
# k = kmeans(na.omit(data), 2, )
# ggplot(data.pca, aes(x=PC1, y=PC2)) +
#         geom_point(aes(size=k$cluster, color=as.factor(labels)), alpha=0.5) + 
#         scale_shape_manual(values=c(3, 25, 16)) + gtex_v8_figure_theme() + ggtitle("length = 80")
# ggsave('coloc_non_coloc_kmean_pLDDT.png', path = "Data/visuals/", height = 5.11, width = 7.92,device='png', dpi=150)
# 
