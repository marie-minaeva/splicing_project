setwd("~/splicing_project/Data/pLDDTs/")
library(ggplot2)
library(stringr)
library(moments)
files = readLines("~/splicing_project/files.txt")
files
mea = numeric(0)
med = numeric(0)

entropy <- function(target) {
        freq <- table(target)/length(target)
        # vectorize
        vec <- as.data.frame(freq)[,2]
        #drop 0 to avoid NaN resulting from log2
        vec<-vec[vec>0]
        #compute entropy
        -sum(vec * log2(vec))
}

ent = numeric(0)
ent1 = numeric(0)
ent2 = numeric(0)
skew = numeric(0)
skew1 = numeric(0)
skew2 = numeric(0)

files1 = readLines("~/splicing_project/files1.txt")
files2 = readLines("~/splicing_project/files2.txt")
print(length(files))
for (file in files){
        file = "Q13530_all_pLDDT.txt"
        data = readLines(file)
        name = unlist(str_split(file, "_"))[1]
        data = as.numeric(unlist(str_split(data, " ")))
        mea = c(mea, mean(data))
        med = c(med, median(data))
        data = data.frame(data)
        ent = c(ent, entropy(data$data))
        skew = c(skew, skewness(data$data))
        # 
         ggplot() + geom_density(data=data, aes(data), alpha = 0.2)
         ggsave(paste(name, ".png", sep=''), height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/pLDDT_densities/", device='png', dpi=700)
}
setwd("~/splicing_project/Data/pLDDTs_non_coloc/")
for (file in files1){
        data = readLines(file)
        name = unlist(str_split(file, "_"))[1]
        data = as.numeric(unlist(str_split(data, " ")))
        #mea = c(mea, mean(data))
        #med = c(med, median(data))
        data = data.frame(data)
        ent1 = c(ent1, entropy(data$data))
        skew1 = c(skew1, skewness(data$data))
        # 
        # ggplot() + geom_density(data=data, aes(data), alpha = 0.2)
        # ggsave(paste(name, ".png", sep=''), height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/pLDDT_densities/", device='png', dpi=700)
}
setwd("~/splicing_project/Data/pLDDTs_non_sQTL/")
med2 = numeric(0)
mea2 = numeric(0)
for (file in files2){
        data = readLines(file)
        name = unlist(str_split(file, "_"))[1]
        data = as.numeric(unlist(str_split(data, " ")))
        mea2 = c(mea2, mean(data))
        med2 = c(med2, median(data))
        data = data.frame(data)
        ent2 = c(ent2, entropy(data$data))
        skew2 = c(skew2, skewness(data$data))
        # 
        ggplot() + geom_density(data=data, aes(data), alpha = 0.2)
        ggsave(paste(name, ".png", sep=''), height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/pLDDT_densities_non_sQTL/", device='png', dpi=700)
}
mea
mea = data.frame(mea)
ggplot() + geom_density(data=mea, aes(mea), alpha = 0.2) + xlab("pLDDT")
ggsave( "mean_pLDDT_density.png", height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/", device='png', dpi=700)


med
med = data.frame(med)
ggplot() + geom_density(data=med, aes(med), alpha = 0.2) + xlab("pLDDT")
ggsave( "median_pLDDT_density.png", height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/", device='png', dpi=700)

med
med = data.frame(med)
ggplot() + geom_density(data=med, aes(med, fill="median"), alpha = 0.2) + xlab("pLDDT") + geom_density(data=mea, aes(mea, fill="mean"), alpha = 0.2) 
ggsave( "mean_median_pLDDT_density.png", height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/", device='png', dpi=700)

mea2
mea2 = data.frame(mea2)
med2
med2 = data.frame(med2)
ggplot() + geom_density(data=med2, aes(med2, fill="median"), alpha = 0.2) + xlab("pLDDT") + geom_density(data=mea2, aes(mea2, fill="mean"), alpha = 0.2) 
ggsave( "mean_median_pLDDT_density_non_sQTL.png", height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/", device='png', dpi=700)

ggplot() + geom_density(data=med2, aes(med2, fill="median"), alpha = 0.2) + xlab("pLDDT") + geom_density(data=mea2, aes(mea2, fill="mean"), alpha = 0.2) + geom_density(data=med, aes(med, fill="median_coloc"), alpha = 0.2) +  geom_density(data=mea, aes(mea, fill="mean_coloc"), alpha = 0.2) 
ggsave( "mean_median_pLDDT_density_non_sQTL_coloc.png", height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/", device='png', dpi=700)

ks.test(mea$mea, mea2$mea2)
ks.test(med$med, med2$med2)

ent
ent = data.frame(ent)
ent1
ent1 = data.frame(ent1)
ent2
ent2 = data.frame(ent2)
ggplot() + geom_density(data=ent, aes(ent, fill="coloc"), alpha = 0.2) + xlab("ent") + geom_density(data=ent1, aes(ent1, fill='non_coloc'), alpha = 0.2) + geom_density(data=ent2, aes(ent2, fill='non_sQTL'), alpha = 0.2)
ggsave( "entropy_pLDDT_density.png", height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/", device='png', dpi=700)
ks.test(ent$ent, ent1$ent1)
ks.test(ent$ent, ent2$ent2)
ks.test(ent1$ent1, ent2$ent2)





skew
skew = data.frame(skew)
skew1
skew1 = data.frame(skew1)
skew2
skew2 = data.frame(skew2)
ggplot() + geom_density(data=skew, aes(skew, fill="coloc"), alpha = 0.2) + xlab("skew") + geom_density(data=skew1, aes(skew1, fill='non_coloc'), alpha = 0.2) + geom_density(data=skew2, aes(skew2, fill='non_sQTL'), alpha = 0.2)
ggsave( "skew_pLDDT_density.png", height = 3.11, width = 7.92,path = "~/splicing_project/Data/visuals/", device='png', dpi=700)
ks.test(skew$skew, skew1$skew1)
ks.test(skew$skew, skew2$skew2)
ks.test(skew1$skew1, skew2$skew2)
