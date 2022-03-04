setwd("~/splicing_project/")
library(dplyr)
library(ggplot2)
library(rfUtilities)

flag_outliers <- function(x){
        quants <- quantile(x, na.rm=TRUE)[c(2,4)]
        iqr <- IQR(x, na.rm = T)
        return((x > quants[2] + 1.5 * iqr) | (x < quants[1] - 1.5 * iqr))
}

data = read.csv("Data/output_top_sQTL.csv")
data2 = read.csv("Data/output_non_sQTL_full.csv")

nrow(data)
nrow(data2)
data = data[,2:ncol(data)]
data2 = data2[,2:ncol(data2)]
#View(data)

data = unique.data.frame(data)
data2 = unique.data.frame(data2)
# View(data)
nrow(data)
nrow(data2)



data = data[!flag_outliers(data$LENGTH),]
data2 = data2[!flag_outliers(data2$LENGTH),]
nrow(data)
nrow(data2)
data = sample_n(data, nrow(data2))
data = cbind(data, 1) ## 1 stands for sQTLs
colnames(data)[ncol(data)] = "GROUP"
View(data)
data2 = cbind(data2, 0) ## 0 stands for non_sQTLs
colnames(data2)[ncol(data2)] = "GROUP"
View(data2)
ncol(data)
ncol(data2)
setdiff(colnames(data), colnames(data2))
data = select(data, -setdiff(colnames(data), colnames(data2)))

total_data = rbind(data, data2)
View(total_data)
total_data = total_data[!is.na(total_data$MIN), ]
View(total_data)
total_data = select(total_data, -c(GENE.NAME))
train = sample_n(total_data, 4000)
test = total_data[! total_data$ALPHAFOLD.NAME %in% train$ALPHAFOLD.NAME, ]



classifier = glm(GROUP~MIN+MAX_pLLDT, family='binomial', data=total_data)
exp(coef(classifier))
summary(classifier)
# anova(classifier, test="LRT")
to_plot = sample_n(total_data, 100)
nrow(data)
nrow(data2)
accuracy(round(predict(classifier, total_data, type="response"), 0), total_data$GROUP)
ggplot() + geom_density(data=total_data, aes(MIN, fill=as.factor(GROUP)), alpha = 0.2)

