setwd("~/splicing_project/")
dev.off()
library(ggplot2)

gtex_v8_figure_theme <- function() {
        return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


data1 = read.csv("Data/median_psi/Brain_Cerebellum_median_psi.tsv", header = T, sep = '\t')
data2 = read.csv("Data/median_psi/Brain_Cortex_median_psi.tsv", header = T, sep = '\t')
data3 = read.csv("Data/median_psi/Brain_Nucleus_accumbens_basal_ganglia_median_psi.tsv", header = T, sep = '\t')



data4 = read.csv("Data/median_psi/Adipose_Subcutaneous_median_psi.tsv", header = T, sep = '\t')
data5 = read.csv("Data/median_psi/Artery_Tibial_median_psi.tsv", header = T, sep = '\t')
data6 = read.csv("Data/median_psi/Cells_Cultured_fibroblasts_median_psi.tsv", header = T, sep = '\t')


data7 = read.csv("Data/median_psi/Cells_EBV.transformed_lymphocytes_median_psi.tsv", header = T, sep = '\t')
data8 = read.csv("Data/median_psi/Colon_Transverse_median_psi.tsv", header = T, sep = '\t')
data9 = read.csv("Data/median_psi/Esophagus_Mucosa_median_psi.tsv", header = T, sep = '\t')

data10 = read.csv("Data/median_psi/Liver_median_psi.tsv", header = T, sep = '\t')
data11 = read.csv("Data/median_psi/Lung_median_psi.tsv", header = T, sep = '\t')
data12 = read.csv("Data/median_psi/Muscle_Skeletal_median_psi.tsv", header = T, sep = '\t')


data13 = read.csv("Data/median_psi/Nerve_Tibial_median_psi.tsv", header = T, sep = '\t')
data14 = read.csv("Data/median_psi/Pituitary_median_psi.tsv", header = T, sep = '\t')
data15 = read.csv("Data/median_psi/Skin_Sun_Exposed_Lower_leg_median_psi.tsv", header = T, sep = '\t')


data16 = read.csv("Data/median_psi/Spleen_median_psi.tsv", header = T, sep = '\t')
data17 = read.csv("Data/median_psi/Thyroid_median_psi.tsv", header = T, sep = '\t')
data18 = read.csv("Data/median_psi/Whole_Blood_median_psi.tsv", header = T, sep = '\t')



ggplot() + gtex_v8_figure_theme() + geom_density(data = data1, aes(mean_psi, fill = "Brain_Cerebellum"), alpha = 0.2) + 
        geom_density(data = data2, aes(mean_psi, fill = "Brain_Cortex"), alpha = 0.2) + 
        geom_density(data = data3, aes(mean_psi, fill = "Brain_Nucleus_accumbens_basal_ganglia"), alpha = 0.2) + 
        geom_density(data = data4, aes(mean_psi, fill = "Adipose_Subcutaneous"), alpha = 0.2) + 
        geom_density(data = data5, aes(mean_psi, fill = "Artery_Tibial"), alpha = 0.2) + 
        geom_density(data = data6, aes(mean_psi, fill = "Cells_Cultured_fibroblasts"), alpha = 0.2) + 
        geom_density(data = data7, aes(mean_psi, fill = "Cells_EBV.transformed_lymphocytes"), alpha = 0.2) + 
        geom_density(data = data8, aes(mean_psi, fill = "Colon_Transverse"), alpha = 0.2) + 
        geom_density(data = data9, aes(mean_psi, fill = "Esophagus_Mucosa"), alpha = 0.2) + 
        geom_density(data = data10, aes(mean_psi, fill = "Liver"), alpha = 0.2) + 
        geom_density(data = data11, aes(mean_psi, fill = "Lung"), alpha = 0.2) + 
        geom_density(data = data12, aes(mean_psi, fill = "Muscle_Skeletal"), alpha = 0.2) + 
        geom_density(data = data13, aes(mean_psi, fill = "Nerve_Tibial"), alpha = 0.2) + 
        geom_density(data = data14, aes(mean_psi, fill = "Pituitary"), alpha = 0.2) + 
        geom_density(data = data15, aes(mean_psi, fill = "Skin_Sun_Exposed_Lower_leg"), alpha = 0.2) + 
        geom_density(data = data16, aes(mean_psi, fill = "Spleen"), alpha = 0.2) + 
        geom_density(data = data17, aes(mean_psi, fill = "Thyroid"), alpha = 0.2) + 
        geom_density(data = data18, aes(mean_psi, fill = "Whole_Blood"), alpha = 0.2) + 
        ggtitle("mean_psi")

ggsave("inclusion_level_across_tissues.png", path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)
