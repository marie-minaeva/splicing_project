library(ggplot2)

dev.off()
setwd("~/splicing_project/")
data = read.table("Data/top_sQTLs_top_coloc.tsv", header = T, sep='\t')
View(data)
ggplot(data, aes(x=mean_01_psi, fill=has_coloc)) + geom_density(alpha=0.2)

ggsave("inclusion_level_sQTLs.png", path = "Data/visuals/", height = 5.11, width =12.0,device='png', dpi=700)
