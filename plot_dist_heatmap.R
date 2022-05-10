library(bio3d)
library(rjson)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
setwd("~/splicing_project/")
gtex_v8_figure_theme <- function() {
        return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

plot_distaces = function(data, data_2, pos1, pos2){
        mat1=matrix((data$xyz), ncol=3, byrow = T)
        mat2 = matrix((data_2$xyz), ncol=3, byrow = T)
        mat1 = mat1[data$calpha, ]
        mat2 = mat2[data_2$calpha, ]
        new_mat1 = matrix(NA, nrow = nrow(mat2), ncol=3)
        print(dim(mat1))
        print(dim(mat2))
        print(dim(new_mat1))
        # new_mat1[-c(pos1:pos2), ] = mat1
        new_mat1 = mat1
        cor_mat = mat2 - new_mat1
        View(cor_mat)
        cor_mat = as.data.frame(cor_mat)
        colnames(cor_mat) = c("X", "Y", "Z")
        View(cor_mat)
        to_draw = data.frame(cbind(1:length(cor_mat$X), sqrt(cor_mat$X**2 + cor_mat$Y**2 + cor_mat$Z**2)))
        View(to_draw)
        colnames(to_draw) = c("position", "distance")
        p = ggplot(to_draw, aes(x=position, y=distance)) + geom_line() + gtex_v8_figure_theme() + ylab("distance, A")
}



prot2 = read.pdb("Data/visuals/AlphaFold_predictions/A6XGL0_0a69f.result/A6XGL0_0a69f_unrelaxed_rank_1_model_3.pdb")
prot1 = read.pdb("Data/visuals/AlphaFold_predictions/A6XGL0_no_182_231_d4286.result/A6XGL0_no_182_231_d4286_unrelaxed_rank_1_model_4.pdb")
plot = plot_distaces(prot1, prot2, 182, 231)
plot
ggsave(plot = plot, filename = "A6XGL0_distance.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


## A6XGL0 MUTANT
prot1 = read.pdb("Data/visuals/AlphaFold_predictions/A6XGL0_0a69f.result/A6XGL0_0a69f_unrelaxed_rank_1_model_3.pdb")
prot2 = read.pdb("Data/visuals/AlphaFold_predictions/A6XGL0_no_182_231_d4286.result/A6XGL0_no_182_231_d4286_unrelaxed_rank_1_model_4.pdb")
plot = plot_distaces(prot1, prot2, 182, 231)
plot
ggsave(plot = plot, filename = "A6XGL0_distance.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)

## O60641 MUTANT
prot1 = read.pdb("Data/visuals/AlphaFold_predictions/O60641_no_567_591_c0aef.result/O60641_no_567_591_c0aef_unrelaxed_rank_1_model_3.pdb")
prot2 = read.pdb("UP000005640_9606_HUMAN/AF-O60641-F1-model_v1.pdb")
plot = plot_distaces(prot1, prot2, 567, 591)
plot
ggsave(plot = plot, filename = "O60641_distance.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)

## P04035 MUTANT
prot1 = read.pdb("Data/visuals/AlphaFold_predictions/P04035_no_522574_d3855.result/P04035_no_522574_d3855_unrelaxed_rank_1_model_1.pdb")
prot2 = read.pdb("UP000005640_9606_HUMAN/AF-P04035-F1-model_v1.pdb")
plot = plot_distaces(prot1, prot2, 522, 574)
plot
ggsave(plot = plot, filename = "P04035_distance.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


## Q00722 MUTANT
prot1 = read.pdb("Data/visuals/AlphaFold_predictions/Q00722_no_1119_1185_a9c16.result/Q00722_no_1119_1185_a9c16_unrelaxed_rank_1_model_3.pdb")
prot2 = read.pdb("UP000005640_9606_HUMAN/AF-Q00722-F1-model_v1.pdb")
plot = plot_distaces(prot1, prot2, 1119, 1185)
plot
ggsave(plot = plot, filename = "Q00722_distance.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


## Q13342 MUTANT
prot1 = read.pdb("Data/visuals/AlphaFold_predictions/Q13342_no_223_247_41432.result/Q13342_no_223_247_41432_unrelaxed_rank_1_model_1.pdb")
prot2 = read.pdb("UP000005640_9606_HUMAN/AF-Q13342-F1-model_v1.pdb")
plot = plot_distaces(prot1, prot2, 223, 247)
plot
ggsave(plot = plot, filename = "Q13342_distance.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


## Q92985 MUTANT
prot1 = read.pdb("Data/visuals/AlphaFold_predictions/Q92985_no_152_226_dad09.result/Q92985_no_152_226_dad09_unrelaxed_rank_1_model_1.pdb")
prot2 = read.pdb("UP000005640_9606_HUMAN/AF-Q92985-F1-model_v1.pdb")
plot = plot_distaces(prot1, prot2, 152, 226)
plot
ggsave(plot = plot, filename = "Q92985_distance.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)




prot1 = read.pdb("Data/visuals/AlphaFold_predictions/A6XGL0_0a69f.result/A6XGL0_0a69f_unrelaxed_rank_1_model_3.pdb")
prot2 = read.pdb("UP000005640_9606_HUMAN/AF-A6XGL0-F1-model_v1.pdb")
plot = plot_distaces(prot1, prot2, NA, NA)
plot


