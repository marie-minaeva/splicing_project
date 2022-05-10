library(rjson)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)


plot_pae = function(data, data_2, pos1, pos2){
        mat1 = as.data.frame(cbind(data[[1]]$residue1, data[[1]]$residue2, data[[1]]$distance))
        mat1 = matrix(mat1$V3, byrow = T, nrow = max(mat1$V1))
        print(max(data_2[[1]]$residue1))
        new_mat1 = matrix(NA, nrow = max(data_2[[1]]$residue1), max(data_2[[1]]$residue1))
        new_mat1[-c(pos1:pos2), -c(pos1:pos2)] = mat1
        to_plot = melt(new_mat1)
        # to_plot = as.data.frame(cbind(data[[1]]$residue1, data[[1]]$residue2, data[[1]]$distance))
        p_heatmap_1 <- ggplot(to_plot, aes(x=Var1, y=Var2)) +
                geom_tile(aes(fill=value)) +
                scale_fill_gradient2("pae", midpoint = 0, low = "#9161A8", high = "#F7941E", mid = "white", na.value = "grey90") +
                theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text.x = element_text(angle = 90)) +
                theme(legend.position = "bottom", legend.direction = "horizontal") +
                guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) + ggtitle("Mutant")
        p_heatmap_1
        
        to_plot = as.data.frame(cbind(data_2[[1]]$residue1, data_2[[1]]$residue2, data_2[[1]]$distance))
        p_heatmap_2 <- ggplot(to_plot, aes(x=V1, y=V2)) +
                geom_tile(aes(fill=V3)) +
                scale_fill_gradient2("pae", midpoint = 0, low = "#9161A8", high = "#F7941E", mid = "white", na.value = "grey90") +
                theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text.x = element_text(angle = 90)) +
                theme(legend.position = "bottom", legend.direction = "horizontal") +
                guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) + ggtitle("WT")
        p_heatmap_2
        
        
        
        lay <- rbind(c(1,2),
                     c(1,2))
        grid.arrange(grobs = list(p_heatmap_1, p_heatmap_2), layout_matrix = lay)
}


plot_pae_corr = function(data, data_2, pos1, pos2){
        mat1 = as.data.frame(cbind(data[[1]]$residue1, data[[1]]$residue2, data[[1]]$distance))
        mat2 = as.data.frame(cbind(data_2[[1]]$residue1, data_2[[1]]$residue2, data_2[[1]]$distance))
        mat1 = matrix(mat1$V3, byrow = T, nrow = max(mat1$V1))
        mat2 = matrix(mat2$V3, byrow = T, nrow = max(mat2$V1))
        new_mat1 = matrix(NA, nrow = nrow(mat2), ncol=ncol(mat2))
        new_mat1[-c(pos1:pos2), -c(pos1:pos2)] = mat1
        cor_mat = cor(new_mat1, mat2, use="pairwise.complete.obs")
        View(cor_mat)
        melted_cor_mat = melt(cor_mat)
        p_heatmap <- ggplot(melted_cor_mat, aes(x=Var1, y=Var2)) +
                geom_tile(aes(fill=value)) +
                scale_fill_gradient2("Pearson\ncorrelation", midpoint = 0, low = "#9161A8", high = "#F7941E", mid = "white", na.value = "grey90") +
                theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text.x = element_text(angle = 90)) +
                theme(legend.position = "bottom", legend.direction = "horizontal") +
                guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5))
        p_heatmap
}


setwd("~/splicing_project/")
## A6XGL0 MUTANT
data <- fromJSON(file = "Data/visuals/AlphaFold_predictions/A6XGL0_no_182_231_d4286.result/A6XGL0_no_182_231_d4286_predicted_aligned_error_v1.json")
data_2 = fromJSON(file="Data/visuals/AlphaFold_predictions/AF-A6XGL0-F1-predicted_aligned_error_v2.json")
plot = plot_pae(data, data_2, 182, 231)
ggsave(plot = plot, filename = "A6XGL0_no_182_231_PAE.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)
plot = plot_pae_corr(data, data_2, 182, 231)
ggsave(plot = plot, filename = "A6XGL0_no_182_231_PAE_corr.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)

## O60641 MUTANT
data <- fromJSON(file = "Data/visuals/AlphaFold_predictions/O60641_no_567_591_c0aef.result/O60641_no_567_591_c0aef_predicted_aligned_error_v1.json")
data_2 = fromJSON(file="Data/visuals/AlphaFold_predictions/AF-O60641-F1-predicted_aligned_error_v2.json")
plot = plot_pae(data, data_2, 567, 591)
ggsave(plot = plot, filename = "O60641_no_567_591_PAE.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)
plot = plot_pae_corr(data, data_2, 567, 591)
ggsave(plot = plot, filename = "O60641_no_567_591_PAE_corr.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)

## P04035 MUTANT
data <- fromJSON(file = "Data/visuals/AlphaFold_predictions/P04035_no_522574_d3855.result/P04035_no_522574_d3855_predicted_aligned_error_v1.json")
data_2 = fromJSON(file="Data/visuals/AlphaFold_predictions/AF-P04035-F1-predicted_aligned_error_v2.json")
plot = plot_pae(data, data_2, 522, 574)
ggsave(plot = plot, filename = "P04035_no_522_574_PAE.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)
plot = plot_pae_corr(data, data_2, 522, 574)
ggsave(plot = plot, filename = "P04035_no_522_574_PAE_corr.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


## Q00722 MUTANT
data <- fromJSON(file = "Data/visuals/AlphaFold_predictions/Q00722_no_1119_1185_a9c16.result/Q00722_no_1119_1185_a9c16_predicted_aligned_error_v1.json")
data_2 = fromJSON(file="Data/visuals/AlphaFold_predictions/AF-Q00722-F1-predicted_aligned_error_v2.json")
plot = plot_pae(data, data_2, 1119, 1185)
ggsave(plot = plot, filename = "Q00722_no_1119_1185_PAE.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)
plot = plot_pae_corr(data, data_2, 1119, 1185)
ggsave(plot = plot, filename = "Q00722_no_1119_1185_PAE_corr.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


## Q13342 MUTANT
data <- fromJSON(file = "Data/visuals/AlphaFold_predictions/Q13342_no_223_247_41432.result/Q13342_no_223_247_41432_predicted_aligned_error_v1.json")
data_2 = fromJSON(file="Data/visuals/AlphaFold_predictions/AF-Q13342-F1-predicted_aligned_error_v2.json")
plot = plot_pae(data, data_2, 223, 247)
ggsave(plot = plot, filename = "Q13342_no_223_247_PAE.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)
plot = plot_pae_corr(data, data_2, 223, 247)
ggsave(plot = plot, filename = "Q13342_no_223_247_PAE_corr.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)


## Q92985 MUTANT
data <- fromJSON(file = "Data/visuals/AlphaFold_predictions/Q92985_no_152_226_dad09.result/Q92985_no_152_226_dad09_predicted_aligned_error_v1.json")
data_2 = fromJSON(file="Data/visuals/AlphaFold_predictions/AF-Q92985-F1-predicted_aligned_error_v2.json")
plot = plot_pae(data, data_2, 152, 226)
ggsave(plot = plot, filename = "Q92985_no_152_226_PAE.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)
plot = plot_pae_corr(data, data_2, 152, 226)
ggsave(plot = plot, filename = "Q92985_no_152_226_PAE_corr.png", path = "Data/visuals/AlphaFold_predictions/", height = 5.11, width = 7.92,device='png', dpi=700)

