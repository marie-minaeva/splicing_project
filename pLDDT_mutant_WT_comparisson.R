library(rjson)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

setwd("~/splicing_project/")

gtex_v8_figure_theme <- function() {
        return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


plot_pLDDT = function(data_mut, data_wt, pos1, pos2){
        data_mut$plddt = append(data_mut$plddt, rep(NA, length(pos1:pos2)), after = pos1)
        to_plot = cbind(1:length(data_wt$plddt), data_wt$plddt, "WT")
        to_plot = rbind(to_plot, cbind(1:length(data_wt$plddt), data_mut$plddt, "Mutant"))
        to_plot = as.data.frame((to_plot))
        colnames(to_plot) = c("position", "pLDDT", "group")
        to_plot$position = as.integer(to_plot$position)
        to_plot$pLDDT = as.numeric(to_plot$pLDDT)
        
        plot = ggplot(to_plot, aes(x=position, y=pLDDT, colour=group)) + 
                geom_line() + gtex_v8_figure_theme()
        return(plot)
}


data_mut = fromJSON(file = "Data/visuals/AlphaFold_predictions/A6XGL0_no_182_231_d4286.result/A6XGL0_no_182_231_d4286_unrelaxed_rank_1_model_4_scores.json")
data_wt = fromJSON(file = "Data/visuals/AlphaFold_predictions/A6XGL0_0a69f.result/A6XGL0_0a69f_unrelaxed_rank_1_model_3_scores.json")
plot = plot_pLDDT(data_mut, data_wt, 182, 231)
plot


data_mut = fromJSON(file = "Data/visuals/AlphaFold_predictions/O60641_no_567_591_c0aef.result/O60641_no_567_591_c0aef_unrelaxed_rank_1_model_3_scores.json")
data_wt = fromJSON(file = "Data/visuals/AlphaFold_predictions/A6XGL0_0a69f.result/A6XGL0_0a69f_unrelaxed_rank_1_model_3_scores.json")
plot = plot_pLDDT(data_mut, data_wt, 182, 231)
plot