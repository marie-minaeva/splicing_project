library(stringr)

setwd("~/splicing_project/")
transeq_parse <- function(input) {
        temp = character(0)
        proteins = character(0)
        for (i in 1:nrow(input)){
                if (startsWith(input[i,], ">")){
                        pars = unlist(strsplit(input[i,],split = ":"))[2]
                        pars = unlist(strsplit(pars,split = "-"))
                        proteins = c(proteins, temp)
                        temp = character(0)
                }else{
                        temp = c(temp, input[i,])
                        temp = paste(temp, collapse="")
                }
                if (i == nrow(input)){
                        proteins = c(proteins, temp) 
                }
        }
        return(proteins)
}


## get aligned sequence
get_seq = function(file){
        data = readLines(file)
        protein = character(0)
        print(length(data))
        coords = numeric(0)
        for (j in seq(32, length(data), by=4)){
                for (i in unlist(strsplit(data[j], " "))){
                        if (nchar(i) >= 3){
                                protein = c(protein, i)
                        }
                }
        }

        
        
        for (j in seq(34, length(data), by=4)){
                coord = gsub("[^0-9.]", " ", data[j])
                coords = c(coords, unlist(strsplit(coord, " ")))
                        }
        coords = coords[coords!= ""]
        coords = as.numeric(coords)
        protein = paste(protein, collapse="")
        if (str_detect(protein, "-")){
                protein = str_remove_all(protein, "-")
                write("EXON NOT FOUND", "best_aligned_protein.txt")
                write(c("", "", ""), "best_align_meta.txt")
                protein = gsub("[0-9]", "", protein)
                write(protein, "protein_seq.txt")
                write('', "align_coords.txt")
                break
        } else {
                protein = gsub("[0-9]", "", protein)
                write(protein, "best_aligned_protein.txt")
                write(paste(c(coords[1], coords[length(coords)]), collapse = "-"), "align_coords.txt")
        }
}

output = read.table("compl_protein_antisense_seq.txt")
proteins_anti = transeq_parse(output)
output = read.table("protein_sense_seq.txt")
proteins = transeq_parse(output)
proteins = proteins[!grepl("[*]", proteins)]
proteins_anti = proteins_anti[!grepl("[*]", proteins_anti)]
files = c("alignment_sense_1.water", "alignment_sense_2.water", "alignment_sense_3.water", "alignment_sense_4.water", "alignment_sense_5.water", "alignment_sense_6.water", "alignment_antisense_1.water", "alignment_antisense_2.water", "alignment_antisense_3.water", "alignment_antisense_4.water", "alignment_antisense_5.water", "alignment_antisense_6.water")
scores = numeric(0)
if (file.exists("alignment_sense_1.water") || file.exists("alignment_sense_2.water") || file.exists("alignment_sense_3.water") || file.exists("alignment_sense_4.water") || file.exists("alignment_sense_5.water") || file.exists("alignment_sense_6.water") || file.exists("alignment_antisense_1.water") || file.exists("alignment_antisense_2.water") || file.exists("alignment_antisense_3.water") || file.exists("alignment_antisense_4.water") || file.exists("alignment_antisense_5.water") || file.exists("alignment_antisense_6.water")){
        ## Extracting all alignment scores
        
        for (file in files){
                if (file.exists(file)){
                        align = readLines(file)
                        for (i in align){
                                if (startsWith(i, "# Score:")){
                                        scores = c(scores,as.numeric(unlist(strsplit(i, " "))[3]))
                                }
                        }
                } else {
                        scores = c(scores, 0.0)
                }
        }

        print(scores)
        ## Getting maximum score
        sc = which.max(scores)
        ## Selecting best aligned protein sequence
        get_seq(files[sc])
        
        ## Selecting metadata for the best aligned protein sequence
        if (sc<=6){
                if (sc==1){
                        write(c("5'-3'", "sense", "frame_1"), "best_align_meta.txt")
                        write(proteins[sc], "protein_seq.txt")
                        
                }
                if (sc==2){
                        write(c("5'-3'", "sense", "frame_2"), "best_align_meta.txt")
                        write(proteins[sc], "protein_seq.txt")
                }
                if (sc==3){
                        write(c("5'-3'", "sense", "frame_3"), "best_align_meta.txt")
                        write(proteins[sc], "protein_seq.txt")
                }
                if (sc==4){
                        write(c("3'-5'", "sense", "frame_1"), "best_align_meta.txt")
                        write(proteins[sc], "protein_seq.txt")
                }
                if (sc==5){
                        write(c("3'-5'", "sense", "frame_2"), "best_align_meta.txt")
                        write(proteins[sc], "protein_seq.txt")
                }
                if (sc==6){
                        write(c("3'-5'", "sense", "frame_3"), "best_align_meta.txt")
                        write(proteins[sc], "protein_seq.txt")
                }
        } else {
                if (sc-6==1){
                        write(c("3'-5'", "antisense", "frame_1"), "best_align_meta.txt")
                        write(proteins_anti[sc-6], "protein_seq.txt")
                }
                if (sc-6==2){
                        write(c("3'-5'", "antisense", "frame_2"), "best_align_meta.txt")
                        write(proteins_anti[sc-6], "protein_seq.txt")
                }
                if (sc-6==3){
                        write(c("3'-5'", "antisense", "frame_3"), "best_align_meta.txt")
                        write(proteins_anti[sc-6], "protein_seq.txt")
                }
                if (sc-6==4){
                        write(c("5'-3'", "antisense", "frame_1"), "best_align_meta.txt")
                        write(proteins_anti[sc-6], "protein_seq.txt")
                }
                if (sc-6==5){
                        write(c("5'-3'", "antisense", "frame_2"), "best_align_meta.txt")
                        write(proteins_anti[sc-6], "protein_seq.txt")
                }
                if (sc-6==6){
                        write(c("5'-3'", "antisense", "frame_3"), "best_align_meta.txt")
                        write(proteins_anti[sc-6], "protein_seq.txt")
                }
        }
        
} else {
        write("NOT IN ALPHAFOLD", "best_aligned_protein.txt")
        write(c("", "", ""), "best_align_meta.txt")
        write("", "protein_seq.txt")
        write('', "align_coords.txt")
}

