library(stringr)

setwd("~/splicing_project/")


## get aligned sequence
get_seq = function(file){
        data = readLines(file)
        protein = character(0)
        for (i in unlist(strsplit(data[32], " "))){
                if (nchar(i) >= 3){
                        protein = c(protein, i)
                }
        }
        for (i in unlist(strsplit(data[36], " "))){
                if (nchar(i) >= 3){
                        protein = c(protein, i)
                }
        }
        protein = paste(protein, collapse="")
        print("This is alignment")
        print(protein)
        write(protein, "best_aligned_protein.txt")
}

files = c("alignment_sense_1.water", "alignment_sense_2.water", "alignment_sense_3.water", "alignment_sense_4.water", "alignment_sense_5.water", "alignment_sense_6.water", "alignment_antisense_1.water", "alignment_antisense_2.water", "alignment_antisense_3.water", "alignment_antisense_4.water", "alignment_antisense_5.water", "alignment_antisense_6.water")
scores = numeric(0)
if (file.exists("alignment_sense_1.water")){
        ## Extracting all alignment scores
        
        for (file in files){   
                align = readLines(file)
                for (i in align){
                        if (startsWith(i, "# Score:")){
                                scores = c(scores,as.numeric(unlist(strsplit(i, " "))[3]))
                        }
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
                }
                if (sc==2){
                        write(c("5'-3'", "sense", "frame_2"), "best_align_meta.txt")
                }
                if (sc==3){
                        write(c("5'-3'", "sense", "frame_3"), "best_align_meta.txt")
                }
                if (sc==4){
                        write(c("3'-5'", "sense", "frame_1"), "best_align_meta.txt")
                }
                if (sc==5){
                        write(c("3'-5'", "sense", "frame_2"), "best_align_meta.txt")
                }
                if (sc==6){
                        write(c("3'-5'", "sense", "frame_3"), "best_align_meta.txt")
                }
        } else {
                if (sc-6==1){
                        write(c("5'-3'", "antisense", "frame_1"), "best_align_meta.txt")
                }
                if (sc-6==2){
                        write(c("5'-3'", "antisense", "frame_2"), "best_align_meta.txt")
                }
                if (sc-6==3){
                        write(c("5'-3'", "antisense", "frame_3"), "best_align_meta.txt")
                }
                if (sc-6==4){
                        write(c("3'-5'", "antisense", "frame_1"), "best_align_meta.txt")
                }
                if (sc-6==5){
                        write(c("3'-5'", "antisense", "frame_2"), "best_align_meta.txt")
                }
                if (sc-6==6){
                        write(c("3'-5'", "antisense", "frame_3"), "best_align_meta.txt")
                }
        }
        
} else {
        write("NOT FOUND", "best_aligned_protein.txt")
        write(c("", "", ""), "best_align_meta.txt")
}

