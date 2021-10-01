library(seqinr)
seq_1 = "atccgtattgcagtaca"
seq = unlist(strsplit(seq_1,split = '')) 
seq
for (i in 0:2){
        print(translate(seq, frame = i, sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE))
}




library(spgs)

## parsing transeq output
transeq_parse <- function(input) {
        
        proteins = character(0)
        for (i in 1:nrow(input)){
                if (startsWith(input[i,], ">_")){
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

## parsing bedtools output
bedtools_parse = function(input){
        sta = character(0)
        sto = character(0)
        temp = character(0)
        exo = character(0)
        for (i in 1:nrow(input)){
                if (startsWith(input[i,], ">chr")){
                        pars = unlist(strsplit(input[i,],split = ":"))[2]
                        pars = unlist(strsplit(pars,split = "-"))
                        exo = c(exo, temp)
                        sta = c(sta, pars[1])
                        sto = c(sto, pars[2])
                        temp = character(0)
                }else{
                        temp = c(temp, input[i,])
                        temp = paste(temp, collapse="")
                }
                if (i == nrow(input)){
                        exo = c(exo, temp) 
                }
        }
        return(list(sta, sto, exo))
}







## Setting directory
setwd("~/splicing_project/")

## reading data
data = read.table("significant_sQTL_GWAS_coloc_events.tsv", header=T)

## main cycle
for (j in 1:1){
        gene_name = needed_coords = unlist(strsplit(data$phenotype_id[j],split = "_"))[1]
        bash = c("bash retrieve_nucleotide_sequence_from_genename.sh -g", gene_name)
        bash =  paste(bash, collapse=" ")
        print(bash)
        ## all gene exons extraction
        system(bash)
        #file.show("~/splicing_project/output.txt")
        
        ## bedtools parsing
        output = read.table("output.txt", sep="\n")
        res = bedtools_parse(output)
        starts = unlist(res[1])
        stops = unlist(res[2])
        exons = unlist(res[3])
        
        ## needed exon extraction
        print(data$Exon_coord[j])
        needed_coords = unlist(strsplit(data$Exon_coord[j],split = "_"))[2:3]
        print(needed_coords)
        print(stops)
        print(match(needed_coords[2], stops))
        #rm(output)

        needed_exon = exons[match(needed_coords[2], stops)]
        print(needed_exon)
        write(needed_exon, "exon_seq.fa")
        
        ## Translation
        system("transeq -frame 6 exon_seq.fa protein_sense_seq.txt")
        system("cat protein_sense_seq.txt")
        
        ## finding complement sequence
        comp_needed_exon = complement(needed_exon)
        
        ## complement translation
        write(comp_needed_exon, "compl_exon_seq.fa")
        system("transeq -frame 6 compl_exon_seq.fa compl_protein_sense_seq.txt")
        system("cat compl_protein_sense_seq.txt")
        
        ## parsing ALphaFold database
        
        trial_ref_seq = "PIGRLETPVVGDVLPQGVPPVHRLPVHAVVAVLLYHALGLTLEGLHG*VLPPWPEVPVLVILPS"
        write(trial_ref_seq, "trial_prot_ref.fa")
        
        ## Preprocessing for alignment
        output = read.table("protein_sense_seq.txt")
        proteins = transeq_parse(output)
        
        ## Aligning sense proteins
        for (k in 1:length(proteins)){
                print("Sense")
                print(proteins[k])
                write(proteins[k], "temp_prot.fa")
                out_path = paste(c("alignment_sense_", as.character(k),".water"), collapse="")
                water = c("water -gapopen 10.0 -gapextend 0.5", "temp_prot.fa", "trial_prot_ref.fa", out_path)
                water =  paste(water, collapse=" ")
                system(water)
                
        }
        
        ## Preprocessing for alignment
        output = read.table("compl_protein_sense_seq.txt")
        proteins = transeq_parse(output)
        print(proteins)
        ## Aligning antisense proteins
        for (k in 1:length(proteins)){
                print("Anti")
                print(proteins[k])
                write(proteins[k], "temp_prot.fa")
                out_path = paste(c("alignment_antisense_", as.character(k),".water"), collapse="")
                water = c("water -gapopen 10.0 -gapextend 0.5", "temp_prot.fa", "trial_prot_ref.fa", out_path)
                water =  paste(water, collapse=" ")
                system(water)
                
        }
        
}

