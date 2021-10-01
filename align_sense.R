## parsing transeq output
transeq_parse <- function(input) {
        temp = character(0)        
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

# Preprocessing for alignment
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