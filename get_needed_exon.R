args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])
data = read.table("non_colocalizing_sQTLs.tsv", header=T)

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

## bedtools parsing
output = read.table("output.txt", sep="\n")
res = bedtools_parse(output)
starts = unlist(res[1])
stops = unlist(res[2])
exons = unlist(res[3])
pos = as.numeric(unlist(strsplit(data$phenotype_id[j],split = "_"))[2])
pos
## needed exon extraction
#print(data$start[j])
#needed_coords = unlist(strsplit(data$Exon_coord[j],split = "_"))[2:3]
print(pos)
print(starts[pos])
print(stops[pos])
#rm(output)
#write(match(needed_coords[2], stops), "exon_meta.txt")
write(pos, "exon_meta.txt")
#needed_exon = exons[match(needed_coords[2], stops)]
needed_exon = exons[pos]
print(needed_exon)
write(needed_exon, "exon_seq.fa")