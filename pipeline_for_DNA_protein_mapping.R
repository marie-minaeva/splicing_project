library(seqinr)
seq_1 = "atccgtattgcagtaca"
seq = unlist(strsplit(seq_1,split = '')) 
seq
for (i in 0:2){
        print(translate(seq, frame = i, sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE))
}







## reading data
data = read.table("~/splicing_project/significant_sQTL_GWAS_coloc_events.tsv", header=T)

## main cycle
for (j in 1:2){
gene_name = needed_coords = unlist(strsplit(data$phenotype_id[j],split = "_"))[1]
bash = c("bash ~/splicing_project/retrieve_nucleotide_sequence_from_genename.sh -g", gene_name)
bash =  paste(bash, collapse=" ")
print(bash)
## exon extraction
system(bash)
output = read.table("~/splicing_project/output.txt", sep="\n")

## parsing bedtools output
starts = character(0)
stops = character(0)
temp = character(0)
exons = character(0)
n=length(rownames(output))
for (i in 1:n){
        if (startsWith(output[i,], ">chr")){
                pars = unlist(strsplit(output[i,],split = ":"))[2]
                pars = unlist(strsplit(pars,split = "-"))
                exons = c(exons, temp)
                starts = c(starts, pars[1])
                stops = c(stops, pars[2])
                temp = character(0)
        }else{
                temp = c(temp, output[i,])
                temp = paste(temp, collapse="")
        }
        if (i == n){
                exons = c(exons, temp) 
        }
}
print(data$Exon_coord[j])
needed_coords = unlist(strsplit(data$Exon_coord[j],split = "_"))[2:3]
print(needed_coords)
print(starts)
print(stops)
print(match(needed_coords[1], starts))
print(match(needed_coords[2], stops))

}

