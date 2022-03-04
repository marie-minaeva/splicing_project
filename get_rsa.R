library(bio3d)
library(stringr)



data = read.csv("Data/output_non_sQTL_full.csv")
names = data$ALPHAFOLD.NAME
poss = data$ALIGN.COORDS
poss
mean_acc = numeric(0)
nrow(data)
turn = logical(0)
she = logical(0)
hel = logical(0)
for (i in 1:nrow(data)){
        print(i)
        start = as.numeric(unlist(str_split(poss[i], pattern = '-'))[1])
        start
        stop = as.numeric(unlist(str_split(poss[i], pattern = '-'))[2])
        file = paste("~/splicing_project/UP000005640_9606_HUMAN/AF-", names[i], "-F1-model_v1.pdb", sep='') 
        if (file.exists(file) & !is.na(start) & !is.na(stop)){
                pdb = read.pdb(file)
                out = dssp(pdb, exefile = "mkdssp", resno=TRUE, full=FALSE, verbose=FALSE)
                if (!is.na(start) || !is.na(stop)){
                        mean_acc = c(mean_acc, mean(out$acc[start:stop], na.rm=T))       
                } else {
                        mean_acc = c(mean_acc, NA)
                }
        
                flag = F
                if (!is.null(out$helix$start)){
                for (i in 1:length(out$helix$start)){
                        temp = out$helix$start[i]:out$helix$end[i]
                        if (start:stop %in% temp){
                                flag = T
                        }
                }}
                hel = c(hel, flag)
                
                flag = F
                if (!is.null(out$sheet$start)){
                for (i in 1:length(out$sheet$start)){
                        temp = out$sheet$start[i]:out$sheet$end[i]
                        if (start:stop %in% temp){
                                flag = T
                        }
                }}
                she = c(she, flag)
                
                flag = F
                if (!is.null(out$turn$start)){
                for (i in 1:length(out$turn$start)){
                        temp = out$turn$start[i]:out$turn$end[i]
                        if (start:stop %in% temp){
                                flag = T
                        }
                }}
                turn = c(turn, flag)
        } else {
                mean_acc = c(mean_acc, NA)
                turn = c(turn, F)
                she = c(she, F)
                hel = c(hel, F)
        }
}
she
hel
turn
mean_acc
data$MEAN.ACC = mean_acc 
data$HELIX = hel
data$SHEET = she
data$TURN = turn
#View(data)
write.csv(data, "Data/output_non_sQTL_full.csv", row.names = F)
