#!/bin/bash


for i in {1..3790}
        do
        Rscript get_pdb.R $i
	file=$(<name.txt)
	freesasa --format=rsa $file | grep RES > try_freesas.txt
        Rscript try_freeSASA.R $i
	rm name.txt
        rm try_freesas.txt
        done


