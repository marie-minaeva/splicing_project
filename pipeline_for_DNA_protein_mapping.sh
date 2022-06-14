#!/bin/bash


module load R_packages/4.1.1
module load bioinfo-tools
module load blast

for i in {1..5232}
	do
	echo This is round $i
	# Exon extraction
	echo Exon extraction
	Rscript get_needed_exons.R $i
	Rscript get_needed_exon.R $i

	# blastx
	echo blastx
	blastx -query exon_seq.fa -db MANE_prot -out try_blastx.txt -outfmt "6 qseqid sacc sstart send sseq length evalue score"

	# Combine final table
	echo Combine final table
	Rscript get_final_table.R

	rm try_blastx.txt
	rm exon_seq.fa
	done

