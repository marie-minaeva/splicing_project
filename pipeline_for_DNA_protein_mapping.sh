#!/bin/bash


# pretraining biomaRt
Rscript get_biomart_pretrained.R

for i in {1..1}
	do
	echo This is round $i
	# Exon extraction
	echo Exon extraction
	Rscript get_needed_exons.R $i
	Rscript get_needed_exon.R $i

	# Translation
	echo Translation
	transeq -frame 6 exon_seq.fa protein_sense_seq.txt

	# Complement sequence
	echo Complement sequence
	revseq exon_seq.fa -noreverse -complement -outseq compl_exon_seq.fa

	# Translation of the ccomplement sequence
	echo Translation of the ccomplement sequence
	transeq -frame 6 compl_exon_seq.fa compl_protein_antisense_seq.txt

	# Parsing AlphaFold
	echo Parsing AlphaFold
	Rscript get_uniprot_id_from_ensembl.R $i

	# Align sense
	echo Align sense
	Rscript align_sense.R

	# Align antisense
	echo Align antisense
	Rscript align_antisense.R

	# Extracting best protein alignment
	echo Extracting best protein alignment
	Rscript get_best_score.R

	# Extract protein properties
	echo Extract protein properties
	Rscript get_protein_properties.R

	# Combine final table
	echo Combine final table
	Rscript get_final_table.R


	rm alignment_*
	rm align_coords*
	rm best*
	rm compl*
	rm exon*
	rm gene_*
	rm output.*
	rm protein*
	rm ref*
	rm temp*
	done

#python3.9 try_pymol.py
