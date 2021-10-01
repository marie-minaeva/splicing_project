#!/bin/bash

# Exon extraction
Rscript get_needed_exons.R 1
Rscript get_needed_exon.R 1

# Translation
transeq -frame 6 exon_seq.fa protein_sense_seq.txt

# Complement sequence
revseq exon_seq.fa -noreverse -complement -outseq compl_protein_sense_seq.fa

# Translation of the ccomplement sequence
transeq -frame 6 compl_exon_seq.fa compl_protein_sense_seq.txt

# Parsing AlphaFold
echo PIGRLETPVVGDVLPQGVPPVHRLPVHAVVAVLLYHALGLTLEGLHG*VLPPWPEVPVLVILPS > trial_prot_ref.fa

# Align sense
Rscript align_sense.R

# Aliggn antisense
Rscript align_antisense.R
