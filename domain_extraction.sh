#!/bin/bash


for i in {1..941}
        do
	Rscript get_seq.R $i
	curl -L -H 'Expect:' -H 'Accept:text/xml' -F hmmdb=pfam -F seq='<temp.fa' https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan | grep alihmmdesc > domain.txt
	Rscript get_domain.R $i
	rm temp.fa
	rm domain.txt
	#cat domain.txt
	done
