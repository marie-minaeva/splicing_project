#! /usr/bin/perl 
#use strict;

# this script is written by Ulrike Löber (contact uloeber@uni-potsdam.de)
# function: parsing "multiple" EMBOSS pairwise (Smith-Waterman) sequence alignments
# output is a tabular summary of the results containing information of identity, start and end points,...
# this is missing in default outputs, complete information not contained in tabular output by EMBOSS

#start with read in of the data

#open (INFILE, "$ARGV[0]") or die "no such file $ARGV[0]\n";

#initialize variables
my $bool=0;
my $check_odd_even;
my @splitstring;
my $id1;
my $id2;
my $length;
my $identity;
my $similarity;
my $gaps;
my $score;
my $start1;
my $end1;
my $seq1;
my $start2;
my $end2;
my $seq2;
my $pattern_id1;
my $pattern_id2;

# initialize outfile header 
my $outfile="$ARGV[0].tab";
open (OUTFILE, ">$outfile") or die "error creating file $outfile !\n";
print OUTFILE "id1\tid2\tlength\tidentity\tsimilarity\tgaps\tscore\tstart1\tseq1\tend1\tstart2\tseq2\tend2\n";
 while (defined(my $line=<INFILE>)) {
	chomp $line;
	$line=~s/\R//g;	#remove linebreaks
# first information to extract are the "method" information
	if ($line	=~	m/\#{10,}/g){	#there should be only two lines with at least 10 # in it, before and after the description at the beginning of the file
		$bool	=	$bool+1;
		$check_odd_even=$bool%2;
	}
	if	($check_odd_even  	==	1){
#		print "$line\n";
	}
# following the alignment information, written in the new outfile
	elsif	($check_odd_even	==	0){
		if($line	=~	m/\#\ 1\:/g){			#if line structure is like "# 1: seqname" print out third element of the line 
			@splitstring	=	split(/ +/,$line);	#split sting by whitespaces
			$id1		=	$splitstring[2];	
			print OUTFILE "$id1\t";
			$pattern_id1=substr($id1,0,13);			# attention! alignments: seqids are croped to 13 signs! 
		}	
		elsif($line	=~	m/\#\ 2\:/g){			#if line structure is like "# 2: seqname" print out third element of the line
			@splitstring	=	split(/ +/,$line);	
			$id2		=	$splitstring[2];
			print OUTFILE "$id2\t";
			$pattern_id2=substr($id2,0,13);			# attention! alignments: seqids are croped to 13 signs! 
		}
		elsif($line=~m/\#\ Length\:/ig){			#if line structure is like "# Length: 24" print out third element of the line
			@splitstring	=	split(/ +/,$line);
			$length		=	$splitstring[2];
			print OUTFILE "$length\t";
		}
		elsif($line=~m/\#\ Identity\:/ig){			#if line structure is like "# Identity: 9/10 (90.0%)" 
			@splitstring=split(/ +/,$line);			#split sting by whitespaces
			$identity	=	$splitstring[3];	#fourth element is the identity in percentage
			$identity	=~	s/[\(,\),\%]//g;	#cut off the braces and %-sign
			print OUTFILE "$identity\t";
		}
		elsif($line=~m/\#\ Similarity\:/ig){
			@splitstring=split(/ +/,$line);
			$similarity	=	$splitstring[3];
			$similarity	=~	s/[\(,\),\%]//g;
			print OUTFILE "$similarity\t";
		}
		elsif($line=~m/\#\ Gaps\:/ig){
			@splitstring=split(/ +/,$line);
			$gaps	=	$splitstring[3];
			$gaps	=~	s/[\(,\),\%]//g;
			print OUTFILE "$gaps\t";
		}
		elsif($line=~m/\#\ Score\:/ig){
			@splitstring	=	split(/ +/,$line);
			$score		=	$splitstring[2];
			print OUTFILE "$score\t";
		}
		
		
		elsif($line=~m/^$pattern_id1\ /g){
			@splitstring	=	split(/ +/,$line);
			$start1		=	$splitstring[1];
			$seq1		=	$splitstring[2];
			$end1		=	$splitstring[3];
			print OUTFILE "$start1\t$seq1\t$end1\t";
		}
		
		
		elsif($line=~m/^$pattern_id2\ /g){
			@splitstring	=	split(/ +/,$line);
			$start2		=	$splitstring[1];
			$seq2		=	$splitstring[2];
			$end2		=	$splitstring[3];
			print OUTFILE "$start2\t$seq2\t$end2\t";
		}
		elsif($line=~m/Aligned_sequences/g){			#insert linebreak, when next alignment starts
			print OUTFILE "\n";
		}
	}
}
close OUTFILE;
print "\n Now parsed to tabular format and written in $outfile\n ";

