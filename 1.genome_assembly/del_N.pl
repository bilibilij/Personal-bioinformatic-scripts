#!/usr/bin/envs perl -w

use strict;

my $usage=<<USAGE;

perl $0 fasta > fasta.delfirstN

Some fasta have the N x 500 at the beginning, this script is to delete first (N x 500) which fasta has

USAGE

die "$usage \n " unless @ARGV == 1;

open IN, "$ARGV[0]" or die "";
my ($seqID, %seq);
while (<IN>){
	chomp;
	if (m/^>(\S+)/){$seqID=$1;}
	else{$seq{$seqID}.=$_;}
}

close IN;



foreach my $id (sort keys %seq){
	my $seq=$seq{$id};
	my $N=("N") x 500;
	if ($seq =~ m/$N/ ){ 
	        my $N=("N") x 500;
       	        $seq=~ s/$N//;
     		print ">$id\n$seq\n";
	}else{print ">$id\n$seq\n";}
}

