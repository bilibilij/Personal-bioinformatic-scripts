#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;

perl $0 all.fasta idlist > id.fasta

USAGE

die $usage unless @ARGV == 2;


open IN, "$ARGV[0]" or die "";

my ($seqID, %seq);
while (<IN>){
	chomp;
	if (m/^>(\S+)/){
		$seqID=$1;
	}else{
		$seq{$seqID}.=$_;
	}
}

close IN;

open ID, "$ARGV[1]" or die "";

while (<ID>){
	chomp;
	print ">$_\n$seq{$_}\n";
}

close ID;




