#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;


sub TranslateDNAFile(){#undefined
(my $infile,my $outfile)=@_;
	my $in=Bio::SeqIO->new(-file=>"$infile",-format=>"fasta");
	my $out=Bio::SeqIO->new(-file=>">$outfile",-format=>"fasta");
	while (my $seq=$in->next_seq()){
	$out->write_seq($seq->translate);
}
}
my $DNAfile=$ARGV[0];
my $pepfile=$ARGV[1];
&TranslateDNAFile($DNAfile,$pepfile);
