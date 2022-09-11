#!/usr/bin/evns perl -w

use strict;


open TS, "$ARGV[0]" or die "";
my %pos;
while (<TS>){
	chomp;
	next unless m/PF00931/;
	my @line=split/\t/, $_;
	$pos{$line[0]}="$line[6]\t$line[7]";
}
close TS;

my ($seqID, %seq);
open FA, "$ARGV[1]" or die "";
while (<FA>){
        chomp;
        s/\.CNls//g;
        s/\.NLs//g;
        s/\.RNLs//g;
        s/\.TNLs//g;
	if (m/^>(\S+)/){$seqID=$1;}
	else{$seq{$seqID}.=$_;}
}
close FA;



open FA, "$ARGV[1]" or die "";
while (<FA>){
	chomp;
	next unless m/^>/;
	s/>//g;
	my $id=$_;
	s/\.CNls//g;
	s/\.NLs//g;
	s/\.RNLs//g;
	s/\.TNLs//g;
	my $seq=$seq{$_};
	my ($start, $end) = split /\t/,$pos{$_};
	$start=$start-1;
	my $len=$end-$start;
	my $d_seq=substr($seq, $start, $len);
	print ">$id\n$d_seq\n";
}
close FA;





