#!/usr/bin/perl -w 
use strict;
use List::MoreUtils ':all';

my $usage=<<USAGE;
my usage:
	perl $0 interpro.tsv >extracted.GO;
USAGE


open IPR, "$ARGV[0]" or die "can not open interproscan input file";
my %hash;
while (<IPR>) {
	chomp;
	if (m/GO:/){
		my @line = split /\t/,$_;
		$hash{$line[0]}.="$line[13]\|";
	}
}
close IPR;

foreach my $id (keys %hash){
	my $GO= $hash{$id};
	my @GO= split /\|/, $GO;
	my @GOF=uniq(@GO);
	foreach my $value (@GOF){
		print "$id\t$value\n";
	}
}
