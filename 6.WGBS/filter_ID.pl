#!/usr/bin/perl -w

use strict;

open ID, "$ARGV[0]" or die "";
my %ID;
while (<ID>){
	chomp;
	$ID{$_}=1;
}

close ID;

open IN, "$ARGV[1]" or die "";

while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	my @a=split/_/, $line[3];
	next unless exists $ID{$a[1]};
	print "$_\n";
}

close IN;


