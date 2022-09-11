#!/usr/bin/perl -w

use strict;

my $MC=shift;
my $LI=shift;

open MC, "$MC" or die "";
my %MC;
while (<MC>){
	chomp;
	my @line=split/\t/, $_;
	$MC{$line[0]}{$line[1]}="$line[3]\t$line[4]\t$line[5]";
}
close MC;

#Liftover tad output = LI
open LI, "$LI" or die "";
while (<LI>){
	chomp;
	my @line = split /\t/, $_;
	next if ( ! exists $MC{$line[13]}{$line[14]} ) ;
	print "$MC{$line[13]}{$line[14]}\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[11]\t$line[12]\t$line[13]\t$line[14]\t$line[15]\t$line[16]\n";
}

