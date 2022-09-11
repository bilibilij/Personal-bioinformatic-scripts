#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;

perl $0 blasr.out5


USAGE


my %hash;
open IN, "$ARGV[0]" or die "";

while (<IN>){
	chomp;
	next if m/^qName\s/;
	my @line=split/\s+/,$_;
	my $Match=$line[11];
	my $Mism=$line[12];
	my $Ins=$line[13];
	my $Del=$line[14];
	my $pid=$Match/($Mism+$Match+$Ins+$Del);
	$hash{$line[0]} = 0 unless exists $hash{$line[0]};
	$hash{$line[0]} = $pid  if ($pid > $hash{$line[0]} );
}
close IN;

foreach my $trans (sort keys %hash){
	print "$trans\t$hash{$trans}\n";
}


