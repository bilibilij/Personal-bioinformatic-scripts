#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;

perl $0 idmapping.gz(file format must be gz) swiss/tremble.blast6.tab > swiss.go

USAGE

die $usage if @ARGV == 0;

my $idm=shift;
my $tab=shift;


open TAB, "$tab" or die "";

my %hash;
while (<TAB>){
	chomp;
	my @line=split/\t/, $_;
	$hash{$line[1]}{$line[0]}=1;
}

close TAB;



my %go;

open IN, "gunzip -c $idm | " or die "";

while (<IN>){
	chomp;
	my @line=split/\t/,$_;
	next if  $line[7] eq  "";
	next unless exists $hash{$line[0]};

        my $go = $line[7];
        $go=~s/ //g;
	my @g = split/;/, $go;
	foreach my $gene (keys %{$hash{$line[0]}}){
		foreach my $gg (@g){
			$go{$gene}{$gg}=1;
		}
	}
}

foreach my $gene (sort keys %go){
	foreach my $go (sort keys %{$go{$gene}}){
		print "$gene\t$go\n";
	}
}




