#!/usr/bin/env perl -w

use strict ;

my $usage=<<USAGE;

perl $0 <match_anno.list> <collapsed.trans.iso.fa> <ref.transcript.fa>

USAGE

die "$usage\n" unless @ARGV == 3;


my $list = shift;
my $col = shift;
my $ref = shift;

open IN, "$col" or die "";

my ($seqID,%col);
while (<IN>){
	chomp;
	if (m/^>([^|]*)/){$seqID=$1;}
	else{$col{$seqID}.=$_;}
}

close IN;


open IN, "$list" or die "";
my %un;
my $Novel=1;

while (<IN>){
	chomp;
	my @line=split/\s+/, $_;
	if ($line[2] =~ m/no_genes_found;/){
		my $id= sprintf "%#07s", $Novel;
		$id = "Novelgene$id";
		$Novel++;
		print ">${id}|$line[0]\n$col{$line[0]}\n";
	}else{
		my @a=split/;/, $line[2];
		foreach my $a (@a){
			$un{$a}=1;
		}
		print ">$line[0]\n$col{$line[0]}\n";
	}
}

close IN;

open IN, "$ref" or die "";
my %ref;
while (<IN>){
        chomp;
        if (m/^>(\S+)/){$seqID=$1; $seqID =~ s/T/G/; }
        else{$ref{$seqID}.=$_;}
}

close IN;


foreach my $ID (sort keys %ref){
	#$ID=~s/\..*//g;
	next if $ID =~ m/\./;
	next if exists $un{$ID};
	print ">$ID\n$ref{$ID}\n";
}


