#!/usr/bin/envs perl -w

use strict;

my $usage=<<USAGE;

#usage:

perl $0 contig_name_for_scaffolding.list total_contig.fa > after_joint.fa

#explanation

contig_name_for_scaffolding.list
---start---
HiC_scaffold_21_1 HiC_scaffold_21_2 
HiC_scaffold_22_2 HiC_scaffold_22_4 HiC_scaffold_1202_7
---end---


USAGE

die "$usage\n" unless @ARGV == 2;

my %hash;
my $ad=\%hash;
my @name;

open ID, "$ARGV[0]" or die "";

while (<ID>){
	chomp;
	my @line = split/\s+/, $_;
	my $name= shift @line;
	push @name, $name;
	push @{$ad -> {$name}}, @line;


}

close ID;


my ($seqID, %seq);
open FA, "$ARGV[1]" or die "";

while (<FA>){
	chomp;
	if (m/>(\S+)/){$seqID=$1;}
	else{$seq{$seqID}.=$_;}
}

close FA;

my $N=("N") x 500;
#$seq=~ s/$N/1/g;

foreach my $id (keys %hash){
	my @line = @{$hash{$id}};
	my @seq;
	foreach my $contig ( @line ){
		if ($contig=~ m/:R$/){
			$contig =~s/:R//;
			my $seq = $seq{$contig};
			$seq=~ tr/ATCGagtc/TAGCtcag/;
			$seq=reverse ($seq);
			push @seq, $seq;
		}else{
			my $seq = $seq{$contig};
			push @seq, $seq;
		}
	}
	my $scaffold=join "$N", @seq;
	print ">$id\n$scaffold\n";
}







