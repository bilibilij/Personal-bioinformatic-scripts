#!/usr/bin/envs perl -w

use strict;

my $usage=<<USAGE;
###usage:

perl $0 target_scaffold_name_list.txt scaffolding.fasta > result.fa


###Explanation:

target_scaffold_name_list.txt:
---Start---
HiC_scaffold_19
HiC_scaffold_20
---End---

scaffold.fasta
---Start---
>HiC_scaffold_1
...
...
>HiC_scaffold_19
...
>HiC_scaffold_20
...
---End---

result.fa
---Start---
>HiC_scaffold_19_1
...
>HiC_scaffold_19_2
...
>HiC_scaffold_19_3
...
...
---End---

USAGE

die "$usage\n" unless @ARGV == 2;

my $list=shift;
my $fa=shift;

my %name;

open IN, "$list" or die "";

while (<IN>){
	chomp;
	$name{$_}=1;
}

close IN;

my ($seqID, %seq);
open FA, "$fa" or die "";

while (<FA>){
	chomp;
	if (m/^>(\S+)/){$seqID=$1;}
	else{$seq{$seqID}.=$_;}
}

close FA;

foreach my $tar_ID (sort keys %name){
	die "target_id can't be found in scaffolding.fa, please check Names of target_list\n" unless my $seq=$seq{$tar_ID};
	my $N=("N") x 500;
	$seq=~ s/$N/1/g;
	my @fragment=split /1/, $seq;
	my $num=1;
	foreach my $contig (@fragment){
		next if $contig eq "";
		print ">${tar_ID}_${num}\n$contig\n";
		$num++;
	}
}







