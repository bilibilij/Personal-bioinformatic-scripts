#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;

perl $0 gff  

Script for transforming gff3 to exon.bed and gene.bed
Author: Changfu Jia
Welcome to contact 1020160171\@qq.com

USAGE

die "$usage" if @ARGV == 0;



#这部分，gene提取一个bed, CDS提取一个bed， 全基因组提取一个BED, 全基因组的用软件算；

open GF, "$ARGV[0]" or die "";

open GE, ">gene.bed" or die "";

open CD, ">exon.bed" or die "";

while (<GF>){
	chomp;
	next if (m/^$/);
	my @line =split/\s+/ , $_;
	if ($line[2] eq "gene"){
		die "" unless $line[-1] =~ m/ID=([^;]*);/;
		my $st= $line[3]-1;
		my $end= $line[4] -1;
		print GE "$line[0]\t$st\t$end\t$1\n";
	}
	elsif ($line[2] eq "exon") { 
		die "" unless $line[-1] =~ m/ID=([^;]*);/;
                my $st= $line[3]-1;
                my $end= $line[4] -1;
		print CD "$line[0]\t$st\t$end\t$1\n";
        }
	else {next;}

}

close GF;



