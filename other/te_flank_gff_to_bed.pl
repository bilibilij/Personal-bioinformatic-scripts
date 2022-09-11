#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;
Usage:
	perl $0  gene.gff  flanking_num  >flanking.bed


USAGE

if (@ARGV ==0  ) {die $usage;}

my $gff=shift;
my $flank=shift;
my $prefix=shift;

#open UP, ">$prefix.upstraem.bed" or die "";
#open DO, ">$prefix.downstream.bed" or die "";
open GFF, "$gff" or die "";

while (<GFF>){
	chomp;
	next if m/^#/;
	next if m/^$/;
	die "spliced gff file \n" unless  my @line=split/\t/, $_;
	next if $line[4]-$line[3]+1 <200;
	#next unless $line[2] eq "gene";
	die "err with extarct gene_name from gff\n" unless m/ID=([^;]*);/;
	my $gene_id=$1;
	die "err with extarct gene_name from gff\n" unless m/Classification=([^;]*);/;
	my $class=$1;
	my ($gs, $ge);
	$gs=$line[3]-1;
	$ge=$line[4]-1;
	#print "$line[0]\t$gs\t$ge\t$class:${gene_id}_$line[6]\n";
	if ($line[6] eq "+" ||$line[6] eq "?" || $line[6] eq "." ){
		my ($up_start, $up_end, $down_start, $down_end);
		$up_start=$line[3]-1-1-$flank;
		$up_end=$line[3]-1-1;
		$down_start=$line[4]-1+1;
		$down_end=$line[4]-1+1+2000;
		$up_start=0 if $up_start <=0;
		$up_end=0 if $up_end <=0;
		$down_end=0 if $down_end<=0;
                $down_start=0 if $down_start<=0;
	        print "$line[0]\t$gs\t$ge\t$class:${gene_id}_+\n";

		print "$line[0]\t$up_start\t$up_end\t$class:${gene_id}_+_up\n";
		print "$line[0]\t$down_start\t$down_end\t$class:${gene_id}_+_down\n";
	}
	elsif($line[6] eq "-"){
                my ($up_start, $up_end, $down_start, $down_end);
		$up_start=$line[4]-1+1;
		$up_end=$line[4]+1+2000-1;
		$down_start=$line[3]-1-1-2000;
		$down_end=$line[3]-1-1;
                $up_start=0 if $up_start <=0;
                $up_end=0 if $up_end <=0;
		$down_end=0 if $down_end<=0;
		$down_start=0 if $down_start<=0;
	        print "$line[0]\t$gs\t$ge\t$class:${gene_id}_-\n";

                print "$line[0]\t$up_start\t$up_end\t$class:${gene_id}_-_up\n";
                print "$line[0]\t$down_start\t$down_end\t$class:${gene_id}_-_down\n";
	}
	else { die "error with bed";}
}

close GFF;
#close UP;
#close DO;
#Chr01   .       gene    2894    5929    0.23    -       .       ID=Ovio00001;Name=Ovio00001;

