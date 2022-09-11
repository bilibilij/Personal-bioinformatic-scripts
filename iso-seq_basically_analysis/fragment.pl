#!/usr/bin/perl -w

use strict;

my $usage =<<USAGE;

perl $0 <bam> <sample_id> >stat

USAGE

die "$usage\n" unless @ARGV == 2;

#Dependencies check and parameters read;
my $software_info = `samtools 2>&1`;
print STDERR "samtools: PASS\n " if  ($software_info =~ m/samtools/) ;
die "samtools not in PATH\n" unless ($software_info =~ m/samtools/);

my $bam = $ARGV[0];
my $sample= $ARGV[1];

#processing

my ($Total_reads_num, $Total_reads_len, $MinLen, $MaxLen, $MeanLen, $N50_len)=(0,0,0,0,0,0);

my (%hash, %qv);

open IN, "samtools view $bam | " or die "";

while (<IN>){
	chomp;
	next if m/^#/;
	next if m/^@/;
	die "can not detect quality \n " unless m/rq:f:(\S+)\s/;
	$qv{$1}++;

	my @line=split/\t/, $_;
	my $len=length($line[9]);
	$hash{$len}++;
	$Total_reads_num++;
	$Total_reads_len=$Total_reads_len+$len;
	$MaxLen =$len  if $len > $MaxLen; 
	print $_ if $len == 1;
}

close IN;


open OUT, ">$sample.fragment.txt" or die "";
my $half = int($Total_reads_num/2);
my $zero=0;
foreach my $num (sort {$a <=> $b} keys %hash) {
	$MinLen = $num if $MinLen == 0 ;
	print OUT "$num\t$hash{$num}\n";
	$N50_len = $num if (  $zero < $half <= ($zero+$hash{$num})) ;
	$zero+=$hash{$num};
}
close OUT;

$MeanLen = $Total_reads_len/$Total_reads_num;
$Total_reads_len = sprintf "%.2f",$Total_reads_len/1000000000;

print "SampleID\tTotal_number\tTotal_bases(Gbp)\tMinimum_length\tAverage_length\tMaximum_length\tN50\n";
print "$sample\t$Total_reads_num\t$Total_reads_len\t$MinLen\t$MeanLen\t$MaxLen\t$N50_len\n";

open OUT, ">$sample.qv.txt" or die "";
foreach my $qv (sort {$a <=> $b}  keys %qv){
	print OUT "$qv\t$qv{$qv}\n";
}

close OUT;






