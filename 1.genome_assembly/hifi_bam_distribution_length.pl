#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;

perl $0 hifi.fastq 10000000 

first print HIFI CCS reads N50, and the second parameter is to random sampling how much base to test hifi assembly quality

USAGE


die "$usage\n" unless @ARGV == 2;


open FQ,"gzip -dc $ARGV[0]|" or die "";

my $len;
#my %hash;

my ($seqID, %hash, %seq);

while (<FQ>)
{
	chomp;
	#next if m/^>/;
	if (m/^>(\S+)/ or m/^@(\S+)/ ){
		$seqID=$1;
	}
	else{
		$hash{$seqID}+=length($_);	
		$seq{$seqID}.=$_;
		$len+=length($_);
	}
}

close FQ;



print "total length is $len\n";


my @len = sort {$hash{$a} <=> $hash{$b}}   keys %hash;

my $num = scalar @len;

my ($med, $aver);

if ($num%2==0){
	my ($small, $big) = ($len[($num/2)-1], $len[$num/2]);
	$med = ($hash{$small} + $hash{$big}) / 2 ;
}
else{
	$med = $hash{$len[($num+1)/2]};
}

print "HIFI reads N50 is $med\n";


my $limit=$ARGV[1];
my $check=0;
open OUT, ">$ARGV[0].$ARGV[1].fa" or die "";

foreach my $id (keys %seq){
	next unless $check <= $limit ;

	$check+=$hash{$id};
	print OUT ">$id\n$seq{$id}\n";
}

close OUT;








