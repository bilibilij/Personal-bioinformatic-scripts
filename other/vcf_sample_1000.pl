#!/usr/bin/envs perl -w

use strict;
use List::Util qw/shuffle/;


my $usage=<<USAGE;

perl $0 small_vcf big_vcf out_dir sampling_num(1000)

USAGE

die "$usage\n" unless @ARGV == 4;


#restore header and mark small vcf for avoid repeat sampling
my $header;
my %small;
open IN, "$ARGV[0]" or die "";
my $len=0;
while (<IN>){
	chomp;
	if (m/^#/){$header.="$_\n";}
	else{
		my @line=split/\t/, $_;
		#$small{$line[2]}=1;
		$len++;
	}
}

close IN;
#my @num=keys %small;
#my $len=scalar @num;




my %big;
open IN, "$ARGV[1]" or die "";

while (<IN>){
	chomp;
	next if m/^#/;
	my @line=split/\t/, $_;
	next if exists $small{$line[2]};
	$big{$line[2]}=$_;
}
close IN;


for (my $i=1;$i<=$ARGV[3];$i++){
	open OUT, ">$ARGV[2]/sample_1000.$i.vcf" or die "";
	my $check=1;
        print OUT "$header";

	my @big_key = keys %big;
	my @shuff = shuffle @big_key;
	foreach my $site (@shuff){
		last if $check > $len;
		print OUT "$big{$site}\n";
		$check++;
	}
	####shuffle replace 1000 sampleing ;
	#foreach my $site (keys %big){
		#last if $check >$len;
		#print OUT "$big{$site}\n";
		#$check++;
	#}
	close OUT;
}





