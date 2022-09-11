#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;

perl $0 matchannot.out collapsed.group.txt > structral_id.list 

USAGE

die "$usage\n" unless  @ARGV == 2;

my %MA;
open MA, "$ARGV[0]" or die "";
while (<MA>){
	chomp;
	next unless m/^result/;
	my @line=split/\s+/, $_;
	$MA{$line[1]}=$line[2];
}
close MA;


open CO, "$ARGV[1]" or die "";
while (<CO>){
	chomp;
	my @line=split/\s+/,$_;
	my @a=split/,/, $line[1];
	my %temp;
	foreach my $c ( @a ){
		$temp{$MA{$c}}=1 unless exists $temp{$MA{$c}};
	}

	my @temp = keys %temp;
        my $pass;
	my $check=0;
        if (@temp == 1){$pass="Unique";}
        elsif(@temp == 2){
		if ($temp[0] eq "no_genes_found" || $temp[1] eq "no_genes_found"){
			$pass="Unique";
			$check=1;
		}else{
			$pass="multi";
		}
	}else{
		$pass="multi";
	}

	print "$line[0]\t$pass\t";

	foreach my $key (keys %temp){
		next if $check==1 && $key eq "no_genes_found";
		print "$key;";
	}
	print "\t$line[1]\n";
}
close CO;




