#!/usr/local/bin/perl -w

use strict;

my $bed=shift;
my $pass=shift;
my $domain_length=0;
my $qry;
my $chr;
open IN, "$bed" or die "";
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	$domain_length=$line[2]-$line[1];
	$qry=$_;
	$chr=$line[0];
}

close IN;

my %hash;
open IN, "$pass" or die "";
while (<IN>)
{
	chomp;
	my @line=split/\t/, $_;
	next if ("$chr" ne "$line[0]");
	my $lift_domain_len=$line[2]-$line[1];
	my $bound_len=$line[6]-$line[5];
	my $keys="$line[0]\t$line[1]\t$line[2]";
	$hash{$keys}{$line[3]}="$_\t$lift_domain_len\t$bound_len\t$qry";
	#print "$hash{$keys}{$line[3]}\n";
}

foreach my $key1 (keys %hash)
{

	#my $A="fail";
	#my $B="fail";
	my @check;
	foreach my $keys2 (keys %{$hash{$key1}} ) 
	{
	
		my @line=split/\t/, $hash{$key1}{$keys2};
		if ( 0.5< $line[8]/$domain_length <1.5 && $line[9] < 40000)
		{
			push @check, $keys2;
		}
	}
	if ($#check == 1){print "$hash{$key1}{1}\n$hash{$key1}{2}\n";}

}


