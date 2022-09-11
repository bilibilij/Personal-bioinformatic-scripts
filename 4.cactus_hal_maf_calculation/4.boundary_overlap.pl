#!/usr/local/bin/perl -w

use strict;
 


#23 910 
my %hash;

open IN, "$ARGV[0]" or die "";

while (<IN>){
	chomp;
	my @line=split/\t/,$_;
	my $boundA1=$line[1]-20000;
	my $boundA2=$line[1]+20000;
	my $boundB1=$line[2]-20000;
	my $boundB2=$line[2]+20000;
	my $key1="$line[3]-"."$line[4]-"."$line[5]";
	if (  ( $boundA1  <= $line[8] <= $boundA2) || ($boundA1 <=$line[9] <= $boundA2) || ($line[8] <= $boundA1 <=$boundA2 <= $line[9] ) ) {
		$hash{$key1}{$line[6]}="$_\tApass";
		#print "$_\tApass";
	}
	elsif ( ($boundB1 <= $line[8] <= $boundB2) || ($boundB1 <= $line[9] <= $boundB2) ||  ($line[8] <= $boundB1 <=$boundB2 <= $line[9])) 
		{
			#print "$_\tBpass";
		$hash{$key1}{$line[6]}="$_\tBpass";
		}
	else
       	{#print "$_\tfail";
		$hash{$key1}{$line[6]}="$_\tfail";
	}
}

foreach my $key1 (keys %hash){
	
	my $checkA="fail";
	my $checkB="fail";
	foreach my $key2 (keys %{$hash{$key1}}){
		my @line=split/\t/, $hash{$key1}{$key2};
		if ($line[17] eq "Apass" ) {$checkA="pass";}
		if ($line[17] eq "Bpass" ) {$checkB="pass";}
	}

	if ($checkA eq "pass" && $checkB eq "pass") {print "$hash{$key1}{1}\n$hash{$key1}{2}\n";}
}
