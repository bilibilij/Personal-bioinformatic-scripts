#!/usr/bin/perl -w

use strict;

open IN, "$ARGV[0]" or die "";

my %hash;
while (<IN>)
{
	chomp;
	my @line=split/\t/, $_;
	my $mark=$line[0].$line[1].$line[2];
	$hash{$mark}{$line[3]}="$_\n";
}


foreach my $keys (keys %hash)
{
	my $A="fail";
	my $B="fail";
	foreach my $bound (keys %{$hash{$keys}})
	{
		if ($bound==1) { $A="pass";}
		if ($bound==2) { $B="pass";}
	}
	if ($A eq "pass" && $B eq "pass"){print "$hash{$keys}{1}$hash{$keys}{2}";}
}
