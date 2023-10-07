#!/usr/bin/perl -w

use strict;

open IN, "$ARGV[0]" or die "";
my ($num, $id);
$id="";
while (<IN>){
	chomp;
	my @line=split/\s+/, $_;
	if ($id eq $line[8] ) {
		$num++;
		print "$line[0]\t$line[1]\t$line[2]\t$line[8]:$num:$line[-1]\n";
	}else{
		$id = $line[8];
		$num=1;
                print "$line[0]\t$line[1]\t$line[2]\t$line[8]:$num:$line[-1]\n";
	}

}
close IN;

