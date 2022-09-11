#!/usr/bin/perl -w

open IN, "$ARGV[0]" or die "";
my %hash;
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	$hash{$line[0]}=$line[1];
}
close IN;


open IN, "$ARGV[1]" or die "";
while (<IN>){
	chomp;
	print "$_\n" if $hash{$_} > 200 ;

}

close IN;
