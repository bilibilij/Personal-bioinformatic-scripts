#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;
usage:
	perl $0 id.list DIR> out.matrix

USAGE

if (@ARGV != 2 ) {die $usage;}

my (%gene, %hash,@name);

open ID, "$ARGV[0]" or die "";

while (<ID>){
	chomp;
	$gene{$_}=1;
}

close ID;

opendir DIR, "$ARGV[1]"  or die "";

my @dir=readdir DIR;

foreach my $file (sort @dir){
	next unless $file=~ m/gtf$/;
	push @name, $b;
	open IN , "$ARGV[1]/$file" or die "";
	while (<IN>){
		next if (m/^#/);
		chomp;
		my @line= split /\t/,$_;
		next unless ($line[2] eq "transcript");
		die unless $line[-1] =~ /TPM "([^"]*)";/;
		my $tpm=$1;
		die unless $line[-1] =~ /transcript_id "([^"]*)";/;
		my $geneid=$1;
		#$hash{$geneid}=$tpm;
		next unless exists $gene{$geneid};
		$hash{$geneid}{$file}=$tpm;
	}
	close IN;
}



my $check=0;

foreach my $geneid (sort keys %hash)
{
	next if $check>0;
	print "geneid";
	foreach my $sample (sort keys %{$hash{$geneid}})
	{
		next if $check>0;
		$sample=~ s/.out.gtf//g;	
		print "\t$sample";
	}
	print "\n";
	$check++;
}




foreach my $geneid (sort keys %hash)
{
	print "$geneid";
	foreach my $sample (sort keys %{$hash{$geneid}})
	{

		print "\t$hash{$geneid}{$sample}";
	}
	print "\n";
}







