#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;
usage:
	perl $0  1.stringtie.gtf 2.stringtie.gtf 3.stringtie.gtf ...   > out.matrix

USAGE


my (%gene, %hash,@name);


opendir DIR, "$ARGV[0]"  or die "";

my @dir=readdir DIR;
#my @dir=@ARGV;

foreach my $file (sort @dir){
	next unless $file =~ m/gtf$/;
	#push @name, $b;
	open IN , "$file" or die "";
	while (<IN>){
		next if (m/^#/);
		chomp;
		my @line= split /\t/,$_;
		next unless ($line[2] eq "transcript");
		die unless $line[-1] =~ /TPM "([^"]*)";/;
		my $tpm=$1;
		#gene_id "BolC1t00004H";
		die unless $line[-1] =~ /gene_id "([^"]*)";/;
		my $geneid=$1;
		$geneid=~s/.TAIR10//g;
		#$hash{$geneid}=$tpm;
		#next unless exists $gene{$geneid};
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
                $sample=~ s/.gtf//g;
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
		my $a = $sample	;
		$a=~s/.gtf//g;
		#print "$a,$hash{$geneid}{$sample}\n";
		print "\t$hash{$geneid}{$sample}";
	}
	print "\n";
}







