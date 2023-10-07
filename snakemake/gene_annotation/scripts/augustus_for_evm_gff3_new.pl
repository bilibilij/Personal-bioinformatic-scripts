#!/usr/bin/perl -w
use strict;

my $usage=<<USAGE;
usage:
	perl $0 augustus.gff3 > evm.gff3
	this script is for augustus output gff3 to transform to evm input denovo gff3;

USAGE


open GFF, "$ARGV[0]" or die "can not open gff";
open ME, ">$ARGV[0].med" or die "can not write medium.gff3";

while (<GFF>){
        chomp;
        next if (m/^\#/);
        next if (m/\tintron/);
        next if (m/codon/);
        if (m/CDS/){
                my $CDS=$_;
                my $exon=$_;
                $exon=~ s/CDS/exon/g;
                print ME "$exon\n$CDS\n";
        }else{
                print ME "$_\n";
        }
}
open ME, "$ARGV[0].med" or die "can't open medium.gff3";





my ( @mRNA_length,@mRNA,  $line,  %hash);
my $mRNA_length=0;
my @gene=();
while (<ME>) {
        chomp;
        next if (m/^\#/);
        next if (m/CDS/);
        next if (m/\tintron/);
        next if (m/codon/);
	if (m/gene/) { 
		if (@gene){
			$hash{$mRNA_length}=$line;
			push @mRNA_length, $mRNA_length;
			my @mRNA_length1 = sort { $b <=> $a } @mRNA_length;
			print "$hash{$mRNA_length1[0]}";
			@mRNA_length1=();
			@mRNA_length=();
			$mRNA_length=0;
			%hash=();
			@mRNA=();
			$line=undef;
		}else{
			unshift @gene, $_;
			$mRNA_length=0;
		}
                print "$_\n";

	}
	if (m/mRNA/){
		if (@mRNA){
			push @mRNA_length, $mRNA_length;
			$hash{$mRNA_length}=$line;
			$mRNA_length=0;
			$line="$_\n";

		}else{
			$line="$_\n";
			unshift @mRNA, $_;
		}
	}
	if (m/exon/){
		$line.="$_\n";
		my $CDS=$_;
		$CDS=~  s/exon/CDS/g;
		$line.="$CDS\n";
	        my @length=split /\t/, $_;
		$mRNA_length=int($mRNA_length)+abs(int($length[4]-$length[3]+1));
	}
}


$hash{$mRNA_length}=$line;
push @mRNA_length, $mRNA_length;
my @mRNA_length1 = sort { $b <=> $a } @mRNA_length;
print "$hash{$mRNA_length1[0]}\n";











