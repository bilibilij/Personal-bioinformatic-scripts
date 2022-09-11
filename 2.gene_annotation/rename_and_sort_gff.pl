#!/usr/bin/perl -w 

my $usage=<<USAGE;

perl $0 gff  prefix > sorted.gff3

Author: Changfu Jia
Welcome to contact 1020160171\@qq.com

USAGE


use strict;

open GF, "$ARGV[0]" or die "";

my $h=$ARGV[1];
my $m=0;
my %hash;
my %gene;
my $start;
while (<GF>){
	chomp;
	next if m/^$/;
	next if m/^#/;
	my @line = split/\t/, $_;
	pop @line;
	die $! if @line != 8;
	my $c=join("\t", @line);
	if ($line[2] eq "gene") {
		$start=$line[3];
		$m=0;
		$gene{$line[0]}{$start}="$c";
	}elsif( $line[2] eq "mRNA"){
		$m++;
		#@{$hash{$line[0]}{$start}->{$m}} ="$c";
		$hash{$line[0]}{$start}{$m}.="$c\n";
		#$m++;
	}else{
		$hash{$line[0]}{$start}{$m}.="$c\n";
	}
}

close GF;


my $gene_num=1;
foreach my $chr (sort keys %hash){
	foreach my $st (sort {$a<=>$b} keys  %{$hash{$chr}}){
		#gene line printf 
		my $s = (sprintf  "%05d", $gene_num);
		my $gene_id="$h$s";
		print "$gene{$chr}{$st}\tID=$gene_id;Name=$gene_id;\n";
		$gene_num++;

		#each line print with its parent mRNA;
		my $mRNA_num=1;
		foreach my $mr (sort {$a <=> $b} keys  %{$hash{$chr}{$st}} ){
			
			my $mRNA="$gene_id.t$mRNA_num";
			die "err\n" unless my @arr =split /\n/, $hash{$chr}{$st}{$mr};
			my %num;
			foreach my $fea (@arr) { 
				my @line=split/\t/ , $fea;

				print "$fea\tID=$mRNA;Parent=$gene_id;Name=$mRNA;\n" if $line[2] eq "mRNA";
				next if $line[2] eq "mRNA";
				$num{$line[2]}++;
				my $N=$num{$line[2]};
				print "$fea\tID=$mRNA.$line[2]$N;Parent=$mRNA;\n";
			}
			$mRNA_num++;
		}
		print "\n";
	}
}
