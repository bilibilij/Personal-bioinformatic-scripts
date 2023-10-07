#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;

perl $0 out(OV_ATAC_S.geneDistance.closest.ACR) acr OV_ATAC_S.geneDistance.closest.ACR.bed 

perl 2.Private_but_Conserved.pl.review OV_ATAC_S.geneDistance.closest.ACR.bed - ov.gff /usr_storage/jcf/2.zhugecai/00.diploid/10.ks/2.wgdi/at_ov_0.2.bi.cor.csv /usr_storage/jcf/2.zhugecai/00.diploid/Ovio_1200.genome.fasta OV_ATAC_S.geneDistance.closest.ACR
USAGE

die "$usage\n" unless @ARGV==3;

my $out=$ARGV[0];
my $acr=$ARGV[1];
open IN, "$acr" or die "";
open COM, ">$out.ACR.intersect.com" or die "";
while (<IN>){
        chomp;
	#print OUT "$_\n" if m/share/;
	next if m/share/;
	#next unless m/private/;
        my @line=split/\t/, $_;
	
	#my $last="private";
        open BL, "$out.conserved.tmp/$line[3]:$line[5].out6" or die "";
        open BE, ">$out.conserved.tmp/$line[3]:$line[5].realn.bed" or die "";

	while (<BL>){
                chomp;
                my @arr=split/\t/, $_;
		#$last="conserved";
		#::Chr05:26578394-26590493 
		die "match err \n " unless m/^::([^:]*):(\d+)-(\d+)\t/;
		my ($chr, $ss, $end) = ($1, $2,$3);

		if ($arr[6] < $arr[7] ) { $ss+=$arr[6];$end+=$arr[7];}
		else{ $ss+=$arr[7];$end+=$arr[6];}
		
		print BE "$chr\t$ss\t$end\n";

        }
        close BL;
	close BE;
	print COM "sort -k1,1 -k2,2n -k3,3n $out.conserved.tmp/$line[3]:$line[5].realn.bed | bedtools intersect  -a - -b $ARGV[-1] > $out.conserved.tmp/$line[3]:$line[5].realn.bed.intersect 2>/dev/null \n";
	#`sort -k1,1 -k2,2n -k3,3n $out.conserved.tmp/$line[3]:$line[5].realn.bed | bedtools intersect  -a - -b $ARGV[-1] > $out.conserved.tmp/$line[3]:$line[5].realn.bed.intersect 2>/dev/null `;
	#my $num=`wc -l $out.conserved.tmp/$line[3]:$line[5].realn.bed.intersect|awk '{print \$1}' `; 
	#$last="shared" if $num>0;
	#$line[-1]=$last;
	#my $out=join("\t", @line);
	#print OUT "$out\n";

}
close IN;
close COM;

`ParaFly -c $out.ACR.intersect.com -CPU 100 `;


open OUT, ">$out.ACR.conserved.Final.bed" or die "";
open IN, "$acr" or die "";
while (<IN>){
        chomp;
        print OUT "$_\n" if m/share/;
	next if m/share/;
	#next unless m/private/;
        my @line=split/\t/, $_;
	my $last="private";
        open BL, "$out.conserved.tmp/$line[3]:$line[5].out6" or die "";

	while (<BL>){
                chomp;
                my @arr=split/\t/, $_;
                $last="conserved";
	}
	close BL;
	my $num=`wc -l $out.conserved.tmp/$line[3]:$line[5].realn.bed.intersect|awk '{print \$1}' `;
	$last="shared" if $num>0;
	$line[-1]=$last;
	my $out=join("\t", @line);
	print OUT "$out\n";

}
close OUT;
close IN;



