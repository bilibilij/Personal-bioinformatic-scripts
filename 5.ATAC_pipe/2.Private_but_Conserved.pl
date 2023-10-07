#!/usr/bin/perl -w

use strict;

my $usage=<<USAGE;

perl $0 ACR.bed  dup.tsv wgdi.gff  at_ov_0.2.bi.cor.csv  genome.fasta out


===: gene
...: syntenic gene pairs
+++: region for blast 

----====---ACR----===------
   .   .  +   +   .  .              
  .   . +       +  .  .             
 .   .+           + .   .   
-====------===-------====--------
    |                |
    ——————————————————
    This area is for non-shared ACR to blast ; 


ACR.bed is be like this 
Chr04   344806  345004  Ovio04459_Ovio13439     Ovio13439       Ovio13439:2     0       -       -       private
Chr04   311688  311925  Ovio04462_Ovio13435     Ovio13435       Ovio13435:1     1200    -       -       private
Chr01   99466742        99467173        Ovio04466_Ovio13433     Ovio04466       Ovio04466:1     747     Ovio13433:1     1       shared
Chr04   305528  305928  Ovio04466_Ovio13433     Ovio13433       Ovio13433:1     1038    Ovio04466:1     1       shared
Chr01   99487382        99488307        Ovio04472_Ovio13427     Ovio04472       Ovio04472:1     0       -       -       private
Chr04   258834  259359  Ovio04473_Ovio13426     Ovio13426       Ovio13426:1     1638    -       -       private


wgdi.gff
Chr01   Ovio00001.t1    2894    5929    -       1       Ovio00001
Chr01   Ovio00002.t1    86775   87242   -       2       Ovio00002
Chr01   Ovio00003.t1    138249  140154  -       3       Ovio00003
Chr01   Ovio00004.t1    214567  215730  +       4       Ovio00004
Chr01   Ovio00005.t1    251592  252816  -       5       Ovio00005
Chr01   Ovio00006.t1    255358  256584  -       6       Ovio00006
Chr01   Ovio00007.t1    287913  289580  -       7       Ovio00007
Chr01   Ovio00008.t1    290111  291205  +       8       Ovio00008
Chr01   Ovio00009.t1    302784  305610  -       9       Ovio00009


USAGE

die "$usage\n" unless @ARGV == 6;

my $acr=shift;
my $dup=shift;
my $gff=shift;
my $cor=shift;
my $fa=shift;
my $out=shift;



my %gff;
my %gene_bed;
open IN, "$gff" or die "";
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	$gff{$line[-1]}="$line[0]\t$line[-2]";
	$gene_bed{"$line[0]\t$line[-2]"}="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[-1]";
	$gene_bed{$line[-1]}="$line[0]\t$line[2]\t$line[3]\t$line[4]";
}
close IN;

my (%dup_ind, %dup_check);
open IN, "$dup" or die "";
while (<IN>){
	chomp;
	next if m/^ref/;
	next unless m/wgd/;
	my @line=split/\t/, $_;
	$dup_ind{$gff{$line[1]}}=1;
        $dup_ind{$gff{$line[2]}}=1;
	$dup_check{"$line[1]_$line[2]"}=1;
}
close IN;


#my (%forward, %reverse);
my %cor;
open IN, "$cor" or die "";
while (<IN>){
	chomp;
	my @line=split/,/, $_;
	$cor{$line[0]}="$line[1],$line[2],";
	my @gene_chain= split/_/, $line[-6];
	my @temp;
	for (my $i=0;$i<=$#gene_chain;$i++){
		next unless exists $dup_ind{"$line[2]\t$gene_chain[$i]"};
		push @temp, $gene_chain[$i];
	}
	$cor{$line[0]}.=join("_", @temp);
	#print "$cor{$line[0]}\n";
	#my $a=join("\n", @temp);	
	#print "$a\n";
}
close IN;


my (%up, %down);
foreach my $block (sort keys %cor ){
        my ($rchr, $qchr, $chain)=split/,/ , $cor{$block};
        my @gene_chain=split/_/, $chain;
        #next if scalar @gene_chain <= 1;
        for (my $i=0;$i<=$#gene_chain;$i++){
                #$up;
                my ($Chr, $ss, $end, $strand, $gene)= split/\t/, $gene_bed{"$qchr\t$gene_chain[$i]"};
                if ( 0 < $i < $#gene_chain){
                        my ($upChr, $upStart, $upEnd, $upStrand, $upGene)= split/\t/, $gene_bed{"$qchr\t$gene_chain[$i-1]"};
                        #my ($Chr, $ss, $end, $strand, $gene)= split/\t/, $gene_bed{"$qchr\t$gene_chain[$i]"};
                        my ($dwChr, $dwStart, $dwEnd, $dwStrand, $dwGene)= split/\t/, $gene_bed{"$qchr\t$gene_chain[$i+1]"};
                        $upEnd-=10000 if ($upEnd >= $ss - 100 ); #if two gene was two close ,we consider a much bigger promoter;
                        $dwStart+=10000 if ($end >= $dwStart - 100 );
                        $up{$gene}="$Chr\t$upEnd\t$ss";
                        $down{$gene}="$Chr\t$end\t$dwStart";
                }elsif ($i == 0 ){
                        my $upS=$ss-10000;
                        if ($#gene_chain == 0){ $up{$gene}="$Chr\t$upS\t$ss"; my $dwStart = $end+10000; $down{$gene}="$Chr\t$end\t$dwStart";   }
                        else{
                        #print "$i\t$#gene_chain\n";
                        #print "$qchr\t$gene_chain[$i+1]\n";
                        my ($dwChr, $dwStart, $dwEnd, $dwStrand, $dwGene)= split/\t/, $gene_bed{"$qchr\t$gene_chain[$i+1]"};
                        $upS-=10000 if ($upS >= $ss - 100 );
                        $dwStart+=10000 if ($end >= $dwStart - 100 );
                        $up{$gene}="$Chr\t$upS\t$ss";
                        $down{$gene}="$Chr\t$end\t$dwStart";
                        }
                }elsif ($i == $#gene_chain){
                        my $dwE=$end+10000;
                        my ($upChr, $upStart, $upEnd, $upStrand, $upGene)= split/\t/, $gene_bed{"$qchr\t$gene_chain[$i-1]"};
                        $upEnd-=10000 if ($upEnd >= $ss - 100 );
                        $dwE+=10000 if ($end >= $dwE - 100 );
                        $up{$gene}="$Chr\t$upEnd\t$ss";
                        $down{$gene}="$Chr\t$end\t$dwE";
                }


        }
}



`mkdir $out.conserved.tmp` unless -e "$out.conserved.tmp";

my %ACR;

open BED, ">$out.get_conservedFa.com" or die "";
open COM, ">$out.blast.Conserved.com" or die "";
open IN, "$acr" or die "";
while (<IN>){
	chomp;
	next unless m/private/;
	my @line=split/\t/, $_;
	next unless exists $dup_check{$line[3]};
	my ($g1, $g2) = split/_/, $line[3];
	my $counter_gene;

	$ACR{$line[3]}{$line[5]}=$_;

	if ($g1 eq $line[4]){ $counter_gene=$g2;}
	else{$counter_gene=$g1;}
	#$gene_bed{$line[-1]}="$line[0]\t$line[2]\t$line[3]\t$line[4]";
	my ($chr, $ss, $end, $strand) = split/\t/, $gene_bed{$line[4]}; 


	open OUT, ">$out.conserved.tmp/$line[3]:$line[5].bed" or die "";
	open OUT2, ">$out.conserved.tmp/$line[3]:$line[5].counter.bed" or die "";
	my $counter_bed="";
	my ($Cchr, $Css, $Cend, $Cstrand) = split/\t/, $gene_bed{$counter_gene};
	if ($line[6] == 0 ) {
		#my ($Cchr, $Css, $Cend, $Cstrand) = split/\t/, $gene_bed{$counter_gene};
		$Css-=5000;
		$Cend+=5000;
		$counter_bed="$Cchr\t$Css\t$Cend\n";
	}else{

			if (not exists $up{$counter_gene}){
                                $Css-=10000;
				$Css = 0 if ($Css<0);
				#$Cend+=10000;
                                $counter_bed.="$Cchr\t$Css\t${Cend}\n";
                        }else{
			$counter_bed.="$up{$counter_gene}\n";	
			}
                        if (not exists $down{$counter_gene}){
				#$Css-=10000;
                                $Cend+=10000;
				$Css = 0 if ($Css<0);

                                $counter_bed.="$Cchr\t$Css\t$Cend\n";}
                        else{
	
				$counter_bed.="$down{$counter_gene}\n";
			}
		#die "$line[0]\t$line[1]\t$line[2]\t$chr\t$ss\t$end\n$_\n" if not exists $down{$counter_gene};
	}
	print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]:$line[5]\n";
	
	#if ($counter_bed eq ""){print "$_\n"; }#print $up{$line[4]}; }

	print OUT2 "$counter_bed";
	close OUT;
	close OUT2;
	print BED "bedtools getfasta -fi $fa -bed $out.conserved.tmp/$line[3]:$line[5].bed -name >  $out.conserved.tmp/$line[3]:$line[5].fa\n";
        print BED "bedtools getfasta -fi $fa -bed $out.conserved.tmp/$line[3]:$line[5].counter.bed -name >  $out.conserved.tmp/$line[3]:$line[5].counter.fa\n";
	print COM "makeblastdb -in $out.conserved.tmp/$line[3]:$line[5].fa -dbtype nucl -out $out.conserved.tmp/$line[3]:$line[5]  2>&1>/dev/null ;    blastn -db  $out.conserved.tmp/$line[3]:$line[5] -query  $out.conserved.tmp/$line[3]:$line[5].counter.fa  -outfmt 6 -out $out.conserved.tmp/$line[3]:$line[5].out6  -task blastn-short -evalue 1e-10   2>&1>/dev/null \n";
}
close IN;

`ParaFly -c $out.get_conservedFa.com -CPU 40 `;
`ParaFly -c $out.blast.Conserved.com -CPU 40 `;

open OUT, ">$out.ACR.conserved.bed" or die "";
open IN, "$acr" or die "";
while (<IN>){
        chomp;
	print OUT "$_\n" if m/share/;
        next unless m/private/;
        my @line=split/\t/, $_;
        next unless exists $dup_check{$line[3]};
        my ($g1, $g2) = split/_/, $line[3];
	my $last="private";
	open BL, "$out.conserved.tmp/$line[3]:$line[5].out6" or die "";
	while (<BL>){
		chomp;
		my @arr=split/\t/, $_;
		#next unless $line[2] > 80;
		$last="conserved";
	}
	close BL;	
	$line[-1]=$last;
	my $out=join("\t", @line);
	print OUT "$out\n";

}
close IN;
close OUT;
