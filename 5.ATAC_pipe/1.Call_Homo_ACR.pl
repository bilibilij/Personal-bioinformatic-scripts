#!/usr/bin/envs perl -w  

use lib("/home/jcf/scripts/jcf_perlib/");
use lsyPig;
#use strict;

my $usage=<<USAGE;

perl $0 dup_tbl ACR_gene_assign genome.fasta out all

USAGE

die "$usage\n" unless @ARGV == 5;

#export PATH="$PATH:/software//ncbi-blast-2.10.1+/bin/"
#check dependencies
print STDERR "\n============================================\n";
print STDERR "Detecting the dependency softwares:\n";
check_dependency( "blastn -h 2>&1" , "blastn"  );
check_dependency( "ParaFly 2>&1 ", "ParaFly"   );

print STDERR "Dependencies:\tOK\n";
print STDERR "====================================\n\n\n";

my $dup=shift;
my $ACR=shift;
my $fa=shift;
my $out=shift;
my $step=shift;

#if ($step==1 || $step eq "all"|| $step eq "debug" ){


#if ($step ne "debug" ){
#if ( -e "./$out.tmp"){
#	`/usr/bin/rm -rf ./$out.tmp ;  `;
#	`mkdir ./$out.tmp`;

#}else{
	`mkdir ./$out.tmp`;
#}
#}

my %ACR;
my %ACR_dis;
open IN, "$ACR" or die "";
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	#my $pos;
	#if ($step eq "debug"){$pos="$line[0],$line[1],$line[2]";}
	#else{$pos="$line[0]\t$line[1]\t$line[2]";}
	my $pos="$line[0],$line[1],$line[2]";
	my @acr=split/:/, $line[-1];
	my $gene=$acr[0];
	$ACR{$gene}{$pos}="$acr[0]:$acr[1]";
	$ACR_dis{$gene}{$pos}=$acr[-1];
}

close IN;

my $seqID;
my %seq;
open IN, "$fa" or die "";
while (<IN>){
	chomp;
	if (m/^>(\S+)/){$seqID=$1;}
	else{$seq{$seqID}.=$_;}
}

close IN;

#open COM, ">blast.com" or die "";

open OUT, ">$out.dupACR.txt" or die $!;
open IN, "$dup" or die "";
while (<IN>){
	chomp;
	next if m/^ref/;
	my @line=split/\t/, $_;
	my $task="$line[1]_$line[2]";
	#if ($step eq "debug"){
	my @ACR1=sort keys %{$ACR{$line[1]}};
	my @ACR2=sort keys %{$ACR{$line[2]}};
	my $ACR1=join(";", @ACR1);
	my $ACR2=join(";", @ACR2);
	my $nACR1=scalar @ACR1;
	my $nACR2=scalar @ACR2;
	print OUT "$line[1]\t$ACR1\t$line[2]\t$ACR2\t$nACR1\t$nACR2\n";
}
close IN;
close OUT;

open IN, "$out.dupACR.txt" or die $!;
open COM, ">$out.blast.com" or die "";
open BED, ">$out.bedtools.com" or die "";

while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	next unless ( $line[4] > 0 && $line[5] >0);
        $line[1] = $line[0];
	open OUT1, ">./$out.tmp/$line[1]_ACR.bed" or die "";
	open OUT2, ">./$out.tmp/$line[2]_ACR.bed" or die "";
	#open STAT, ">./$out.tmp/$task.ACR.stat" or die "";
	foreach my $pos (sort keys %{$ACR{$line[1]}}){
		my ($chr, $start, $end)= split/,/, $pos;
		my $len=$end-$start;
		#my $seq=substr($seq{$chr}, $start, $len);
		#print OUT1 ">$line[1]_${chr}_${start}_${end};$ACR_dis{$line[1]}{$pos}\n$seq\n";
                #print OUTB "$chr\t$s\t$end\t${geneid}_${id}_${j}_$strand\t1\t$strand\n";
		print OUT1 "$chr\t$start\t$end\t$ACR{$line[1]}{$pos}\t$ACR_dis{$line[1]}{$pos}\n";

	}

	foreach my $pos (sort keys %{$ACR{$line[2]}}){
                my ($chr, $start, $end)= split/,/, $pos;
                my $len=$end-$start;
		print OUT2 "$chr\t$start\t$end\t$ACR{$line[2]}{$pos}\t$ACR_dis{$line[2]}{$pos}\n";
		#my $seq=substr($seq{$chr}, $start, $len);
		#print OUT2 ">$line[2]_${chr}_${start}_${end};$ACR_dis{$line[2]}{$pos}\n$seq\n";
        }
	close OUT1;
	close OUT2;
	print BED "bedtools getfasta -fi $fa -bed ./$out.tmp/$line[1]_ACR.bed -name > ./$out.tmp/$line[1]_ACR.fasta\n";
	print BED "bedtools getfasta -fi $fa -bed ./$out.tmp/$line[2]_ACR.bed -name > ./$out.tmp/$line[2]_ACR.fasta\n";
	print COM "makeblastdb -in ./$out.tmp/$line[1]_ACR.fasta -dbtype nucl -out ./$out.tmp/$line[1]_ACR  -parse_seqids 2>&1>/dev/null ;    blastn -db ./$out.tmp/$line[1]_ACR -query ./$out.tmp/$line[2]_ACR.fasta -outfmt 6 -out ./$out.tmp/$line[1]_$line[2].out6  -task blastn-short -evalue 1e-3 2>&1>/dev/null \n";
}
	close COM;
	close BED;
	`ParaFly -c $out.bedtools.com -CPU 20`;
	`ParaFly -c $out.blast.com -CPU 20`;
	close IN;

open IN, "$out.dupACR.txt" or die $!;
open BED, ">$out.ACR.bed" or die "";
#Store All dup ACR;
my %all;
#store shared ACR
my %share;
#determined one to multi;
my (%dup1, %dup2);
while (<IN>){
        chomp;
        my @line=split/\t/, $_;
	next if $line[4] == 0 && $line[5] == 0;
	if ($line[4] == 0 && $line[5] >0){
		&export_bed( "$line[0]_$line[2]", $line[2], $line[3],\%ACR,\%ACR_dis   );
	}elsif ($line[4] >0 && $line[5] ==0){
                &export_bed( "$line[0]_$line[2]", $line[0], $line[1],\%ACR,\%ACR_dis   );
	}else{
		&definite_status( "$line[0]_$line[2]", $line[1], $line[3], \%ACR, \%ACR_dis );
	}
}
close IN;
close BED;


sub definite_status{

	my ($wgd, $pos1, $pos2, $ACR, $ACR_dis )= @_;
	my ($g1, $g2)= split/_/, $wgd;
	
	my %temp;
	print "./$out.tmp/$wgd.out6\n";
	if (-e "./$out.tmp/$wgd.out6"){
	open INN, "./$out.tmp/$wgd.out6" or die "";
	while (<INN>){
		chomp;
		my @l=split/\t/, $_;
		$l[0]=~s/::.*//g;
		$l[1]=~s/::.*//g;
		$temp{$g2}{$l[0]}{$l[1]}=1;
		$temp{$g1}{$l[1]}{$l[0]}=1;

	}
	close INN;
	}


	my @pos1=split/;/, $pos1;
	for (my $i=0; $i<=$#pos1;$i++){
		my $p=$pos1[$i];
                my $cor=$pos1[$i];
                $cor=~s/,/\t/g;
		my $acr=$$ACR{$g1}{$p};
		if (exists $temp{$g1}{$acr}) {
			my @counter = keys %{$temp{$g1}{$acr}};
			my $counter= join(";", @counter);
			my $multi;
			if (scalar @counter == 1) {
				$multi=1;
			}else{
				$multi = scalar @counter;

			}
			print BED "$cor\t$wgd\t$g1\t$$ACR{$g1}{$p}\t$$ACR_dis{$g1}{$p}\t$counter\t$multi\tshared\n";

		}else{

			print BED "$cor\t$wgd\t$g1\t$$ACR{$g1}{$p}\t$$ACR_dis{$g1}{$p}\t-\t-\tprivate\n";

		}
	}

	my @pos2=split/;/, $pos2;
        for (my $i=0; $i<=$#pos2;$i++){
                my $p=$pos2[$i];
                my $cor=$pos2[$i];
                $cor=~s/,/\t/g;
                my $acr=$$ACR{$g2}{$p};
                if (exists $temp{$g2}{$acr}) {
                        my @counter = keys %{$temp{$g2}{$acr}};
                        my $counter= join(";", @counter);
                        my $multi;
                        if ( scalar @counter == 1) {
                                $multi=1;
                        }else{
                                $multi = scalar @counter;

                        }
                        print BED "$cor\t$wgd\t$g2\t$$ACR{$g2}{$p}\t$$ACR_dis{$g2}{$p}\t$counter\t$multi\tshared\n";

                }else{

                        print BED "$cor\t$wgd\t$g2\t$$ACR{$g2}{$p}\t$$ACR_dis{$g2}{$p}\t-\t-\tprivate\n";

                }
        }

}


sub export_bed{
	my $wgd=shift;
	my $gene=shift;
	my $pos=shift;
	my $ACR=shift;
	my $ACR_dis=shift;
	#my @status=@_;
	
	my @pos=split/;/, $pos;
	
	for (my $i=0; $i<=$#pos;$i++){
		my $p=$pos[$i];
		my $cor=$pos[$i];
		$cor=~s/,/\t/g;
		print BED "$cor\t$wgd\t$gene\t$$ACR{$gene}{$p}\t$$ACR_dis{$gene}{$p}\t-\t-\tprivate\n";
	}
}

#}




