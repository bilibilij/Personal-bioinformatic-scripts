#!/usr/bin/envs perl -w

use strict;

my $usage=<<USAGE;

perl $0 WGD flair_collapsed(filtered).gtf genome.fasta outDir command cpu (re)

3.10685598017575        0.210817717206133       Ovio04457       Ovio13440       0.843856366167796       2-4     0.2-0.3 0.666
2.84468586650301        0.318339100346021       Ovio04459       Ovio13439       0.0359413300890133      2-4     0.3-0.4 0-0.2
3.72187549261204        0.397248113626276       Ovio04461       Ovio13438       0.310746884228657       2-4     0.3-0.4 0.2-0.4

USAGE


die "$usage\n" unless @ARGV ==6 ||@ARGV == 7;

&check_dependency("ParaFly", "ParaFly");
&check_dependency("blastn --help", "blastn");

my ($com, $dir, $cpu)= ($ARGV[4], $ARGV[3], $ARGV[5]);

my (@WGD, %WGD);

open IN, "$ARGV[0]" or die "";
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	push @WGD, "$line[2]_$line[3]";
	$WGD{"$line[2]_$line[3]"}=$_;
}

close IN;


#my %trans;
my %exon;
open IN, "$ARGV[1]" or die "";
while (<IN>){
	chomp;
	next unless m/\texon\t/;
	my @line=split/\t/, $_;
	die "exon match error in $ARGV[1] $$\n" unless m/gene_id "([^"]*)"; transcript_id "([^"]*)"; exon_number "([^"]*)";/;
	my ($gid, $tid, $eid)= ($1, $2, $3);
	$exon{$gid}{$tid}{$eid}="$line[0]\t$line[3]\t$line[4]\t$line[6]";

}
close IN;



my ($seqID , %seq);
open IN, "$ARGV[2]" or die "";
while (<IN>){
	chomp;
	if (m/^>(\S+)/){$seqID=$1;}
	else{$seq{$seqID}.=$_;}
}
close IN;

if ($ARGV[-1] eq "re"){
	`/usr/bin/rm -rf $dir`
}

foreach my $wgd (@WGD){
	my ($g1, $g2)= split/_/, $wgd;
	next unless exists $exon{$g1}&& $exon{$g2};
	&exon_blast_com($wgd, \%exon, \%seq, $com, $dir);
}

`ParaFly -c $com -CPU $cpu`;

sub exon_blast_com{
	my $wgd=shift;
	my $exon=shift;
	my $seq=shift;
	my $com=shift;
	my $dir=shift;
	open COM, ">>$com" or die "";
	my ($g1, $g2)= split/_/, $wgd;
	`mkdir -p $dir/$wgd `;
	#open G1, ">$dir/$wgd/${g1}_exon.fa" or die "";
	#open G2, ">$dir/$wgd/${g2}_exon.fa" or die "";
	&extract_fa($g1, $exon, $seq, "$dir/$wgd/${g1}_exon.fa", "$dir/id_match.txt");
        &extract_fa($g2, $exon, $seq, "$dir/$wgd/${g2}_exon.fa", "$dir/id_match.txt");
        print COM "makeblastdb -in $dir/$wgd/${g1}_exon.fa -dbtype nucl -out $dir/$wgd/${g1}_exon  -parse_seqids 2>&1>/dev/null ;    blastn -db  $dir/$wgd/${g1}_exon -query $dir/$wgd/${g2}_exon.fa -outfmt 6 -out $dir/$wgd/$wgd.exon.blast.tab -evalue 1e-6 2>&1>/dev/null \n";
	close COM;
}

sub extract_fa{
	my $gene=shift;
	my $exon=shift;
	#	my $trans=shift;
	my $seq=shift;
	my $out=shift;
	my $match=shift;
	open MA, ">>$match" or die "";
	open FA, ">$out" or die "";
	my $number=1;
	foreach my $trans (keys %{${$exon}{$gene}}){
                foreach my $ex (keys %{${$exon}{$gene}{$trans}}){
                        my ($chr, $start, $end, $strand) = split/\t/, $$exon{$gene}{$trans}{$ex};
                        my $seq = $$seq{$chr};
				my $len=$end-$start+1;
                        if ($strand eq "+"){
                                $seq=substr($seq,$start-1, $len);
                        }else{
				
                                $seq=substr($seq, $start-1, $len);
                                $seq=~tr/ATCGatcg/TAGCtagc/;
                                $seq=reverse($seq);
                        }
			print FA ">${gene}_trans${number}_${ex}_${strand}_$chr:$start-$end\n$seq\n";
			print MA "${gene}_${trans}_${ex}_${strand}_$chr:$start-$end\t${gene}_trans${number}_${ex}_${strand}_$chr:$start-$end\n";
                }
		$number++;
	}
	close FA;
	close MA;
}



sub check_dependency {
        my $software_info =shift;
        my $match = shift;
        if ($software_info =~ m/$match/) {
            print STDERR "$match:\tOK\n";
        }
        else {
            die "$match:\tFailed\n\n";
        }
}

