#!/usr/bin/envs perl -w

#use strict;

my $usage=<<USAGE;

perl $0 WGD flair_collapsed(filtered).gtf genome.fasta outDir command cpu step (re)

3.10685598017575        0.210817717206133       Ovio04457       Ovio13440       0.843856366167796       2-4     0.2-0.3 0.666
2.84468586650301        0.318339100346021       Ovio04459       Ovio13439       0.0359413300890133      2-4     0.3-0.4 0-0.2
3.72187549261204        0.397248113626276       Ovio04461       Ovio13438       0.310746884228657       2-4     0.3-0.4 0.2-0.4

USAGE


die "$usage\n" unless @ARGV ==7 ||@ARGV == 8;

&check_dependency("ParaFly 2>&1", "ParaFly");
&check_dependency("blastn --help 2>&1", "blastn");
&check_dependency("bedtools getfasta --help 2>&1", "getfasta");



my ($com, $dir, $cpu, $step)= ($ARGV[4], $ARGV[3], $ARGV[5], $ARGV[6]);

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
my %gene_len;
open IN, "$ARGV[1]" or die "";
while (<IN>){
	chomp;
	if (m/\texon\t/){
	my @line=split/\t/, $_;
	die "exon match error in $ARGV[1] $$\n" unless m/gene_id "([^"]*)"; transcript_id "([^"]*)"; exon_number "([^"]*)";/;
	my ($gid, $tid, $eid)= ($1, $2, $3);
	$exon{$gid}{$tid}{$eid}="$line[0]\t$line[3]\t$line[4]\t$line[6]";
	}elsif (m/\ttranscript\t/){
		my @line=split/\t/, $_;
		die "exon match error in $ARGV[1] $$\n" unless m/gene_id "([^"]*)"; transcript_id "([^"]*)";/;
		my $gid= $1;
		if (exists $gene_len{$gid}){
			my @info = split/\t/, $gene_len{$gid};
			if ($info[1] >= $line[3]){ $info[1] = $line[3];  }
			if ($info[2] <= $line[4]){ $info[2] = $line[4];  }
			$gene_len{$gid}= join("\t", @info);
		}else{
			$gene_len{$gid}="$line[0]\t$line[3]\t$line[4]\t$line[6]";
		}
	}

}
close IN;

my %exon_organ=%exon;
open IN, "/usr_storage/jcf/2.zhugecai/00.diploid/Ovio_geta_best.gff3" or die "";
while (<IN>){
	chomp;
	if (/\texon\t/){
        my @line=split/\t/, $_;
	die "exon match error in Ovio_geta_best.gff3 $$\n" unless m/ID=([^\.]*)\..*exon(\d+);Parent=([^;]*);/;
	my ($gid, $eid, $tid)= ($1, $2, $3);
	$tid.="_origin";
	$exon_organ{$gid}{$tid}{$eid}="$line[0]\t$line[3]\t$line[4]\t$line[6]";
	}elsif (m/\tgene\t/){
		my @line=split/\t/, $_;
		die "gene match error in Ovio_geta_best.gff3 $$\n" unless m/ID=([^;]*);Name/;
		my $gid=$1;
		if (exists $gene_len{$gid}){
                        my @info = split/\t/, $gene_len{$gid};
                        if ($info[1] >= $line[3]){ $info[1] = $line[3];  }
                        if ($info[2] <= $line[4]){ $info[2] = $line[4];  }
                        $gene_len{$gid}= join("\t", @info);
                }else{
                        $gene_len{$gid}="$line[0]\t$line[3]\t$line[4]\t$line[6]";
                }
	}
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


mkdir "$dir/exon_blast/" unless -e "$dir/exon_blast/";



if ($step ==1 || $step eq "all"){
open COM, ">$com" or die "";
open MA, ">$dir/id_match.txt" or die "";

foreach my $wgd (@WGD){
	my ($g1, $g2)= split/_/, $wgd;
	next unless exists $exon{$g1}&& $exon{$g2};
	#&exon_blast_com($wgd, \%exon, \%seq, $com, $dir);
	&exon_blast_com($wgd, \%exon, \%seq, $dir);
}
close COM;
close MA;

`ParaFly -c $com -CPU $cpu`;
}


if ($step ==2 || $step eq "all"){
open COM, ">$com" or die "";
open MA, ">$dir/id_match_specific_exon.txt" or die "";
`/usr/bin/rm -rf $dir/get_specific_exon_fa.com` if -e "$dir/get_specific_exon_fa.com" ;
`/usr/bin/rm -rf $dir/get_specific_exon_fa.com.completed` if -e "$dir/get_specific_exon_fa.com.completed" ;

foreach my $wgd (@WGD){
        my ($g1, $g2)= split/_/, $wgd;
	&exon_blast_com($wgd, \%exon_organ, \%seq, $dir, \%gene_len);
}
close COM;
close MA;

`ParaFly -c $dir/get_specific_exon_fa.com -CPU $cpu `;
`ParaFly -c $com -CPU $cpu`;

}



###Step2: exon number check
my (%idf, %idr);
#idf: key is Nanopore ID and values is newly renamed transcript ID;
#idr: reverse the key and value in the %idf;
open IN, "$dir/id_match.txt" or die "";
while (<IN>){
	chomp;
	my @line=split/\t/, $_;
	my @old=split/_/, $line[0];
	my @new=split/_/, $line[1];
	$idf{$old[0]}{$old[1]}=$new[1];
	$idr{$old[0]}{$new[1]}=$old[1];
}
close IN;

my %exon_num;
foreach my $wgd (@WGD){
        my ($g1, $g2)= split/_/, $wgd;
        next unless exists $exon{$g1}&& $exon{$g2};
	foreach my $trans (keys %{$exon{$g1}}){
		my $transid=$idf{$g1}{$trans};
		my $exon_number=scalar keys %{$exon{$g1}{$trans}};
		$exon_num{$wgd}{$exon_number}{$g1}{$transid}="${g1}_$trans";
	}

	foreach my $trans (keys %{$exon{$g2}}){
                my $transid=$idf{$g2}{$trans};
                my $exon_number=scalar (keys %{$exon{$g2}{$trans}});
                $exon_num{$wgd}{$exon_number}{$g2}{$transid}="${g2}_$trans";
        }
}

mkdir "$dir/results/" unless -e "$dir/results/" ;
mkdir "$dir/exon_flank/" unless -e " $dir/exon_flank";

open NUM, ">$dir/results/exon_num.txt" or die "";
my %share_trans;
my %share_wgd;
my %trans_exon;
#my %final_shared;
my (%g1_shared, %g2_shared);
foreach my $wgd (sort keys %exon_num){
	foreach my $exon_num (sort {$b<=>$a} keys %{$exon_num{$wgd}}){
		foreach my $geneid (keys %{$exon_num{$wgd}{$exon_num}}){
			foreach my $transid (sort keys %{$exon_num{$wgd}{$exon_num}{$geneid}}){
				if (scalar %{$exon_num{$wgd}{$exon_num}} == 1 ){
					print NUM "$wgd\t$exon_num\t$geneid\t$transid\t$exon_num{$wgd}{$exon_num}{$geneid}{$transid}\tprivate\n";
				}else{
					#&exon_blast_com($wgd, \%exon, \%seq, $dir, \%idf);
					print NUM "$wgd\t$exon_num\t$geneid\t$transid\t$exon_num{$wgd}{$exon_num}{$geneid}{$transid}\tshared\n";
					#@$share_trans{$geneid->$transid};
					#push @{$share_trans -> {$geneid}}, $transid;
					my $old_id=(split/_/, $exon_num{$wgd}{$exon_num}{$geneid}{$transid})[1];
					$share_trans{$geneid}{$old_id}=1;
					$share_wgd{$wgd}=1;
					$trans_exon{$wgd}{$exon_num}{$geneid}{$old_id}=1;

					my ($g1, $g2)=split/_/, $wgd;
					if ($geneid eq $g1){$g1_shared{$wgd}{$exon_num}{$old_id}=1;   }
					elsif ($geneid eq $g2){$g2_shared{$wgd}{$exon_num}{$old_id}=1; }
					else{die "do not match g1 and g2 $$\n";}

				}
			}
		}
	}
}
close NUM;


#open COM, ">$dir/exon_flanking_blast.com" or die "";
###Step3 : exon blast check
#my %blast_check;
#foreach my $wgd (sort keys %share_wgd){
	#&exon_blast_com($wgd, \%exon, \%seq, $dir, \%idf, \%share_trans);
#}
#close COM;



if ($step == 3 || $step eq "all"){
mkdir "$dir/1.blast_flank/"  unless -e "$dir/1.blast_flank/";
#`/usr/bin/rm -rf $dir/exon_sort/*`;
open OUT, ">$dir/exon_sort.com" or die "";
open FLA, ">$dir/exon_flank.com" or die $$;
open BED, ">$dir/get_fasta_new.com" or die $$;
foreach my $wgd (sort keys %trans_exon){
	foreach my $exon_num ( sort {$b<=>$a} keys  %{$trans_exon{$wgd}}){
		foreach my $geneid (keys %{$trans_exon{$wgd}{$exon_num}}){
			#`/usr/bin/rm -rf $dir/exon_sort/${wgd}_${exon_num}_${geneid}.$i.fa` if -e "$dir/exon_sort/${wgd}_${exon_num}_${geneid}.$i.fa";
			foreach my $trans ( keys %{$trans_exon{$wgd}{$exon_num}{$geneid}}){
				my $ss = (split/\t/, $exon{$geneid}{$trans}{"0"})[3];
				print "$trans\t$ss\n";
					for (my $i=0;$i<=$exon_num-1;$i++){
						my $j;
						if ($ss eq "-"){$j=$exon_num-1-$i;}
						else{$j=$i;}
						my ($chr, $start, $end, $strand) = split/\t/, $exon{$geneid}{$trans}{$j};
						##my $seq = $seq{$chr};
						my $len=$end-$start+1;
						my $s=$start-1;
						my $id=$idf{$geneid}{$trans};
						#open OUTB, ">$dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$i.bed" or die "";
						#print OUTB "$chr\t$s\t$end\t${geneid}_${id}_${j}_$strand\t1\t$strand\n";
						#close OUTB;
						#print BED "bedtools getfasta -fi $ARGV[2] -bed $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$i.bed -s -name |sed 's/::.*//g'  >> $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$i.fa \n";

						if ($ss eq "+"){
						$s=$start-1-45;
						my $e=$start+45;
                                                open OUTB, ">$dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$i.flank.up.bed" or die "";
                                                print OUTB "$chr\t$s\t$e\t${geneid}_${id}_${j}_$strand\t1\t$strand\n";
                                                close OUTB;

						print BED  "bedtools getfasta -fi $ARGV[2] -bed $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$j.flank.up.bed -s -name |sed 's/::.*//g'  >> $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$i.flank.up.fa \n";
						$s=$end-1-45;
						$e=$end+45;

                                                open OUTB, ">$dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$i.flank.down.bed" or die "";
                                                print OUTB "$chr\t$s\t$e\t${geneid}_${id}_${j}_$strand\t1\t$strand\n";
                                                close OUTB;

						print BED  "bedtools getfasta -fi $ARGV[2] -bed $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$j.flank.down.bed -s -name |sed 's/::.*//g'  >> $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$i.flank.down.fa \n";
						}else{
                                                my $e=$end+45;
						$s=$end-1-45;
                                                open OUTB, ">$dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$i.flank.up.bed" or die "";
                                                print OUTB "$chr\t$s\t$e\t${geneid}_${id}_${j}_$strand\t1\t$strand\n";
                                                close OUTB;

                                                print BED  "bedtools getfasta -fi $ARGV[2] -bed $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$j.flank.up.bed -s -name |sed 's/::.*//g'  >> $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$i.flank.up.fa \n";
                                                $s=$start-1-45;
                                                $e=$start+45;

                                                open OUTB, ">$dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$i.flank.down.bed" or die "";
                                                print OUTB "$chr\t$s\t$e\t${geneid}_${id}_${j}_$strand\t1\t$strand\n";
                                                close OUTB;

                                                print BED  "bedtools getfasta -fi $ARGV[2] -bed $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$id.$j.flank.down.bed -s -name |sed 's/::.*//g'  >> $dir/1.blast_flank/${wgd}_${exon_num}_${geneid}.$i.flank.down.fa \n";

						}
				}
			}

		}
	my ($g1, $g2)= split/_/, $wgd;
	for (my $i=0;$i<=$exon_num-1;$i++){
		print OUT "makeblastdb -in $dir/1.blast_flank/${wgd}_${exon_num}_${g1}.$i.fa  -dbtype nucl -out $dir/1.blast_flank/${wgd}_${exon_num}_${g1}.$i  2>&1>/dev/null ; blastn -db $dir/1.blast_flank//${wgd}_${exon_num}_${g1}.$i  -query $dir/1.blast_flank/${wgd}_${exon_num}_${g2}.$i.fa -outfmt 6 -out $dir/1.blast_flank/${wgd}_${exon_num}.$i.blast.tab -evalue 1e-3 -task blastn-short 2>&1>/dev/null\n";
		print FLA "makeblastdb -in $dir/1.blast_flank/${wgd}_${exon_num}_${g1}.$i.flank.down.fa  -dbtype nucl -out $dir/1.blast_flank/${wgd}_${exon_num}_${g1}.$i.flank.down  2>&1>/dev/null ;blastn -db $dir/1.blast_flank/${wgd}_${exon_num}_${g1}.$i.flank.down  -query $dir/1.blast_flank/${wgd}_${exon_num}_${g2}.$i.flank.down.fa -outfmt 6 -out $dir/1.blast_flank/${wgd}_${exon_num}.${i}_flank_down.tab -evalue 1e-3 -task blastn-short 2>&1>/dev/null\n";
		print FLA "makeblastdb -in $dir/1.blast_flank/${wgd}_${exon_num}_${g1}.$i.flank.up.fa  -dbtype nucl -out $dir/1.blast_flank/${wgd}_${exon_num}_${g1}.$i.flank.up  2>&1>/dev/null ;blastn -db $dir/1.blast_flank/${wgd}_${exon_num}_${g1}.$i.flank.up  -query $dir/1.blast_flank/${wgd}_${exon_num}_${g2}.$i.flank.up.fa -outfmt 6 -out $dir/1.blast_flank/${wgd}_${exon_num}.${i}_flank_up.tab -evalue 1e-3 -task blastn-short 2>&1>/dev/null\n";
	}
	}
}
close OUT;
close FLA;
close BED;
`ParaFly -c $dir/get_fasta_new.com -CPU $cpu`;
#`ParaFly -c $dir/exon_sort.com -CPU $cpu`;
`ParaFly -c $dir/exon_flank.com -CPU $cpu`;
}

#my %exon_shared;

if ($step == 4 || $step eq "all"){
my (%all, %private, %shared);
foreach my $wgd (sort keys %g1_shared)
{
	
	foreach my $exon_num (sort keys %{$g1_shared{$wgd}})
	{
		foreach my $g1_transid (sort keys %{$g1_shared{$wgd}{$exon_num}})
		{
			foreach my $g2_transid (sort keys %{$g2_shared{$wgd}{$exon_num}})
			{
				$all{$wgd}{$exon_num}{$g1_transid}{$g2_transid}=0;
			}
		
		}

	}
}


%shared=%all;


foreach my $wgd (sort keys %all)
{
	foreach my $exon_num (sort keys %{$all{$wgd}})
	{
		#g1 = rexon ; g2 = qexon;
		my ($temp_shared, %temp_private, %rexon_len, %qexon_len);
		$temp_shared = \%{$all{$wgd}{$exon_num}};		

		#my $id=$idf{$geneid}{$trans};
		my ($g1, $g2)= split/_/, $wgd;
	        for (my $i=0;$i<=$exon_num-1;$i++)
		{
                    #Step1.1: get each exon length of each transcritps;
                    &fa_len("$dir/1.blast_flank/${wgd}_${exon_num}_${g2}.$i.fa", \%qexon_len);
                    &fa_len("$dir/1.blast_flank/${wgd}_${exon_num}_${g1}.$i.fa", \%rexon_len);

                    #Step2.2: extract align result from blast.tab
                    #&blast_extract("$dir/exon_sort/${wgd}_${exon_num}.$i.blast.tab", \%rblast, \%rblast_len);
		    &func( "$dir/1.blast_flank/${wgd}_${exon_num}.$i.blast.tab",$wgd,  \%rexon_len, \%qexon_len, $temp_shared, \%temp_private, \%idf); #private for add non shared information; shared for calculate score;

                }

                #deal with private
                foreach my $trans_g1 (sort keys %temp_private)
                {
                        foreach my $trans_g2 (sort keys %{$temp_private{$trans_g1}})
                        {
                                $private{$wgd}{$exon_num}{$trans_g1}{$trans_g2}=1;
				#delete $shared{$wgd}{$exon_num}{$trans_g1}{$trans_g2};
                        }
                }

		#adding shared score
		foreach my $trans_g1 (sort keys %$temp_shared)
		{
			foreach my $trans_g2 (sort keys %{$$temp_shared{$trans_g1}})
			{
				#next if exists $private{$wgd}{$exon_num}{$trans_g1}{$trans_g2};
				 if (exists $private{$wgd}{$exon_num}{$trans_g1}{$trans_g2}) {my $num=$$temp_shared{$trans_g1}{$trans_g2}/$exon_num; print "$wgd\t$exon_num\t$idf{$g1}{$trans_g1}\t$idf{$g2}{$trans_g2}\t$$temp_shared{$trans_g1}{$trans_g2}\t$num\tprivate\n";delete $shared{$wgd}{$exon_num}{$trans_g1}{$trans_g2};}
				 else{
				$shared{$wgd}{$exon_num}{$trans_g1}{$trans_g2}=$$temp_shared{$trans_g1}{$trans_g2};}
			}
		}
	}
}

my %sort;

foreach my $wgd (sort keys %shared){
	my ($g1, $g2) =split/_/, $wgd;
	foreach my $exon_num (sort keys %{$shared{$wgd}}){
		foreach my $trans1 (sort keys %{$shared{$wgd}{$exon_num}}){
			foreach my $trans2  (sort keys %{$shared{$wgd}{$exon_num}{$trans1}}){
				#print "$wgd\t$exon_num\t$trans1\t$trans2\t$shared{$wgd}{$exon_num}{$trans1}{$trans2}\n";
				my $num=$shared{$wgd}{$exon_num}{$trans1}{$trans2}/$exon_num;
				#print "$wgd\t$exon_num\t$idf{$g1}{$trans1}\t$idf{$g2}{$trans2}\t$shared{$wgd}{$exon_num}{$trans1}{$trans2}\t$num\n";
				$sort{$wgd}{$exon_num}{$num}="$trans1\t$trans2";
			}
		}
	}
}

open OUT, ">$dir/results/sorted_shared.txt" or die "";


my %final;
foreach my $wgd (sort keys %sort){
	my ($g1, $g2) =split/_/, $wgd;
	foreach my $exon_num (sort keys %{$sort{$wgd}}){
		my (%temp1, %temp2);
		foreach my $score (sort {$b<=>$a} keys %{$sort{$wgd}{$exon_num}}){
			my ($trans1, $trans2) = split/\t/, $sort{$wgd}{$exon_num}{$score};
			next if exists $temp1{$trans1} || exists $temp2{$trans2};
			$temp1{$trans1}=$trans2;
			$temp2{$trans2}=$trans1;
		}
		foreach my $trans1 (sort keys %temp1){
			my $trans2=$temp1{$trans1};
			print OUT  "$wgd\t$exon_num\t$idf{$g1}{$trans1}\t$trans1\t$idf{$g2}{$trans2}\t$trans2\tshared\n";
		}
	}
}

close OUT;

}
#Ovio17090_trans2_6

if ($step == 5 || $step eq "all")
{
	my %hash;
	open IN, "$dir/results/sorted_shared.txt" or die "";
	while (<IN>){
		chomp;
		my @line=split/\t/, $_;
		my ($g1, $g2) = split/_/, $line[0];
		my ($g1t, $g2t)= ("${g1}_$line[2]", "${g2}_$line[4]");
		$hash{$g1t}=1;
		$hash{$g2t}=1;
	}
	close IN;
	open OUT, ">$dir/results/private.txt" or die "";
	open IN, "$dir/results/exon_num.txt" or die "";
	while (<IN>){
		chomp;
		my @line=split/\t/, $_;
		pop @line;
		if (exists $hash{$line[4]}){ push @line, "shared"; }
		else{ push @line ,"private";    }
		my $out=join("\t", @line); 
		print OUT "$out\n";
	}
	close IN;
}



sub func
{
	my $blast=shift;
	my $wgd=shift;
	#my $down=shift;
	#my $up=shift;
	my $rlen=shift;
	my $qlen=shift;
	my $temp_shared=shift;
	my $temp_private=shift;
	my $idf=shift;
	die "can not match exon$$\n" unless $blast =~ m/[^\.]*.(\d+).blast.tab/;
	my $exon=$1;
	my (%gene, %perfect);
	&blast_extract($blast, \%gene, \%perfect, $temp_private);

	#my (%up, %up_perfect);
	#&blast_extract($up, \%up, \%up_perfect, $temp_private);

	#my (%down, %down_perfect);
	#&blast_extract($down, \%down, \%down_perfect, $temp_private);

	#print "$down\n";
	#print %down;
	#print %down_perfect;
	#
	my ($gene1, $gene2) =split/_/, $wgd;
	foreach my $g1 ( keys %{$temp_shared})
	{
		foreach my $g2 (keys %{$$temp_shared{$g1}})
		{
			#my ($g1t, $g2t)= ("${gene1}_$$idf{$gene1}{$g1}_$exon", "${gene2}_$$idf{$gene2}{$g2}_$exon");
			my ($g1t, $g2t)= ("${gene1}_$$idf{$gene1}{$g1}", "${gene2}_$$idf{$gene2}{$g2}");
			my ($g1_len,$g2_len)=($$rlen{$g1t}, $$qlen{$g2t} );
			#print "$gene{$g1t}{$g2t}\t$perfect{$g1t}{$g2t}/$g1_len/$g2_len\n";
			if (not exists $gene{$g1t}{$g2t}){$$temp_private{$g1}{$g2}=1;}
			elsif( $perfect{$g1t}{$g2t}/$g1_len <0.5 || $perfect{$g1t}{$g2t}/$g2_len<0.5) {$$temp_private{$g1}{$g2}=1;}
			else{
				#my $coverage=($perfect{$g1t}{$g2t}/$g1_len+$perfect{$g1t}{$g2t}/$g2_len)/2;
				#print "$wgd\t$g1\t$g2\t$coverage\n";
				#print "$blast\t$wgd\t$g1t\t$g2t\t$g1_len\t$g2_len\t";
				#$perfect{$g1t}{$g2t}\t$perfect{$g1t}{$g2t}\n";
				#print %perfect;
				#print "\n";
				$$temp_shared{$g1}{$g2}+=($perfect{$g1t}{$g2t}/$g1_len+$perfect{$g1t}{$g2t}/$g2_len)/2;
			}

		}

	}
}



sub blast_extract{
	my $file=shift;
	my $mat=shift;
	my $perfect=shift;
	my $temp_private=shift;
	open IN, "$file" or die $$;
	while (<IN>){
		chomp;
		my ($qid, $rid, $iden, $aln, $mis, $gap, $qstart, $qend, $rstart, $rend, $evalue, $score)=split/\t/, $_;
		my $qlen=$qend-$qstart+1;
		my $rlen=$rend-$rstart+1;
		my $match=$aln-$mis-$gap;
		die "match err $$\n" unless $qid=~m/^([^_]*)_([^_]*)_(\d+)_.$/;
		$qid="$1_$2";
		die "match err $$\n" unless $rid=~m/^([^_]*)_([^_]*)_(\d+)_.$/;
                $rid="$1_$2";
		#next if exists $$temp_private{$rid};
                $$mat{$rid}{$qid}=1;
                $$perfect{$rid}{$qid}+=$match;

	}
	close IN;
}


sub fa_len{
	my $file=shift;
	my $trans_len=shift;
	my ($geneid, $transid, $exonID);
	open IN, "$file" or die "";
	while (<IN>){
		chomp;
		if (m/^>([^_]*)_([^_]*)_(\d+)_.$/){
		#if (m/^>(\S+)_.$/){
			$seqID="$1_$2";
		}else{
			#$$trans_len{$geneid}{$transid}{$exonID}+=length($_);
			$$trans_len{$seqID}=length($_);
		}
	}
	close IN;
}





sub exon_blast_com{
	 my ($wgd, $exon, $seq, $dir, $idf, $share_trans, $gene, $gene_len);
        if (scalar @_ == 5){
                ($wgd, $exon, $seq, $dir, $gene_len) = @_;

        }elsif (scalar @_ == 6){
                ($wgd, $exon, $seq, $dir, $idf, $share_trans) = @_;

        }

	my ($g1, $g2)= split/_/, $wgd;
	if (scalar @_ == 5){
		if ($step == 1 ||$step == 2|| $step eq "all"){
		&bedtools_fa($$gene_len{$g1}, $g1, "$dir/exon_blast/${g1}_gene.bed", "$dir/get_specific_exon_fa.com", "$dir/exon_blast/${g1}_gene.fa");
                &bedtools_fa($$gene_len{$g2}, $g2, "$dir/exon_blast/${g2}_gene.bed", "$dir/get_specific_exon_fa.com", "$dir/exon_blast/${g2}_gene.fa");
		open CBED, ">>$dir/get_specific_exon_fa.com" or die "";
		print CBED "sort -k1,1n -k2,2n -k3,3n $dir/exon_blast/${g1}_gene.bed |uniq | bedtools getfasta -fi $ARGV[2] -bed - -s -name >>$dir/exon_blast/${g1}_gene.fa\n";
                print CBED "sort -k1,1n -k2,2n -k3,3n $dir/exon_blast/${g2}_gene.bed |uniq | bedtools getfasta -fi $ARGV[2] -bed - -s -name >>$dir/exon_blast/${g2}_gene.fa\n";
		close CBED;
		&extract_fa($g1, $exon, $seq, "$dir/exon_blast/${g1}_exon.fa");
		&extract_fa($g2, $exon, $seq, "$dir/exon_blast/${g2}_exon.fa");}
		print COM "makeblastdb -in $dir/exon_blast/${g1}_exon.fa -dbtype nucl -out $dir/exon_blast/${g1}_exon  -parse_seqids 2>&1>/dev/null ;    blastn -db  $dir/exon_blast/${g1}_exon -query $dir/exon_blast/${g2}_gene.fa -outfmt 6 -out $dir/exon_blast/$wgd.${g2}.exon.blast.tab -evalue 1e-3 -task blastn-short  2>&1>/dev/null \n";
		print COM "makeblastdb -in $dir/exon_blast/${g2}_exon.fa -dbtype nucl -out $dir/exon_blast/${g2}_exon  -parse_seqids 2>&1>/dev/null ;    blastn -db  $dir/exon_blast/${g2}_exon -query $dir/exon_blast/${g1}_gene.fa -outfmt 6 -out $dir/exon_blast/$wgd.${g1}.exon.blast.tab -evalue 1e-3 -task blastn-short  2>&1>/dev/null \n";
	}
	if (scalar @_ == 6){
		if ($step == 3 || $step eq "all"){
                &extract_fa($g1, $exon, $seq, "$dir/exon_flank/${g1}_exon_flank.fa", $idf, $share_trans);
                &extract_fa($g2, $exon, $seq, "$dir/exon_flank/${g2}_exon.flank.fa", $idf, $share_trans);}
		print COM "makeblastdb -in $dir/exon_flank/${g1}_exon_flank.fa -dbtype nucl -out $dir/exon_flank/${g1}_exon_flank  -parse_seqids 2>&1>/dev/null ;    blastn -db  $dir/exon_flank/${g1}_exon_flank -query $dir/exon_flank/${g2}_exon.flank.fa -outfmt 6 -out $dir/exon_flank/${wgd}_exon.flank.tab -evalue 1e-3 2>&1>/dev/null \n";
	}
}	

sub bedtools_fa{
	#my $chr=shift;
	#my $start=shift;
	#my $end=shift;
	#my $strand=shift;
	my $cor=shift;
	my $name=shift;
	my $bed=shift;
	#my $com=shift;
	my $out=shift;
	my ($chr,$start,$end,$strand)=split/\t/, $cor;
	#open CBED, ">>$com" or die "";
	open BED, ">>$bed" or die "";
	print BED "$chr\t$start\t$end\t$name\t1\t$strand\n";
	close BED;
	#print CBED "sort -k1,1n -k2,2n -k3,3n $bed |uniq | bedtools getfasta -fi $ARGV[2] -bed - -s -name >>$out\n";
	#close CBED;

}

sub extract_fa{
	my ($gene, $exon, $seq, $out, $idf, $share_trans);
	if (scalar @_ == 4){
		($gene, $exon, $seq, $out) = @_;
	
	}elsif (scalar @_ == 6){
		($gene, $exon, $seq, $out, $idf, $share_trans) = @_;
	}
	#my $match=shift;
	#open MA, ">>$match" or die "";Ovio04499_Ovio13395.exon.blast.tab
	if (scalar @_ == 4){
	open FA, ">$out" or die "";
	open OO, ">$out.bed" or die "";
	close OO;
	my $number=1;
	foreach my $trans (sort keys %{${$exon}{$gene}}){
                foreach my $ex (sort keys %{${$exon}{$gene}{$trans}}){
			&bedtools_fa( $$exon{$gene}{$trans}{$ex}, "${gene}_trans${number}_${ex}", "$out.bed", $out );
			my ($chr, $start, $end, $strand) = split/\t/, $$exon{$gene}{$trans}{$ex};
			#my $seq = $$seq{$chr};
			#	my $len=$end-$start+1;
			#if ($strand eq "+"){
			#        $seq=substr($seq,$start-1, $len);
			#}else{
				
			#        $seq=substr($seq, $start-1, $len);
			#        $seq=~tr/ATCGatcg/TAGCtagc/;
			#        $seq=reverse($seq);
			#}
			#print FA ">${gene}_trans${number}_${ex}_${strand}_$chr:$start-$end\n$seq\n";
			print MA "${gene}_${trans}_${ex}_${strand}_$chr:$start-$end\t${gene}_trans${number}_${ex}_${strand}_$chr:$start-$end\n";
                }

		$number++;
	}
	open CBED, ">>$dir/get_specific_exon_fa.com" or die "";
	print CBED "sort -k1,1n -k2,2n -k3,3n $out.bed |uniq | bedtools getfasta -fi $ARGV[2] -bed - -s -name >> $out\n";
	close CBED;
	close FA;
	}elsif (scalar @_ == 6){
	open FA, ">$out" or die "";
	foreach my $trans (sort keys %{${$exon}{$gene}}){
		print "$trans\t$$share_trans{$gene}{$trans}\n";
		#next unless $trans eq $trans_homo;
		next unless exists $$share_trans{$gene}{$trans};
		print FA ">${gene}_$$idf{$gene}{$trans}\n";
		#print FA ">${gene}_$trans\n";
		foreach my $ex (sort keys %{${$exon}{$gene}{$trans}}){
			my ($chr, $start, $end, $strand) = split/\t/, $$exon{$gene}{$trans}{$ex};
			my $Chr_seq = $$seq{$chr};
			my $seq;
			#extract fasta
			if ($strand eq "+"){
				#$seq=substr($$seq{$chr},$start-1-45, 90);
				#$seq.=substr($$seq{$chr},$end-1-45, 90);
				$seq=substr($Chr_seq,$start-1-45, 90);
                                $seq.=substr($Chr_seq,$end-1-45, 90);

			}else{
				$seq=substr($Chr_seq,$start-1-45, 90);
				$seq.=substr($Chr_seq,$end-1-45, 90);
                                $seq=~tr/ATCGatcg/TAGCtagc/;
                                $seq=reverse($seq);
                        }
			print FA "$seq";
		}
		print FA"\n";
	}
	close FA;
	}
}

sub check_dependency {
        my $com =shift;
        my $match = shift;
	my $software_info=`$com`;
        if ($software_info =~ m/$match/) {
            print STDERR "$match:\tOK\n";
        }
        else {
            die "$match:\tFailed\n\n";
        }
}

