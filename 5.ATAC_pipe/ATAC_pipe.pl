#!/usr/bin/envs perl 

use strict;
use lib("/home/jcf/scripts/jcf_perlib/");
use lsyPig;
use threads;
use Getopt::Long;

my $usage=<<USAGE;
usage:

This personal script contain four step
1.ATAC.index
2.fastp
3.Bowtie
4.Genrich

perl $0

--genome genome.fasta  --fq fq.list --out outprefix --cpu cpu --size chrom_size --tss tss.bed

fq.list>>>>
/usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_F.clean.paired.R1.fq.gz /usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_F.clean.paired.R2.fq.gz
/usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_L.clean.paired.R1.fq.gz /usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_L.clean.paired.R2.fq.gz
/usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_R.clean.paired.R1.fq.gz /usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_R.clean.paired.R2.fq.gz
/usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_S.clean.paired.R1.fq.gz /usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_S.clean.paired.R2.fq.gz
/usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_Seed.clean.paired.R1.fq.gz /usr_storage/data/zhugecai_ATAC/clean/OV_ATAC_Seed.clean.paired.R2.fq.gz
>>>>>

Alternative parameters : --step (default :all) ,you can choose which step to run ;
such as if you want to run step 4.genrich  without rerun step 1 to 3, you could use \"--step 4 \" to specify;
USAGE


#die "$usage\n" if @ARGV == 0;

my ($genome, $fq, $outprefix, $cpu, $size, $tss_bed);
my $step="all";

GetOptions(
    'genome:s' => \$genome,
    'fq:s' => \$fq,
    'out:s' => \$outprefix,
    'cpu:s' => \$cpu,
    'step:s' => \$step,
    'tss:s' => \$tss_bed,
    'size:s' => \$size,
) or die $!;

die "$usage\n" if ($genome eq "" or $fq eq "" or $outprefix eq "" );

#check dependencies
print STDERR "\n============================================\n";
print STDERR "Detecting the dependency softwares:\n";
##mafft
#my $software_info = `mafft --help 2>&1 `;
check_dependency( "macs2 2>&1" , "macs2"  );
check_dependency( "fastp --help 2>&1", "fastp");
check_dependency( "samtools", "samtools"   );
check_dependency( "picard ", "icard"   );

print STDERR "Dependencies:\tOK\n";
print STDERR "====================================\n\n\n";



chomp(my $pwd=`pwd`);
$pwd="$pwd/$outprefix/";

my (%fq_list,%fq_fastp);
open IN, "$fq" or die "";
while (<IN>){
	chomp;
	my @line=split/\s+/, $_;
	#foreach my 
	die "fq nor pired in fastq.list \n" unless @line == 2;
	chomp(my $p=`basename $line[0]`);
	#.clean.paired.R1.fq.gz
	$p =~ s/.clean.paired.R..fq.gz$//g;
	$fq_list{$p}=$_;
	$fq_fastp{$p}="$pwd/2.fastp/$p.fastp.R1.fq.gz\t$pwd/2.fastp/$p.fastp.R2.fq.gz";
}
close IN;



#`mkdir $outprefix/1.ATAC.index` unless -e "$outprefix/1.ATAC.index" && -d "$outprefix/1.ATAC.index";

my ($seqID, %seq);
open IN, "$genome" or die "";
while (<IN>){
        if (m/^>(\S+)/){$seqID = $1;}
        else{$seq{$seqID}.=$_;}
}
close IN;



`mkdir $outprefix/1.ATAC.index` unless -e "$outprefix/1.ATAC.index" && -d "$outprefix/1.ATAC.index";
chomp(my $a=`basename $genome`) ;
open OUT, ">$outprefix/1.ATAC.index/$a" or die "";
my $chr_len;
foreach my $id (keys %seq){
        next unless $id =~m/^Chr/;
	next if $id =~ m/ChrM/;
	next if $id =~ m/ChrC/;
        print OUT ">$id\n$seq{$id}\n";
        $chr_len += length($seq{$id}) ;
}


chdir "$outprefix/";
chomp($genome=`basename $genome`);

#if ($step eq "all")
my @ok= ( "1.ATAC.index.ok", "2.fastp.ok", "3.Bowtie.ok", "4.Downstream.ok" );
unless ($step eq "all"){
        if ( -d "pre_conf" & -e "pre_conf" ) {`/usr/bin/rm -rf pre_conf`;}
        else {  `mkdir pre_conf && mv *.ok pre_conf/ `; }
        for (my $i=1;$i<=4;$i++){
                die "\$step : $step not numeric\n" unless $step =~ m/^\d+$/;
                if ($step == $i){ my $dir=$ok[$i-1];  $dir=~ s/.ok//g; `/usr/bin/rm -rf $dir`; }
                else{open OUT, ">$ok[$i-1]" or die "";close OUT;}
        }
}


#Step 1 ATAC index
print STDERR "\n============================================\n";
print STDERR "Step 1: ATAC genome index  " . "(" . (localtime) . ")" . "\n";
my $index_run;
unless (-e "1.ATAC.index.ok"){
        `mkdir 1.ATAC.index` unless -e "1.ATAC.index" && -d "1.ATAC.index";
        $index_run = threads->create(\&com, "bowtie2-build -q  1.ATAC.index/$genome  1.ATAC.index/index "  );

        #command ( "bowtie2-build -q  1.ATAC.index/$genome  1.ATAC.index/index "  );
        #open OUT ,">1.ATAC.index/index.ok" or die $!;  close OUT;
}else{
        print STDERR "Skip Step 0 for the file 1.ATAC.index/index.ok exists \n";
}


#Step 2 fastp
my $fastp_run;
print STDERR "\n============================================\n";
print STDERR "Step 2: fastp  " . "(" . (localtime) . ")" . "\n";
unless ( -e "2.fastp.ok") {
        `mkdir 2.fastp` unless -e "2.fastp" && -d "2.fastp";

        open OUT, ">2.fastp/command" or die "";
        foreach my $p (sort keys %fq_list){
                my ($f1, $f2) = split/\s+/, $fq_list{$p};
                print OUT "  fastp  -i $f1  -I $f2  -o 2.fastp/$p.fastp.R1.fq.gz -O 2.fastp/$p.fastp.R2.fq.gz  -5 10 -3 10  -w $cpu &> 2.fastp/$p.log  \n";
                print STDERR "  fastp  -i $f1  -I $f2  -o 2.fastp/$p.fastp.R1.fq.gz -O 2.fastp/$p.fastp.R2.fq.gz  -5 10 -3 10  -w $cpu  &> 2.fastp/$p.log \n";
        }
        close OUT;
        #command("ParaFly -c 2.fastp/command -CPU $cpu");
        $fastp_run = threads->create(\&com, "ParaFly -c 2.fastp/command -CPU $cpu");

        #open OUT, ">2.fastp.ok" or die ""; close OUT;
}else{
        print STDERR "Skip Step 2 for the file 2.fastp.ok exists \n";
}




#thread join
unless ($fastp_run eq ""){
        $fastp_run->join();
        open OUT, ">2.fastp.ok" or die ""; close OUT;


}

unless ($index_run eq ""){
        $index_run->join();
        open OUT ,">1.ATAC.index.ok" or die $!;  close OUT;

}


###ATAC command followed by zxx;
#     cd $INDIR1;$BOWTIE2 -N 1 -p 8 -q -I 10 -X 1000 --dovetail --no-unal --very-sensitive-local --no-mixed --no-discordant -x sc -1 $DATA/${i}.R1.clean.fq.gz -2 $DATA/${i}.R2.clean.fq.gz 2>${i}.map.log > ${i}.map.bam
#     &&$SAMTOOLS sort -@ 6 -O BAM -o ${i}.bam ${i}.map.bam&&rm ${i}.map.bam
#     &&$SAMTOOLS view -F 4 -u  -b -f 2 -q 30 -o ${i}.q30.bam ${i}.bam
#     &&$PICARD/picard MarkDuplicates PG=null VERBOSITY=ERROR QUIET=true CREATE_INDEX=false INPUT=${i}.q30.bam OUTPUT=${i}.rmdup.bam M=${i}.rmdup.log
#     &&$SAMTOOLS flagstat ${i}.bam > ${i}.stat
#     &&$SAMTOOLS depth ${i}.bam > ${i}.depth&&cat ${i}.depth | awk '{FS=" "}{sum+=$3} {if ($3>0)sites+=1} END {print "out",sum/NR,sites/419540624}' > ${i}.depth.cover.txt
#     &&$SAMTOOLS sort --threads 4 -o ${i}.final.sort.bam ${i}.rmdup.bam
#     &&$SAMTOOLS index ${i}.final.sort.bam;$MACS callpeak -f BAMPE -g 500000000 --keep-dup all -t ${i}.final.sort.bam -n ${i} --outdir $INDIR/macs --call-summits --SPMR -B 
#     &&$PICARD/picard CollectInsertSizeMetrics I=${i}.q30.bam O=${i}.q30.insersize.txt H=${i}.q30.insersize.pdf M=0.5
#Step3 bowtie alin
print STDERR "\n============================================\n";
print STDERR "Step 3: Bowtie align  " . "(" . (localtime) . ")" . "\n";
unless ( -e "3.Bowtie.ok") {	
        `mkdir 3.Bowtie` unless -e "3.Bowtie" && -d "3.Bowtie";
	chdir "3.Bowtie";
	open OUT, ">align.com" or die "";
	foreach my $p (sort keys %fq_fastp){
		my ($f1, $f2)= split/\t/, $fq_fastp{$p};
		print OUT "bowtie2    -N 1 -p $cpu -q -I 10 -X 1000 --dovetail --no-unal --very-sensitive-local --no-mixed --no-discordant -x   ../1.ATAC.index/index -1 $f1 -2 $f2 2> $p.align.log |samtools sort -@ 4 -O BAM -o $p.bam ; samtools view -F 4 -u -b -f2 -q 30 -o $p.algin.Q30bam  $p.bam \n";
		print STDERR "bowtie2    -N 1 -p $cpu -q -I 10 -X 1000 --dovetail --no-unal --very-sensitive-local --no-mixed --no-discordant -x   ../1.ATAC.index/index -1 $f1 -2 $f2 2> $p.align.log |samtools sort -@ 4 -O BAM -o $p.bam ; samtools view -F 4 -u -b -f2 -q 30 -o $p.algin.Q30bam  $p.bam \n";

	close OUT;
	command("ParaFly -c align.com -CPU 100");
	}


	unless( -e "picard.ok"){
		open PC, ">picard.com" or die $!;
		open DP, ">depth.com" or die "";
		open SO, ">sort.com" or die "";
		open MA, ">macs.com" or die "";
		open DU, ">du.com" or die "";
		open QF, ">q30.com" or die $!;
		foreach my $p (sort keys %fq_fastp){
			#picard MarkDuplicates PG=null VERBOSITY=ERROR QUIET=true CREATE_INDEX=false INPUT=${i}.q30.bam OUTPUT=${i}.rmdup.bam M=${i}.rmdup.log
			#samtools view -h MeJA-1.pairs.dedup.sort.bam|grep -v 'chrM'|grep -v 'chrC'|samtools sort -O bam  -@ 5 -o ->MeJA-1.final.bam
			print QF "samtools sort -@ 6 $p.bam | samtools view -F 4 -u -b -f2 -q 30 -o $p.align.Q30bam \n";
			print PC "picard MarkDuplicates PG=null VERBOSITY=ERROR QUIET=true CREATE_INDEX=false INPUT= $p.align.Q30bam OUTPUT=$p.align.Q30.rmdup.bam M=$p.rmdup.log\n";
			print DP "samtools flagstat  $p.bam > $p.stat  && samtools depth $p.bam > $p.depth  \n ";	
			print SO "samtools view -h  $p.align.Q30.rmdup.bam|grep -v 'ChrM'|grep -v 'ChrC' | samtools  sort --threads $cpu -o $p.final.sort.bam ;  samtools index $p.final.sort.bam;\n" ;
			`mkdir macs` unless -e "macs" && -d "macs";
			print MA "macs3 callpeak -f BAMPE -g $chr_len --keep-dup all -t $p.final.sort.bam -n $p --outdir macs  --call-summits --SPMR -B \n";
			print DU "picard CollectInsertSizeMetrics I=$p.algin.Q30bam O=$p.algin.Q30.insersize.txt H=$p.algin.Q30.insersize.pdf M=0.5\n";

		}
		close QF;
		close PC;
		close DP;
		close SO;
		close MA;
		close DU;
		command("ParaFly -c q30.com -CPU 100 ");
		command("ParaFly -c picard.com -CPU 100");
                command("ParaFly -c depth.com -CPU 100");
		foreach my $p (sort keys %fq_fastp){
                        open IN, "$p.depth" or die $!;
                        open OUT , ">$p.depth.cover.txt" or die "";
                        my ($num, $len,$total) =  (0,0,0);
                        while (<IN>){
                                chomp;
                                my @line=split/\t/, $_;
                                #$chr_len
                                $total++;
                                next unless $line[2] > 0;
                                $len+=$line[2];
                                $num++;
                        }
                        close IN;
                        #print OUT "out\t".$len/$total."\t".$num/$chr_len."\n";
			#$chr_len
                        close OUT;
		}

                command("ParaFly -c sort.com -CPU 100");
                command("ParaFly -c macs.com -CPU 100");
                command("ParaFly -c du.com  -CPU 100");

		open OUT, ">picard.ok" or die $!; close OUT;

	}
	


	chdir "../";
	open OUT, ">3.Bowtie.ok" or die ""; close OUT;
}else{
	print STDERR "Skip Step 3 for the file 3.Bowtie.ok exists \n ";
}

print STDERR "\n============================================\n";
print STDERR "Step 4: Downstream  " . "(" . (localtime) . ")" . "\n";
unless ( -e "4.Downstream.ok") {
#$chr_len;
	`mkdir 4.Downstream` unless -d "4.Downstream" or -e "4.Downstream";
	die $! unless chdir "4.Downstream";

#OV_ATAC_F_treat_pileup.bdg
#OV_ATAC_F_treat_pileup.bw
#OV_ATAC_F_treat_pileup.matrix_TSS.gz
#OV_ATAC_F_treat_pileup.matrix_TSS.PNG
#OV_ATAC_F_treat_pileup.sort_bdg
	unless (-e "heatmap.ok"){
		open OUT, ">heatmap.com" or die $!;
                foreach my $p (sort keys %fq_fastp){
			print OUT "bedSort  ../3.Bowtie/macs/${p}_treat_pileup.bdg stdout | bedClip  stdin ../../$size stdout | perl -ane 'print if(\$F[1]<\$F[2]) '> $p.tmp.bdg;bedGraphToBigWig $p.tmp.bdg ../../$size ${p}_treat_pileup.bw; computeMatrix reference-point --referencePoint TSS  -p 6 -a 3000 -b 3000 -R ../../$tss_bed -S ${p}_treat_pileup.bw -o ${p}_matrix_TSS.gz --missingDataAsZero ; plotHeatmap --heatmapHeight 16 --heatmapWidth 4  --colorMap Reds  --legendLocation none  --samplesLabel $p -m ${p}_matrix_TSS.gz  -o ${p}.heatmap.PNG  \n ";
		
		}
		close OUT;
		command("ParaFly -c heatmap.com -CPU $cpu");
		open OUT, ">heatmap.ok" or die $!; close OUT;
	}


}else{
        print STDERR "Skip Step 4 for the file 4.Downstream.ok exists \n ";
}





sub com {
	my $com=shift;
	command("$com");
	return(" $com successful!\n");
}
