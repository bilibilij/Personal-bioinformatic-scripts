#software dependecies
bedtools="/usr_storage/software/bedtools2/bin/bedtools"


#parameters or input file 
### the gene body and flanking bed file is extracted by gff3 file with a home-make script, using gff3 with a different 9th line may encouter some error, please checking the gff if is correct. 
##The example gff format ; gene_id =~ m/ID=([^;]*)/;
#Chr01   .       gene    2894    5929    0.23    -       .       ID=Ovio00001;Name=Ovio00001;



gff=$1
flank=$2
bin=$3
prefix=$4
meth=$5
out=$6

#step1: bin examination
#$bedtools  makewindows -b $bed -n $bin   -i srcwinnum > $prefix.bin.bed

perl /home/jcf/scripts/gene_flank_gff_to_bed.pl $gff $flank >$prefix.flank.bed
$bedtools  makewindows -b $prefix.flank.bed -n $bin -i srcwinnum   > $prefix.flank.bin.bed
$bedtools intersect -a $prefix.flank.bin.bed -b $meth -wb 2>/dev/null>$prefix.flank.bin.meth
perl ./filter_ID.pl tpm1_pass_ID.list $prefix.flank.bin.meth > $prefix.flank.bin.meth.filtered
#step2: calculate meth level

perl /home/jcf/scripts/meth_bin_calculate.pl $prefix.flank.bin.meth.filtered  $bin ${prefix%%_*} |sort -k1,1n >>  $out




