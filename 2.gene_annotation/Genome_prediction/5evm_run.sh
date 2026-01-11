##0.load environment
source ~/lsy/bashrc
#weight="~/software/EVidenceModeler/testing/changfu_weights.txt"
weight="/SciBorg/array0/changfu/software/EVidenceModeler/testing/changfu_weights.txt"

###weight file like this, second line should be consistent with the second line of corresponse gff
#PROTEIN GENEWISE    5
#TRANSCRIPT   pasa_transcript_alignments   10
#ABINITIO_PREDICTION   AUGUSTUS    1


##1.Read input

#Initial genome assembly without repeatmasked
genome=$1
#repeatmasker gff initial files from geta2.7
repeat_gff=$2
#augustus file from geta2.7
ab=$3
#homologous file from geta2.7
homo=$4
#pasa file assembly.gff3
pasa=$5

##2.preprocess of input
mkdir input/;
dir=`realpath input`
#Repeat process
perl -pe 's/.*Low_complexity.*//; s/.*Simple_repeat.*//; s/^#.*//; s/^\s*S//;' $repeat_gff > input/Repeat.gff3
#PASA process
perl -pe 's/\t\S+/\tpasa_transcript_alignments/' $pasa >input/PASA.gff3
#ab process
augustus_for_evm_gff3_new.pl $ab > input/ab.gff3 
#Just copy genewise file
cp $homo input/homo.gff3

##3.split genomic region into small pieces and write down the commands
mkdir split;
ref=`realpath $1`
#partition_EVM_inputs.pl --genome $ref  --gene_predictions $dir/ab.gff3 --protein_alignments $dir/homo.gff3   --repeats $dir/Repeat.gff3  --transcript_alignments $dir/PASA.gff3  --segmentSize 5000000 --overlapSize 10000 --partition_listing split/partition.list   --partition_dir  split;

write_EVM_commands.pl --genome $ref  --gene_predictions $dir/ab.gff3 --protein_alignments $dir/homo.gff3   --repeats $dir/Repeat.gff3  --transcript_alignments $dir/PASA.gff3 --output_file_name evm.out --weights $weight  --partitions split/partition.list   > evm.command.list

##4.run and obtain final gff3
ParaFly -c evm.command.list -CPU 48
recombine_EVM_partial_outputs.pl --partitions split/partition.list  --output_file_name evm.out 
convert_EVM_outputs_to_GFF3.pl --partitions split/partition.list  --output_file_name evm.out --genome $ref
cat split/*/evm.out.gff3 > evm.out.gff3
rename_and_sort_gff.pl evm.out.gff3 evm > evm.gff3



