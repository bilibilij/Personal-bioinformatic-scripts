
interpro=$1
total_pep=$2


grep PF00931 $interpro | awk '{print $1}' |sort |uniq > ${interpro%.*}.nlr.list

~/scripts/extract_fa_from_id.pl $total_pep ${interpro%.*}.nlr.list > ${interpro%.*}.nlr.fasta


java -jar /usr_storage/software/NLR-Annotator/NLR-Parser3.jar  -i ${interpro%.*}.nlr.fasta  -x meme.xml -y /usr_storage/jcf/.conda/envs/meme/bin/mast -o ${interpro%.*}.nlr.txt -g ${interpro%.*}.nlr.gff -b ${interpro%.*}.nlr.bed -c ${interpro%.*}.nlr.xml -t 20




