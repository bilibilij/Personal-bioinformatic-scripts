export PATH="$PATH:/usr_storage/software/minimap2"

genome=$1
hifi=$2

minimap2 -ax   map-hifi $genome $hifi  -t 80  --secondary=no --split-prefix ref |samtools sort -@ 12 -o $genome.aligned.bam -T tmp.ali.${genome}


#minimap2 -ax map-pb -t 50  $genome $clr  --secondary=no --split-prefixref | samtools sort -@ 12 -m 1G -o aligned.bam -T tmp.ali
