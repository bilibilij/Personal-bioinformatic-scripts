ls 0.data/* |while read id ;do echo "nanoQC $id -o ${id%/*}/${id##*/}_QCout ";done |sed 's/fastq.gz_//g' >0.nanoQC.com

