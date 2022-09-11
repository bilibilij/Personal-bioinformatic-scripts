
bed=$1
size=$2
total=`cat $bed|wc -l `
name=$3

 awk '{print $1"-"$2}' OFS="\t" $bed  | uniq -c |sed "s/\ /\t/g" | awk '{print $NF,$(NF-1)}' OFS="\t"|  awk '{print $1, '$size'*$2/'$total'}'  OFS="\t"  |sort -k1 |  awk '{seen[$1]+=$2}END{for (i in seen) print i, seen[i]}' OFS="\t" |sed "s/-/\t/g" | sort -n -k2 > $name.bg


