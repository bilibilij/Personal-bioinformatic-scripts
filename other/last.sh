genome1=$1
genome2=$2
basename1=`basename $genome1`
index=${basename1%%.*}
basename2=`basename $genome2`
index2=${basename2%%.*}


source /usr_storage/jcf/2.zhugecai/4.last/.bashrc
#lastdb -P0 -uNEAR -R01 $index $genome1
#ln -s $genome1
last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 $index $genome2  > $index-$index2.mat
lastal -m50 -E0.05 -C2 -p $index-$index2.mat  $index  $genome2 |last-split -m1 > $index-$index2.maf


maf-swap   $index-$index2.maf  |awk '/^s/ {$2 = (++s % 2 ? "" : "") $2}1'|last-split -m1 |maf-swap > $index-$index2-2.maf
#maf-swap  $index-$index2.maf |awk '/^s/ { = (++s % 2 ? Ovio. : Ovio.) }1'|last-split -m1 |maf-swap > O.vio-self2.maf

last-postmask $index-$index2-2.maf |maf-convert -n tab |awk -F'=' '$2 <= 1e-5'  > $index-$index2.tab
#less Ovio-self.tab|grep -v scaffold > Ovio-self-chr.tab
#sed 's/Ovio\.//g' Ovio-self-chr.tab >Ovio-self-chr-noOv.tab
last-dotplot  $index-$index2.tab  $index-$index2.PNG

