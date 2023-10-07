ORTHOFINDER_HOME=$1
ORTHOGROUPS=$2
NBTHREADS=$3
TEMPFOLDER=$4
ORTHOFINDER_RESULTS=$(dirname $(dirname $ORTHOGROUPS))
TEMP=$(mktemp -d -p $TEMPFOLDER)
for F in ${@:1,2,3,4,5}; do
 cp $F ${TEMP}
done
${ORTHOFINDER_HOME} -T fasttree -M msa -f ${TEMP} -S diamond -t ${NBTHREADS} -a ${NBTHREADS} 
mv ${TEMP}/OrthoFinder/Results*/* ${ORTHOFINDER_RESULTS}
rm -r ${TEMP}
