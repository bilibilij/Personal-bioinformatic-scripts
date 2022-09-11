export PATH=$PATH:/usr_storage/software/cactus/bin
export PYTHONPATH=/usr_storage/software/cactus/submodules/:$PYTHONPATH
source /usr_storage/software/cactus/cactus_env/bin/activate
export PATH="$PATH:/usr_storage/software/bedtools2/bin/"
#export PATH="$PATH:/usr_storage/jcf/3.Pangenome/3.TAD_liftover/Liftover_one_by_one/filter_chr/script/"
export PATH="$PATH:/usr_storage/jcf/3.Pangenome/3.TAD_liftover/Liftover_one_by_one/filter_chr/wjl"


qry=$1
ref=$2
hal=$3
cpu=$4


qry_domain=${qry}.TAD_domain.bed
ref_domain=${ref}.TAD_domain.bed


[ -d ${qry}_${ref} ] && rm -rf ${qry}_${ref}
mkdir -p ${qry}_${ref}/command
#1.split TAD domain and two boundaries of this one TAD domain for one to one feature liftover. TAD boundaries size default as 20kb.
[ -f ${qry}_${ref}/command/${qry}_qry_$ref.liftover ] && rm -rf ${qry}_${ref}/command/${qry}_qry_$ref.liftover  ${qry}_${ref}/command/${qry}_qry_$ref.liftover.completed

1.split_bed_for_liftover.pl.wjl $qry_domain  $qry $ref $hal $cpu

ParaFly -c ${qry}_${ref}/command/${qry}_qry_$ref.liftover  -CPU $cpu

#2. each dataset (per one TAD domain two liftover of boundaries and one liftover of TAD body) contains several small liftover features because of speci-specific indel , TE transposon and structral variation, so we merge distance of features of boundaries  under 10kb, body under 100kb.

[ -f ${qry}_${ref}/command/${qry}_qry_$ref.merge_domain.com.completed ] && rm -f  ${qry}_${ref}/command/${qry}_qry_$ref.merge_domain.com.completed   ${qry}_${ref}/command/${qry}_qry_$ref.merge_domain.com
[ -f ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com.completed ] && rm -f ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com.completed ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com

awk '{print "cat "$NF"|sort -k1,1 -k2,2n |bedtools merge -d 10000 > "$NF".merge10k" }'  ${qry}_${ref}/command/${qry}_qry_$ref.liftover |grep bound> ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com

ParaFly -c ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com -CPU $cpu


[ -f ${qry}_${ref}.boundary.pass ] && rm ${qry}_${ref}.boundary.pass

TAD_num=`wc -l $qry_domain | awk '{print $1}' `

for  ((i=1;i<=$TAD_num;i++)) ;do 4.boundary_overlap.pl.wjl ${qry}_${ref}/$qry.TAD_domain.${i}.bed  ${qry}_${ref}/$qry.TAD_boundaryA.$i.liftover.merge10k ${qry}_${ref}/$qry.TAD_boundaryB.$i.liftover.merge10k >>${qry}_${ref}.boundary.pass0 ;done

#Pade_Pkor/Pade.TAD_boundaryA.510.liftover.merge10k


sort -k1 ${qry}_${ref}.boundary.pass0 |uniq >${qry}_${ref}.boundary.pass

#3. there are several liftovered TAD domain of one query feature, this step we choose TAD domain which have also support by liftovered boundary A and B to perform next analysis.
#TAD_num=`wc -l $qry_domain | awk '{print $1}' `

#for ((i=1;i<=$TAD_num;i++)) ;do paste  ${qry}_${ref}/$qry.TAD_domain.$i.bed ${qry}_${ref}/$qry.TAD_boundaryA.$i.bed ${qry}_${ref}/$qry.TAD_boundaryB.$i.bed $qry_domain >${qry}_${ref}/$qry.TAD.$i  ;done

