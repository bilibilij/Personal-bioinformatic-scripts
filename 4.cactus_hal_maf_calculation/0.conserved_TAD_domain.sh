export PATH=$PATH:/usr_storage/software/cactus/bin
export PYTHONPATH=/usr_storage/software/cactus/submodules/:$PYTHONPATH
source /usr_storage/software/cactus/cactus_env/bin/activate
export PATH="$PATH:/usr_storage/software/bedtools2/bin/"
export PATH="$PATH:/usr_storage/jcf/3.Pangenome/3.TAD_liftover/Liftover_one_by_one/filter_chr/script/"

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

1.split_bed_for_liftover.pl  $qry_domain  $qry $ref $hal $cpu

ParaFly -c ${qry}_${ref}/command/${qry}_qry_$ref.liftover  -CPU $cpu

#2. each dataset (per one TAD domain two liftover of boundaries and one liftover of TAD body) contains several small liftover features because of speci-specific indel , TE transposon and structral variation, so we merge distance of features of boundaries  under 10kb, body under 100kb.

[ -f ${qry}_${ref}/command/${qry}_qry_$ref.merge_domain.com.completed ] && rm -f  ${qry}_${ref}/command/${qry}_qry_$ref.merge_domain.com.completed   ${qry}_${ref}/command/${qry}_qry_$ref.merge_domain.com
[ -f ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com.completed ] && rm -f ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com.completed ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com

awk '{print "cat "$NF"|sort -k1,1 -k2,2n |bedtools merge -d 100000 > "$NF".merge100k"}' ${qry}_${ref}/command/${qry}_qry_$ref.liftover |grep domain > ${qry}_${ref}/command/${qry}_qry_$ref.merge_domain.com
awk '{print "cat "$NF"|sort -k1,1 -k2,2n |bedtools merge -d 10000 > "$NF".merge10k" }'  ${qry}_${ref}/command/${qry}_qry_$ref.liftover |grep bound> ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com

ParaFly -c ${qry}_${ref}/command/${qry}_qry_$ref.merge_domain.com -CPU $cpu
ParaFly -c ${qry}_${ref}/command/${qry}_qry_$ref.merge_bound.com -CPU $cpu

#3. there are several liftovered TAD domain of one query feature, this step we choose TAD domain which have also support by liftovered boundary A and B to perform next analysis.
TAD_num=`wc -l $qry_domain | awk '{print $1}' `

for ((i=1;i<=$TAD_num;i++)) ;do paste  ${qry}_${ref}/$qry.TAD_domain.$i.bed ${qry}_${ref}/$qry.TAD_boundaryA.$i.bed ${qry}_${ref}/$qry.TAD_boundaryB.$i.bed $qry_domain >${qry}_${ref}/$qry.TAD.$i  ;done


[ -f ${qry}_${ref}/command/$qry.bound_pass.com.completed ] && rm -f ${qry}_${ref}/command/$qry.bound_pass.com.completed ${qry}_${ref}/command/$qry.bound_pass.com

for ((i=1;i<=$TAD_num;i++)) ;do echo "bedtools intersect -a ${qry}_${ref}/$qry.TAD_domain.$i.liftover.merge100k -b ${qry}_${ref}/$qry.TAD_boundaryA.$i.liftover.merge10k ${qry}_${ref}/$qry.TAD_boundaryB.$i.liftover.merge10k -wo  2>/dev/null  >${qry}_${ref}/$qry.TAD_domain.$i.liftover.merge100k.intersect_bound; 2.intersect_bound_pass.pl ${qry}_${ref}/$qry.TAD_domain.$i.liftover.merge100k.intersect_bound >${qry}_${ref}/$qry.TAD_domain.$i.liftover.merge100k.intersect_bound_pass" >> ${qry}_${ref}/command/$qry.bound_pass.com;done

#Pdav.bound_pass.com


ParaFly -c ${qry}_${ref}/command/$qry.bound_pass.com -CPU $cpu


#4. we define conserved domain between two speicies should not exceed 50%, boundaries not exceed 100% to get final liftovered TAD domain and boundaries.
for ((i=1;i<=$TAD_num;i++)) ; do  3.TAD_body_size_pass.pl ${qry}_${ref}/$qry.TAD_domain.$i.bed  ${qry}_${ref}/$qry.TAD_domain.$i.liftover.merge100k.intersect_bound_pass >${qry}_${ref}/$qry.TAD_domain.$i.liftover.merge100k.intersect_bound_pass.size_pass ;done


#5. To determine conserved TAD domain, liftovered TAD domain should overlap reference TAD domain exceed 0.8 in both species and also contain any boundary overlap.

[ -d ${qry}_${ref}_ortho_TAD ] && rm -fr ${qry}_${ref}_ortho_TAD
mkdir ${qry}_${ref}_ortho_TAD


ls ${qry}_${ref} |grep size_pass|while read id ;do wc -l ${qry}_${ref}/$id ;done | awk '{if ($1==2) print $2}'  > ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_2_pass.list
ls ${qry}_${ref} |grep size_pass|while read id ;do wc -l ${qry}_${ref}/$id ;done | awk '{if ($1==4) print $2}'  > ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_4_pass.list
ls ${qry}_${ref} |grep size_pass|while read id ;do wc -l ${qry}_${ref}/$id ;done | awk '{if ($1==6) print $2}' > ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_6_pass.list
ls ${qry}_${ref} |grep size_pass|while read id ;do wc -l ${qry}_${ref}/$id ;done | awk '{if ($1==8) print $2}' > ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_8_pass.list


cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_2_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.9 -F 0.9 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.size2.0.9.tad
cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_2_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.8 -F 0.8 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.size2.0.8.tad
cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_2_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.7 -F 0.7 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.size2.0.7.tad
cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_2_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.6 -F 0.6 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.size2.0.6.tad
cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_2_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.5 -F 0.5 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.size20.5.tad


ls ${qry}_${ref} |grep size_pass|while read id ;do wc -l ${qry}_${ref}/$id ;done | awk '{if ($1!=0) print $2}' >${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_pass.list

cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.9 -F 0.9 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.9.tad
cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.8 -F 0.8 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.8.tad
cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.7 -F 0.7 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.7.tad
cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.6 -F 0.6 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.6.tad
cat ${qry}_${ref}_ortho_TAD/${qry}_${ref}.size_pass.list  |while read id; do bedtools intersect -a $ref_domain -b $id -f 0.5 -F 0.5 2>/dev/null -wo ;done >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.5.tad



4.boundary_overlap.pl ${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.9.tad >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.9.tad.boundary_pass
4.boundary_overlap.pl ${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.8.tad >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.8.tad.boundary_pass
4.boundary_overlap.pl ${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.7.tad >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.7.tad.boundary_pass
4.boundary_overlap.pl ${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.6.tad >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.6.tad.boundary_pass
4.boundary_overlap.pl ${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.5.tad >${qry}_${ref}_ortho_TAD/${qry}_${ref}_ortho.0.5.tad.boundary_pass




