
import os

ALL_SAMPLE = [key for key, value in config["ONT"].items()]

ALL_RNA = [key for key, value in config["rna"].items()]
genome=config["genome"]
genome_prefix=os.path.split(genome)[1]
EDTA= config["EDTA"]
sp=config["sp"]
FASTP=config["FASTP"]
SAMTOOLS=config["SAMTOOLS"]
HISAT_DIR=config["HISAT_DIR"]


rule all:
    input:
        "0.EDTA/"+genome_prefix+".mod.EDTA.TEanno.sum",
        #["1.RNA_based/1.ONT_preprocessing/"+str(sample)+".cleaned.fastq" for sample in ALL_SAMPLE],
        expand("3.geta/geta.tmp/2.hisat2/genome.{num}.ht2", num=range(1,8)),
        "3.geta/geta.tmp/2.hisat2/hisat2.sorted.bam",
        "1.RNA_based/2.trinity/trinity_out_dir/Trinity-GG.fasta",
        "4.evm/evm.gff3",
        "5.update/PASA_update.gff3",
#        "3.geta/geta.tmp/1.trimmomatic/{sample}.R1.cleaned.fastq",
#        "3.geta/geta.tmp/1.trimmomatic/{sample}.R2.cleaned.fastq"
        ["3.geta/geta.tmp/1.trimmomatic/"+str(sample)+".R1.cleaned.fastq" for sample in ALL_RNA],
        ["3.geta/geta.tmp/1.trimmomatic/"+str(sample)+".R2.cleaned.fastq" for sample in ALL_RNA],
#        r1t = "3.geta/geta.tmp/1.trimmomatic/{sample}.R1.cleaned.fastq",


#["0.fastp/"+str(sample)+"r2.trimmed.fq.gz" for sample in ALL_SAMPLE],



rule EDTA:
    input:
        genome
    conda:
        "envs/EDTA.yml"
    output:
        "0.EDTA/"+genome_prefix+".mod.EDTA.TEanno.sum",
    threads:
        20
    shell:
        """
        export PERL5LIB=/
        cd 0.EDTA/;
        perl {EDTA} --genome {input} --sensitive 1 --anno 1 --threads {threads}
        """

#rule ONT:
#    input:
#        fq = lambda w: config["ONT"][w.sample],
#    output:
#        trim_fq = "1.RNA_based/1.ONT_preprocessing/{sample}.cleaned.fastq"
#    conda:
#        "envs/pychopper.yaml"
#    threads:
#        20
#    shell:
#        """
#        pychopper -r 1.RNA_based/1.ONT_preprocessing/{wildcards.sample}.report.pdf -u 1.RNA_based/1.ONT_preprocessing/{wildcards.sample}.unclassified.fq -w 1.RNA_based/1.ONT_preprocessing//{wildcards.sample}.rescued.fq -t {threads}  {input.fq} {output.trim_fq}
#       """



#j=" "
#fqL=j.join(["1.RNA_based/1.ONT_preprocessing/"+str(sample)+".cleaned.fastq" for sample in ALL_SAMPLE])


#rule fqfa:
#    input:
#        ["1.RNA_based/1.ONT_preprocessing/"+str(sample)+".cleaned.fastq" for sample in ALL_SAMPLE]
#    params:
#        fq= fqL
#    output:
#        fid="1.RNA_based/1.ONT_preprocessing/flnc.id",
#        ffa="1.RNA_based/1.ONT_preprocessing/ONT_cleaned.fasta"
#    shell:
#        """
#        sh scripts/fq_to_fq.sh   {fqL} >{output.ffa}
#        perl scripts/fa_id.pl  1.RNA_based/1.ONT_preprocessing/ONT_cleaned.fasta  >{output.fid} 
#        """


rule fastp:
    input:
        r1 = lambda w: config["rna"][w.sample]['r1'],
        r2 = lambda w: config["rna"][w.sample]['r2']
    output:
        r1t = "3.geta/geta.tmp/1.trimmomatic/{sample}.R1.cleaned.fastq",
        r2t = "3.geta/geta.tmp/1.trimmomatic/{sample}.R2.cleaned.fastq",
    log:
        "3.geta/geta.tmp/1.trimmomatic/{sample}.log"
    threads: 8 
    shell:
        """
        {FASTP} -i {input.r1} -I {input.r2} -o {output.r1t} -O {output.r2t} --length_required=35 --cut_mean_quality 20 -5 20 -3 20 -w {threads} 2> {log}
        touch 3.geta/geta.tmp/1.trimmomatic.ok
        """



rule merge_fq:
    input:
        ["3.geta/geta.tmp/1.trimmomatic/"+str(sample)+".R1.cleaned.fastq" for sample in ALL_RNA],
        ["3.geta/geta.tmp/1.trimmomatic/"+str(sample)+".R2.cleaned.fastq" for sample in ALL_RNA],
    output:
        r1 = "3.geta/geta.tmp/1.trimmomatic/reads1.fastq",
        r2 = "3.geta/geta.tmp/1.trimmomatic/reads2.fastq"
    shell:
        """
        cat 3.geta/geta.tmp/1.trimmomatic/*R1.cleaned.fastq > 3.geta/geta.tmp/1.trimmomatic/reads1.fastq
        cat 3.geta/geta.tmp/1.trimmomatic/*R2.cleaned.fastq > 3.geta/geta.tmp/1.trimmomatic/reads2.fastq
        """


rule hisat_index:
    input:
        genome
    output:
        expand("3.geta/geta.tmp/2.hisat2/genome.{num}.ht2", num=range(1,8))
    threads: 8
    shell:
        """
        {HISAT_DIR}/hisat2-build -p {threads} {genome} 3.geta/geta.tmp/2.hisat2/genome
        """

rule hisat_align:
    input:
        r1 = "3.geta/geta.tmp/1.trimmomatic/reads1.fastq",
        r2 = "3.geta/geta.tmp/1.trimmomatic/reads2.fastq",
        ind = expand("3.geta/geta.tmp/2.hisat2/genome.{num}.ht2", num=range(1,8))
    output:
        bam = "3.geta/geta.tmp/2.hisat2/hisat2.sorted.bam",
        sam = "3.geta/geta.tmp/2.hisat2/hisat2.sorted.sam"
    threads: 20
    log:
        "3.geta/geta.tmp/2.hisat2/hisat_align.log"
    shell:
        """
        export PATH="/usr/bin/:$PATH"

        {HISAT_DIR}/hisat2 -p {threads} -x 3.geta/geta.tmp/2.hisat2/genome -1 {input.r1} -2 {input.r2} --min-intronlen 20 --max-intronlen 20000 --dta --score-min L,0.0,-0.4 2>{log} | {SAMTOOLS} sort -O BAM  -@ 4 -T 3.geta/geta.tmp/2.hisat2/temp > {output.bam}
        {SAMTOOLS} index {output.bam}
        {SAMTOOLS} view 3.geta/geta.tmp/2.hisat2/hisat2.sorted.bam > 3.geta/geta.tmp/2.hisat2/hisat2.sorted.sam
        touch 3.geta/geta.tmp/2.hisat2.ok
        """

rule trinity_d:
    input:
        r1 = "3.geta/geta.tmp/1.trimmomatic/reads1.fastq",
        r2 = "3.geta/geta.tmp/1.trimmomatic/reads2.fastq"
    output:
        "1.RNA_based/2.trinity/trinity_out_dir_D/Trinity.fasta",
        "1.RNA_based/2.trinity/flnc.id"
    threads:
        80
    conda:
        "envs/trinity.yaml"
    log:
        "1.RNA_based/2.trinity/Trinity_DD.log"
    shell:
        """
        export PERL5LIB=/
        Trinity --seqType fq --max_memory 50G --left 3.geta/geta.tmp/1.trimmomatic/reads1.fastq  --right 3.geta/geta.tmp/1.trimmomatic/reads2.fastq  --CPU {threads} --output 1.RNA_based/2.trinity/trinity_out_dir_D/ --no_normalize_reads --bflyCalculateCPU 
        perl scripts/fa_id.pl  1.RNA_based/2.trinity/trinity_out_dir_D/Trinity.fasta  >  1.RNA_based/2.trinity/flnc.id
        cp 1.RNA_based/2.trinity/trinity_out_dir_D/Trinity.fasta 1.RNA_based/2.trinity/
        cp 1.RNA_based/2.trinity/flnc.id 1.RNA_based/

        """

rule trinity:
    input:
        bam = "3.geta/geta.tmp/2.hisat2/hisat2.sorted.bam"
    output:
        "1.RNA_based/2.trinity/trinity_out_dir/Trinity-GG.fasta"
    conda:
        "envs/trinity.yaml"
    log:
        "1.RNA_based/2.trinity/Trinity_GG.log"
    threads:
        20
    shell:
        """
        export PERL5LIB=/
        Trinity --genome_guided_bam {input.bam} --genome_guided_max_intron 10000 --max_memory 60G --no_normalize_reads --bflyCalculateCPU --CPU {threads} --output 1.RNA_based/2.trinity/trinity_out_dir/

        """
            


rule tdn:
    input:
        GG= "1.RNA_based/2.trinity/trinity_out_dir/Trinity-GG.fasta",
#        ONT= "1.RNA_based/1.ONT_preprocessing/ONT_cleaned.fasta"
        D="1.RNA_based/2.trinity/trinity_out_dir_D/Trinity.fasta"
    output:
        trans = "1.RNA_based/3.PASA/transcripts.fasta",
        fid="1.RNA_based/3.PASA/flnc.id",
    shell:
        """
        cat {input.GG} {input.D} > {output.trans}
        perl scripts/fa_id.pl  {input.D}  >{output.fid}
        """

rule PASA_seq:
    input:
        trans = "1.RNA_based/3.PASA/transcripts.fasta"
    output:
        "1.RNA_based/3.PASA/transcripts.fasta.clean"
    shell:
        """
        source envs/pasa.bashrc
        export PERL5LIB=/

        seqclean {input.trans} -o {output}
        mv transcripts* 1.RNA_based/3.PASA/
        """
#mv transcript* 1.RNA_based/3.PASA/

rule PASA:
    input:
        #trans = "1.RNA_based/3.PASA/transcripts.fasta",
        trans = "1.RNA_based/3.PASA/transcripts.fasta.clean",
        genome=genome,
    output:
        "1.RNA_based/3.PASA/"+sp+".pasa_assemblies.gff3"
    threads:
        40
    log:
        "PASA_aln.log"
    shell:
        """
        source envs/pasa.bashrc
        export PERL5LIB=/
        export PATH="/usr/local/bin:$PATH"
        sp={sp}

        echo "DATABASE=$sp
#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=<__MIN_PERCENT_ALIGNED__>
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=<__MIN_AVG_PER_ID__>

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50
        " > 1.RNA_based/3.PASA/aln.config
        cd 1.RNA_based/3.PASA/;
        #seqclean  transcripts.fasta
        Launch_PASA_pipeline.pl -c  /usr_storage/jcf/2.zhugecai/00.diploid/2.anno/1.RNA_based/3.PASA/aln.config -C -r -R -T -g {genome} -t /usr_storage/jcf/2.zhugecai/00.diploid/2.anno/1.RNA_based/3.PASA//transcripts.fasta.clean -u /usr_storage/jcf/2.zhugecai/00.diploid/2.anno/1.RNA_based/3.PASA//transcripts.fasta --ALIGNERS gmap,blat  --CPU {threads} --TDN /usr_storage/jcf/2.zhugecai/00.diploid/2.anno/1.RNA_based/3.PASA/flnc.id 2>{log}
        cd ../../
        """

geta_dir=config["geta"]
homo=config["homo"]


rule geta:
    input:
        genome=genome,
        config=geta_dir+"conf.txt",
        hisat="3.geta/geta.tmp/2.hisat2/hisat2.sorted.bam"
    output:
        "3.geta/geta.tmp/5.augustus.ok",
        #"3.geta/geta.tmp/5.augustus/augustus.gff3",
        #"3.geta/geta.tmp/4.homolog/genewise.gff3",
        "3.geta/geta.tmp/4.homolog.ok",
        #"3.geta/geta.tmp/0.repeat/genome.repeat.gff3",
#        "3.geta/geta.tmp/0.RepeatMasker.ok"
    conda:
        "envs/geta.yaml"
    log:
        "3.geta/geta.log"
    threads:
        20
    shell:
        """
        source envs/geta.bashrc
        export PERL5LIB=/
        geta_anno.pl --RM_species Embryophyta --out_prefix 3.geta/geta --config {input.config} --cpu {threads} --protein {homo} -genome {genome} -1 a -2 a --augustus_species {sp}  2>{log}

        """

pwd=os.getcwd()

rule EVM:
    input:
        "3.geta/geta.tmp/5.augustus.ok",
        "3.geta/geta.tmp/4.homolog.ok",
#        "3.geta/geta.tmp/0.RepeatMasker.ok"
        "1.RNA_based/3.PASA/"+sp+".pasa_assemblies.gff3"

#        augustus="3.geta/geta.tmp/5.augustus/augustus.gff3",
#        genewise="3.geta/geta.tmp/4.homolog/genewise.gff3",
#        pasa="1.RNA_based/3.PASA/"+sp+".pasa_assemblies.gff3",
#        repeat="3.geta/geta.tmp/0.repeat/genome.repeat.gff3",
    output:
        "4.evm/evm.gff3"
    log:
        "4.evm/evm.log"
    threads:
        50
    shell:
        """
        source envs/evm.bashrc
        export PERL5LIB=/

        #perl scripts/augustus_for_evm_gff3_new.pl {input.augustus} > 4.evm/ab.gff3
        perl scripts/augustus_for_evm_gff3_new.pl 3.geta/geta.tmp/5.augustus/augustus.gff3 > 4.evm/ab.gff3
        #cp {input.genewise} 4.evm/genewise.gff3
        cp 3.geta/geta.tmp/4.homolog/genewise.gff3 4.evm/genewise.gff3
 #       perl -pe 's/.*Low_complexity.*//; s/.*Simple_repeat.*//; s/^#.*//; s/^\s*S//;' {input.repeat} > 4.evm/repeat.gff3
        perl -pe 's/.*Low_complexity.*//; s/.*Simple_repeat.*//; s/^#.*//; s/^\s*S//;'  3.geta/geta.tmp/0.repeat/genome.repeat.gff3  > 4.evm/repeat.gff3
#        perl -p -i -e 's/\t\S+/\tpasa_transcript_alignments/'  {input.pasa} > 4.evm/pasa.gff3
        perl -p -i -e 's/\t\S+/\tpasa_transcript_alignments/'  1.RNA_based/3.PASA/"+sp+".pasa_assemblies.gff3 > 4.evm/pasa.gff3
        mkdir 4.evm/split
        partition_EVM_inputs.pl --genome {genome}  --gene_predictions {pwd}/4.evm/ab.gff3 --protein_alignments {pwd}/4.evm/genewise.gff3 \
                --repeats {pwd}/4.evm/repeat.gff3  --transcript_alignments {pwd}/4.evm/pasa.gff3  --segmentSize 5000000 --overlapSize 10000 --partition_listing 4.evm/split

echo "#-------------weight.txt---------------#
ABINITIO_PREDICTION	AUGUSTUS	1
PROTEIN	GeneWise	5
TRANSCRIPT	pasa_transcript_alignments	10
#-------------weight.txt---------------#
" > 4.evm/weight.txt

        write_EVM_commands.pl --genome {genome}  --gene_predictions {pwd}/4.evm/ab.gff3 --protein_alignments {pwd}/4.evm/genewise.gff3 \
                --repeats {pwd}/4.evm/repeat.gff3 --transcript_alignments {pwd}/4.evm/pasa.gff3 --output_file_name evm.out --weights {pwd}/4.evm/weight.txt >4.evm/command.list
        {ParaFly} -c 4.evm/command.list -CPU {threads}

        recombine_EVM_partial_outputs.pl --partitions 4.evm/split --output_file_name evm.out 
        convert_EVM_outputs_to_GFF3.pl --partitions 4.evm/split --output_file_name evm.out --genome {genome}

        cat 4.evm/split/*/evm.out.gff3 > 4.evm/evm.out.gff3

        perl scripts/rename_and_sort_gff.pl 4.evm/evm.out.gff3 evm > 4.evm/evm.gff3

        """
 
 
rule PASA_update:
    input:
        evm = "4.evm/evm.gff3",
    output:
        "5.update/PASA_update.gff3"
    threads:
        20
    log:
        "5.update/PASA_update.log"
    shell:
        """
        source envs/pasa.bashrc
        sp={sp}
        export PERL5LIB=/

echo"
## templated variables to be replaced exist as <__var_name__>
# database settings
DATABASE=$sp
#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.


#script cDNA_annotation_comparer.dbi
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP=<__MIN_PERCENT_OVERLAP__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING=<__MIN_PERCENT_PROT_CODING__>
cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE=<__MIN_PERID_PROT_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE=<__MIN_PERCENT_LENGTH_FL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE=<__MIN_PERCENT_LENGTH_NONFL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_FL_ORF_SIZE=<__MIN_FL_ORF_SIZE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH=<__MIN_PERCENT_ALIGN_LENGTH__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE=<__MIN_PERCENT_OVERLAP_GENE_REPLACE__>
cDNA_annotation_comparer.dbi:--STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE=<__STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE__>
cDNA_annotation_comparer.dbi:--TRUST_FL_STATUS=<__TRUST_FL_STATUS__>
cDNA_annotation_comparer.dbi:--MAX_UTR_EXONS=<__MAX_UTR_EXONS__>
cDNA_annotation_comparer.dbi:--GENETIC_CODE=<__GENETIC_CODE__>
" > 5.update/update.config
        cd 5.update/;
        Launch_PASA_pipeline.pl -c update.config -A -T -L -g {genome} -t ../1.RNA_based/3.PASA/transcripts.fasta.clean -u ../1.RNA_based/3.PASA/transcripts.fasta --annots ${input.evm}
        ll -tr |grep gff3$  |tail -n 1 |awk '{print $NF}' |while read id ;do Launch_PASA_pipeline.pl -c update.config -A -T -L -g {genome} -t ../1.RNA_based/3.PASA/transcripts.fasta.clean -u ../1.RNA_based/3.PASA/transcripts.fasta --annots $id ;done 
        ll -tr |grep gff3$  |tail -n 1 |awk '{print $NF}' |while read id ;do Launch_PASA_pipeline.pl -c update.config -A -T -L -g {genome} -t ../1.RNA_based/3.PASA/transcripts.fasta.clean -u ../1.RNA_based/3.PASA/transcripts.fasta --annots $id ;done
        ll -tr  |grep gff3$  |tail -n 1 |awk '{print $NF}' |while read id ;do mv $id PASA_update.gff3;done 
        """



