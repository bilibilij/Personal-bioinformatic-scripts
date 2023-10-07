
import os

configfile: "config.yaml"

MINIMAP=config["MINIMAP"]
SAMTOOLS=config["SAMTOOLS"]
ALL_SAMPLE = [key for key, value in config["inputFiles"].items()]


bam=config["bam"]
j=','
bam_str=j.join(bam)

genome=config["genome"]
gff=config["gff"]
gtf=config["gtf"]
#gtf=config["gtf"]
#RSCRIPT=config["Rscript"] 
WGD=config["WGD"]

rule all:
    input:
        #["1.minimap_reference_aln/"+str(sample)+".sorted.bam" for sample in ALL_SAMPLE],
        ["0.data_preprocess/"+str(sample)+".cleaned.fastq" for sample in ALL_SAMPLE],
        ["1.minimap_reference_aln/"+str(sample)+".flair.bed" for sample in ALL_SAMPLE],
        ["1.minimap_reference_aln/"+str(sample)+".flair.bam" for sample in ALL_SAMPLE],
        "1.minimap_reference_aln/ONT_flair_merge.bed",
        "1.minimap_reference_aln/ONT_flair_merge.bam",
        "2.correct/ONT_all_corrected.bed",
        "2.correct/flair_sam_junctions.bed",
        "3.collapse/flair.collapse.isoforms.bed",
        ["4.liqa/"+str(sample)+".isoform_expression_estimates.txt" for sample in ALL_SAMPLE],
        "4.liqa/flair.collapse.isoforms.refgene"
#        "3.2.quantity/flair_quantify.counts.tsv",
#        "3.3.splice/diffsplice.alt5.events.quant.tsv",





rule preprocess:
    input:
        gff,
        genome
    output:
        #bed=gff +".bed",
        bed = "0.data_preprocess/"+ os.path.split(gff)[1] + ".bed",
        #Chr=genome+".Chr"
        Chr=  "0.data_preprocess/"+ os.path.split(genome)[1] +".Chr"
    conda:
        "envs/AGAT.yaml"
    shell:
        """
        export PERL5LIB=/
        grep \> {genome} |sed 's/>//g' |grep -v scaf |grep -v fold |grep -v tig |  sed 's/\tOri.*//g' | perl   scripts/extract_fa_from_id.pl {genome} - > {output.Chr}
        agat_convert_sp_gff2bed.pl --gff {gff}  -o {output.bed}
        """

rule trim:
    input:
        fq = lambda w: config["inputFiles"][w.sample],
    output:
        trim_fq = "0.data_preprocess/{sample}.cleaned.fastq"
    conda:
        "envs/pychopper.yaml"
    threads:
        20
    shell:
        """
        pychopper -r 0.data_preprocess/{wildcards.sample}.report.pdf -u 0.data_preprocess/{wildcards.sample}.unclassified.fq -w 0.data_preprocess/{wildcards.sample}.rescued.fq -t {threads}  {input.fq} {output.trim_fq}
       """


rule flair_aln:
    input:
        fq = "0.data_preprocess/{sample}.cleaned.fastq",
        Chr=  "0.data_preprocess/"+ os.path.split(genome)[1] +".Chr",
        bed = "0.data_preprocess/"+ os.path.split(gff)[1] + ".bed"
    output:
        bam= "1.minimap_reference_aln/{sample}.flair.bam",
        bed= "1.minimap_reference_aln/{sample}.flair.bed"
    conda:
        "envs/flair.yaml"
    threads:
        30 
    shell:
        """
        flair align -g {input.Chr} -r {input.fq} --threads {threads} --junction_bed {input.bed}  -o 1.minimap_reference_aln/{wildcards.sample}.flair
        """

#rule flair_aln:
#    input:
#        fq= lambda w: config["inputFiles"][w.sample],
#    output:
#        bam=fq
    
#    shell:
#        """
#        {MINIMAP}
#        """


rule merge_bam:
    input:
        bam=["1.minimap_reference_aln/"+str(sample)+".flair.bam" for sample in ALL_SAMPLE],
        samtools=SAMTOOLS
    output:
        merge_bam= "1.minimap_reference_aln/ONT_flair_merge.bam"
    log:
        "1.minimap_reference_aln/merge_bam.log"
    script:
        "scripts/merge_bam.sh"

rule merge_bed:
    input:
        ["1.minimap_reference_aln/"+str(sample)+".flair.bed" for sample in ALL_SAMPLE],
    output:
        merge_bed= "1.minimap_reference_aln/ONT_flair_merge.bed"
    shell:
        """
        perl scripts/merge_bed.pl 1.minimap_reference_aln > {output.merge_bed} 
        """

rule junc:
    input:
        bam
    output:
        "2.correct/flair_sam_junctions.bed"
    shell:
        """
        junctions_from_sam -s bam_str -n flair 
        """

rule correct:
    input:
        bed= "1.minimap_reference_aln/ONT_flair_merge.bed",
        merge_bam= "1.minimap_reference_aln/ONT_flair_merge.bam",

        junc= "2.correct/flair_sam_junctions.bed"
    output:
        all_corr = "2.correct/ONT_all_corrected.bed"
    conda:
        "envs/flair.yaml"
    log:
        "2.correct/1_correct.log"
    shell:
        """
        flair correct -q {input.bed} -f {gtf}  --shortread {input.junc} -g {genome} --output 2.correct/ONT --threads {threads} 2>{log}
        """

k=","
fq=k.join(["0.data_preprocess/"+str(sample)+".cleaned.fastq" for sample in ALL_SAMPLE])

rule collapse:
    input:
        genome=genome,
        corrected_bed= "2.correct/ONT_all_corrected.bed",
        ONT= ["0.data_preprocess/"+str(sample)+".cleaned.fastq" for sample in ALL_SAMPLE],
        gtf=gtf ,
    output:
        "3.collapse/flair.collapse.isoforms.bed",
        "3.collapse/flair.collapse.isoforms.gtf",
        "3.collapse/flair.collapse.isoforms.fa"
    threads:
        60
    conda:
        "envs/flair.yaml"
    shell:
        """
        flair collapse -g {input.genome} -q {input.corrected_bed} -r {fq}  -o 3.collapse/flair.collapse -t {threads} -f {input.gtf} 

        """


with open("3.collapse/sample.txt", 'w') as f:
    for sample in ALL_SAMPLE:
        f.write(str(sample)+"\tcondition1\tbatch"+str(sample)+"\t"+"0.data_preprocess/"+str(sample)+".cleaned.fastq\n")



rule quantity:
    input:
        "3.collapse/flair.collapse.isoforms.fa",
        "3.collapse/sample.txt"
    output:
        #"3.2.quantity/flair_quantify.counts_matrix.tsv"
        "3.2.quantity/flair_quantify.counts.tsv"
    threads:
        100
    conda:
        "envs/flair.yaml"
    shell:
        """
        flair quantify  -r 3.collapse/sample.txt -i 3.collapse/flair.collapse.isoforms.fa  --output  3.2.quantity/flair_quantify --threads {threads} --sample_id_only --trust_ends  --stringent --check_splice --isoform_bed 3.collapse/flair.collapse.isoforms.bed

        """

rule splice:
    input:
        "3.2.quantity/flair_quantify.counts.tsv",
        "3.collapse/flair.collapse.isoforms.bed",
    output:
        "3.3.splice/diffsplice.alt5.events.quant.tsv"
    threads:
        40
    conda:
        "envs/flair.yaml"
    shell:
        """

        flair diffSplice  -i 3.collapse/flair.collapse.isoforms.bed -q 3.2.quantity/flair_quantify.counts.tsv --out_dir 3.3.splice --threads {threads}
        """

#rule specifc_exon:
#    input:
#        "3.2.quantity/flair_quantify.counts.tsv"
#    output:
#        "4.exon/exon.com.complete"
#    threads:
#        60
#    shell:
#        """
#        perl scripts/specifc_exon_in_wgd.pl {WGD}  3.collapse/flair.collapse.isoforms_filtered.gtf {genome} 4.exon/exon_blast 4.exon/exon.com {threads} re 
#
#        """



rule liqa_convert:
    input:
        gtf="3.collapse/flair.collapse.isoforms.gtf"
    output:
        refgene="4.liqa/flair.collapse.isoforms.refgene"
    conda:
        "envs/liqa.yaml"
    shell:
        """
        liqa -task refgene -format gtf -ref  {input.gtf}  -out  {output.refgene}  
        """





