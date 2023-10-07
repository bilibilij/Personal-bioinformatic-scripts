#!/usr/bin/bash

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file


bam= $(echo "${${snakemake_input[bam]}[*]}" ) 
samtool=${snakemake_input[samtools]}

bam=$(echo "${snakemake_input[bam]}")

$samtool merge -@ 8  ${snakemake_output[merge_bam]} $bam

#${snakemake_output[merge_bam]}

