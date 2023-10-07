![image](https://github.com/bilibilij/Personal-bioinformatic-scripts/blob/main/5.ATAC_pipe/ACR.jpg)
1. Peak calling
ATAC_pipe.pl is a pipeline which integrate the fastp, bowtie, macs2, samtools, picard software to perform peak calling.
2. ACR Peak filter
We filtered ATAC-peak by Tn5 frequencies using scripts of ACR_filter.
3. ACR classification
We combined gene collinearity information to classify the ACR of WGD dulicates into shared-accessible, shared-inaccessible and not shared.
ACR_distance-processing.pl => 1.Call_Homo_ACR.pl => 2.Private_but_Conserved.pl => 3.Final_check_with_Cshared.pl
