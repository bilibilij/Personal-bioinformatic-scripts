library(readr)
library(dplyr)
library(Biostrings)
library(ape)
options(stringsAsFactors = F)

files_dup <- list.files(snakemake@input[["dup_clans_tree_nwk"]],full.names = T)
files_nondup <- list.files(snakemake@input[["nondup_clans_tree_nwk"]],full.names = T)
allFiles <- rbind(data.frame("file" = files_dup) %>% mutate("dup" = TRUE),
               data.frame("file" = files_nondup) %>% mutate("dup" = FALSE)) %>%
  mutate(OG = gsub("[.][0-9]+.nwk","",basename(file)))
dir.create(snakemake@output[["dup_clans_cds"]])
dir.create(snakemake@output[["nondup_clans_cds"]])
cds_alignments <- paste0(dirname(snakemake@output[["dup_clans_cds"]]),"/nt_aln_for_treebest")
for(og in unique(allFiles$OG)){
  cds <- readBStringSet(paste0(cds_alignments,"/",og,".fasta"))
  if (length(cds)>0){
    d <- data.frame("seq_id"=1:length(cds), "geneID" = names(cds))
    trees_dup <- filter(allFiles, OG == og, dup == TRUE)
    for (t in trees_dup$file){
      ogClans <- gsub('.nwk','',basename(t))
      ogSplit <- filter(d, geneID %in% read.tree(t)$tip.label)
      cds_sub <- cds[ogSplit$seq_id]
      writeXStringSet(cds_sub,filepath = paste0(snakemake@output[["dup_clans_cds"]],"/",ogClans,".fasta"))
    }
    trees_nondup <- filter(allFiles, OG == og, dup == FALSE)
    for (t in trees_nondup$file){
      ogClans <- gsub('.nwk','',basename(t))
      ogSplit <- filter(d, geneID %in% read.tree(t)$tip.label)
      cds_sub <- cds[ogSplit$seq_id]
      writeXStringSet(cds_sub,filepath = paste0(snakemake@output[["nondup_clans_cds"]],"/",ogClans,".fasta"))
    }
  }
}
