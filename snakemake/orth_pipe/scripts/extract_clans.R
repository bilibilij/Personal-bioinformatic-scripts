library(treeio)
library(ape)
library(yaml)
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(Biostrings)
options(stringsAsFactors = F)

getDesInNode <- function(tree, node){
  if (node > length(tree$tip.label)){
    children <- tree$edge[which(tree$edge[,1] == node),2]
    return(c(node,Recall(tree,children[1]),Recall(tree,children[2])))
  } else{
    return(c())
  }
}

OG <- read_tsv(snakemake@input[["og_tbl"]]) 
cds_alignments <- paste0(dirname(snakemake@output[["clans_cds"]]),"/nt_aln_for_treebest")
dir.create(snakemake@output[["clans_cds"]])
dir.create(snakemake@output[["clans_tree_nwk"]])

supportThresholdClan <- snakemake@params[["supportClan"]]*100
allTrees <- readRDS(snakemake@input[["rds_tree"]])
speciesTree <- read.tree(snakemake@input[["species_tree"]])
minClan <- snakemake@params[["min"]]
clans <- c(snakemake@params[["preWGD"]], snakemake@params[["postWGD"]])
clanLabel <- speciesTree$node.label[getMRCA(speciesTree, clans) - length(speciesTree$tip.label)]
OGclans <- OG[!is.na(OG[[clanLabel]]),] %>%
  select(geneID, spc, OG, clanLabel) %>%
  set_colnames(c("geneID","spc","OG","splitNode")) %>%
  mutate(splitNode = as.integer(str_extract(splitNode,"[0-9]+$")))
splitOG <- split(OGclans, OGclans$OG)
allClanTrees <- lapply(splitOG, function(sog){
  og <- sog$OG[1]
  tree <- allTrees[[og]]
  phy <- tree@phylo
  cds <- readBStringSet(paste0(cds_alignments,"/",og,".fasta"))
  d <- data.frame("seq_id"=1:length(cds), "geneID" = names(cds))
  if ("B" %in% names(tree@data)){
    treeData <- select(tree@data, D, B, node) %>%
      mutate(B = as.integer(B))
    phy$node.label <- as.character(treeData$B[(length(phy$tip.label)+1):(nrow(tree@data))])  
    lapply(split(sog, sog$splitNode), function(ssog){
      splitNode = ssog$splitNode[1]
      bs <- filter(treeData,node == splitNode)$B[1]
      if (nrow(ssog) > minClan && bs >= supportThresholdClan){
        clans <- extract.clade(phy, splitNode)
        ogClans <- paste0(og,".",splitNode)
        ogSplit <- filter(d, geneID %in% clans$tip.label)
        cds_sub <- cds[ogSplit$seq_id]
        writeXStringSet(cds_sub,filepath = paste0(snakemake@output[["clans_cds"]],"/",ogClans,".fasta"))
        write.tree(clans, file = paste0(snakemake@output[["clans_tree_nwk"]],"/",ogClans,".nwk"))
        clans
      } else {
        NULL
      } 
    })
  } else {
    lapply(split(sog, sog$splitNode), function(ssog){
      splitNode = ssog$splitNode[1]
      if (nrow(ssog) > minClan){
        clans <- extract.clade(phy, splitNode)
        ogClans <- paste0(og,".",splitNode)
        ogSplit <- filter(d, geneID %in% clans$tip.label)
        cds_sub <- cds[ogSplit$seq_id]
        writeXStringSet(cds_sub,filepath = paste0(snakemake@output[["clans_cds"]],"/",ogClans,".fasta"))
        write.tree(clans, file = paste0(snakemake@output[["clans_tree_nwk"]],"/",ogClans,".nwk"))
        clans
      } else {
        NULL
      }
    })
  } 
}) %>% unlist(recursive = F)
allClanTrees <- allClanTrees[!sapply(allClanTrees,is.null)]
write_rds(allClanTrees,paste0(dirname(snakemake@output[["clans_tree_nwk"]]),"/allClans.RDS"))


