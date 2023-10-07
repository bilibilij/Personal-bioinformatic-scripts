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
cds_alignments <- paste0(dirname(snakemake@output[["dup_clans_cds"]]),"/nt_aln_for_treebest")
dir.create(snakemake@output[["dup_clans_cds"]])
dir.create(snakemake@output[["nondup_clans_cds"]])
dir.create(snakemake@output[["dup_clans_tree_nwk"]])
dir.create(snakemake@output[["nondup_clans_tree_nwk"]])

supportThresholdClan <- snakemake@params[["supportClan"]]*100
supportThresholdDup <- snakemake@params[["supportDup"]]*100
allTrees <- readRDS(snakemake@input[["rds_tree"]])
speciesTree <- read.tree(snakemake@input[["species_tree"]])
minClan <- snakemake@params[["min"]]
clans <- c(snakemake@params[["preWGD"]], snakemake@params[["postWGD"]])
clanLabel <- speciesTree$node.label[getMRCA(speciesTree, clans) - length(speciesTree$tip.label)]
WGDLabel <- speciesTree$node.label[getMRCA(speciesTree, snakemake@params[["postWGD"]]) - length(speciesTree$tip.label)]

OGclans <- OG[!is.na(OG[[clanLabel]]),] %>%
  select(geneID, spc, OG, clanLabel, WGDLabel) %>%
  set_colnames(c("geneID","spc","OG","splitNode","WGDNode")) %>%
  mutate(splitNode = as.integer(str_extract(splitNode,"[0-9]+$")),
         WGDNode = as.integer(str_extract(WGDNode,"[0-9]+$")))
splitOG <- split(OGclans, OGclans$OG)
allClanTrees <- lapply(splitOG, function(sog){
  og <- sog$OG[1]
  tree <- allTrees[[og]]
  phy <- tree@phylo
  if ("B" %in% names(tree@data)){
    treeData <- select(tree@data, D, B, node) %>%
      mutate(B = as.integer(B))
    phy$node.label <- as.character(treeData$B[(length(phy$tip.label)+1):(nrow(tree@data))])  
    cds <- readBStringSet(paste0(cds_alignments,"/",og,".fasta"))
    d <- data.frame("seq_id"=1:length(cds), "geneID" = names(cds))
    lapply(split(sog, sog$splitNode), function(ssog){
      splitNode = ssog$splitNode[1]
      bs <- filter(treeData,node == splitNode)$B[1]
      if (nrow(ssog) > minClan && bs >= supportThresholdClan){
        clantree <- extract.clade(phy, splitNode)
        clantree$dup <- FALSE
        clanNodes <- unique(ssog$WGDNode)
        clanNodes <- clanNodes[!is.na(clanNodes)]
        if (length(clanNodes) == 1){ #the initial WGD node was lost
          clanNodes <- phy$edge[which(phy$edge[,1] == clanNodes),2]
        }
        if (length(clanNodes) == 2){
          dupNode <- phy$edge[phy$edge[,2]==clanNodes[1],1]
          bd <- filter(treeData,node == dupNode)$B[1]
          if (filter(tree@data, node == dupNode)$D[1]=="Y" && (bd >= supportThresholdDup || is.na(bd))){
              if (clanNodes[1] > length(phy$tip.label) && clanNodes[2] > length(phy$tip.label)){
                inNodeClan1 <- getDesInNode(phy, clanNodes[1])
                inNodeClan2 <- getDesInNode(phy, clanNodes[2])
                inDupClan1 <- filter(tree@data, node %in% inNodeClan1, D=="Y", !(S %in% clans))
                inDupClan2 <- filter(tree@data, node %in% inNodeClan2, D=="Y", !(S %in% clans))
                spcClan1 <- unique(str_extract(extract.clade(phy,clanNodes[1])$tip.label,"[^_]+$"))
                spcClan2 <- unique(str_extract(extract.clade(phy,clanNodes[2])$tip.label,"[^_]+$"))
                if (nrow(inDupClan1)==0 && nrow(inDupClan2)==0 && length(intersect(spcClan1,spcClan2))>2){
                  clantree$dup <- TRUE
                }
              }
          }
        }
        ogClans <- paste0(og,".",splitNode)
        ogSplit <- filter(d, geneID %in% clantree$tip.label)
        cds_sub <- cds[ogSplit$seq_id]
        if (clantree$dup){
          writeXStringSet(cds_sub,filepath = paste0(snakemake@output[["dup_clans_cds"]],"/",ogClans,".fasta"))
          write.tree(clantree, file = paste0(snakemake@output[["dup_clans_tree_nwk"]],"/",ogClans,".nwk"))  
        } else {
          writeXStringSet(cds_sub,filepath = paste0(snakemake@output[["nondup_clans_cds"]],"/",ogClans,".fasta"))
          write.tree(clantree, file = paste0(snakemake@output[["nondup_clans_tree_nwk"]],"/",ogClans,".nwk"))
        }
        clantree
      } else {
        NULL
      } 
    })
  } else {
    lapply(split(sog, sog$splitNode), function(ssog){
      splitNode = ssog$splitNode[1]
      if (nrow(ssog) > minClan){
        clantree <- extract.clade(phy, splitNode)
        clantree$dup <- FALSE
        ogClans <- paste0(og,".",splitNode)
        ogSplit <- filter(d, geneID %in% clantree$tip.label)
        cds_sub <- cds[ogSplit$seq_id]
        writeXStringSet(cds_sub,filepath = paste0(snakemake@output[["nondup_clans_cds"]],"/",ogClans,".fasta"))
        write.tree(clantree, file = paste0(snakemake@output[["nondup_clans_tree_nwk"]],"/",ogClans,".nwk"))
        clantree
      } else {
        NULL
      }
    })
  } 
}) %>% unlist(recursive = F)
allClanTrees <- allClanTrees[!sapply(allClanTrees,is.null)]
write_rds(allClanTrees,paste0(dirname(snakemake@output[["dup_clans_tree_nwk"]]),"/allClans.RDS"))
