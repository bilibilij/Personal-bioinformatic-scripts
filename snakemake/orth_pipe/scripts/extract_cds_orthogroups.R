library(readr)
library(magrittr)
library(tidyverse)
library(dplyr)
library(ape)
library(Biostrings)
library(yaml)
options(stringsAsFactors = FALSE)

if (snakemake@params[["clans_cds"]]){
  species <- c(snakemake@params[["postWGD"]], snakemake@params[["preWGD"]])
} else {
  species <- names(snakemake@params[["inputFiles"]])
}
orthogroups <- readr::read_tsv(paste0(snakemake@input[["OG_table_with_version"]])) %>%
  select(OG,species) %>%
  gather(key="spc",value="geneID",-OG) %>%
  drop_na() %>% 
  mutate(geneID=strsplit(geneID,split = ", ")) %>% 
  unnest() 

if (snakemake@params[["splitBefore"]]){
  supportThreshold <- snakemake@params[["support"]]
  speciesTree <- read.tree(snakemake@input[["species_tree"]])
  speciesRoot <- speciesTree$edge[1,1]
  node <- getMRCA(speciesTree, species)
  splitNodes <- node
  while (node != speciesRoot){
    node <- speciesTree$edge[which(speciesTree$edge[,2]==node),1]
    splitNodes <- c(splitNodes, node)
  }
  splitNodes <- speciesTree$node.label[splitNodes - length(speciesTree$tip.label)]
  dup <- read_tsv(snakemake@input[["dup_events"]]) 
  dupOG <- split(dup,dup$Orthogroup)
  allClans <- lapply(dupOG, function(d){
    og = d$Orthogroup[1]
    tree <- read.tree(paste0(snakemake@input[["resolved_trees"]],"/",og,"_tree.txt"))
    og = gsub("OG",paste0("OG",snakemake@params[["version"]]),d$Orthogroup[1]) #Get OG name with version
    d <- d %>% mutate(`Gene Tree Node` = as.integer(gsub("n","",`Gene Tree Node`)) + length(tree$tip.label) +1) %>% #Replace node label by node id of gene tree
      filter(`Species Tree Node` %in% splitNodes, Support >= supportThreshold)                                      #keep only nodes with support >= threhsold
    if (nrow(d) > 0){
      children <- lapply(d$`Gene Tree Node`, function(i){       #get children of all nodes in d 
        tree$edge[which(tree$edge[,1] == i),2]
      }) 
      children <- do.call(rbind, children) %>% set_colnames(c("left","right"))
      d <- cbind(d, children)
      root <- tree$edge[1,1]
      anc <- sapply(d$`Gene Tree Node`,function(i){ #get all ancestors nodes of all nodes in d that are children of some other nodes.
          while (i!=root){
            if (i %in% c(d$left, d$right)) return(i)
            i <- tree$edge[which(tree$edge[,2]==i),1]
          }
          return (NA)
      }) 
      anc <- anc[!is.na(anc)]
      d <- mutate(d, left = ifelse(left %in% anc, NA, left),
               right = ifelse(right %in% anc, NA, right)) %>%
        filter(!is.na(left) | !is.na(right))               
      nodesToSplit <- c(d$left, d$right)
      nodesToSplit <-  nodesToSplit[!is.na(nodesToSplit)]
      clans <- lapply(nodesToSplit, function(node){
        if (node > length(tree$tip.label)) genes <- extract.clade(tree,node)$tip.label
        else genes <- tree$tip.label[node]
        sps <- str_extract(genes,"^[^_]+")
        genes[sps %in% species]
        if (length(genes)>0) {
          data.frame("OG" = paste0(og,".",node), spc = str_extract(genes,"^[^_]+"), geneID = gsub("^[^_]+_","",genes))
        } else NULL
      }) %>% bind_rows()
    } else {
      clans <- NULL
    }
    if (is.null(clans)){
      data.frame("OG" =og, spc = str_extract(tree$tip.label,"^[^_]+"), geneID = gsub("^[^_]+_","",tree$tip.label))
    } else {
      remain <- setdiff(tree$tip.label, paste0(clans$spc,"_",clans$geneID))
      remain <- remain[str_extract(remain,"^[^_]+") %in% species]
      if (length(remain)>0) {
        remain <- data.frame("OG" = paste0(og,".0"), spc = str_extract(remain,"^[^_]+"), geneID = gsub("^[^_]+_","",remain))
        results <- rbind(remain,clans)
      } else {
        clans
      }  
    }
  }) %>% bind_rows()
  missingOG <- filter(orthogroups, !(OG %in% gsub("OG",paste0("OG",snakemake@params[["version"]]),names(dupOG))))
  orthogroups <- rbind(allClans, missingOG)
}

all_cds <- lapply(snakemake@input[["cds_fasta"]],function(sp){
  cds <- readBStringSet(sp)
  spc = gsub(".fna","",basename(sp))
  data.frame("spc"=spc,"geneID"=names(cds),"Sequence"=cds)
}) %>% bind_rows()
all_cds <- left_join(all_cds, orthogroups)
split_orthogroup <- split(all_cds,all_cds$OG)

cds_outFolder_macse <- snakemake@output[["cds_orthogroup_macse"]]
cds_outFolder_mafft <- snakemake@output[["cds_orthogroup_mafft"]]
dir.create(cds_outFolder_macse, recursive = T)
dir.create(cds_outFolder_mafft, recursive = T)
for (i in 1:length(split_orthogroup)){
  seq <- BStringSet(split_orthogroup[[i]]$Sequence)
  names(seq) <- paste0(split_orthogroup[[i]]$geneID,"_",split_orthogroup[[i]]$spc)
  if (length(seq)>1){
    og <- split_orthogroup[[i]]$OG[1]
    if (length(seq) > snakemake@params[["max_nbSeq_macse"]] || max(lengths(seq)) > snakemake@params[["max_lenSeq_macse"]]){
      writeXStringSet(seq,filepath = paste0(cds_outFolder_mafft,"/",og,".fasta"))
    } else {
      writeXStringSet(seq,filepath = paste0(cds_outFolder_macse,"/",og,".fasta"))    
    }
  }
}
