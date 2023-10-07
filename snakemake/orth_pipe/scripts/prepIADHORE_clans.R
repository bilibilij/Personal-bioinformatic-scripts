library(ape)
library(magrittr)
source("scripts/generateIADHORE_input.R")
species <- basename(snakemake@input[["geneListDir"]])
posTbls <- 
  snakemake@input[["genePos"]] %>% 
  set_names(species) %>% 
  map_df(read_tsv, .id="spc", col_types = cols(
    seqname = col_character(),
    start = col_double(),
    end = col_double(),
    strand = col_character(),
    geneID = col_character()
  )) 
# save gene lists
geneListsDir <- set_names(snakemake@input[["geneListDir"]],species)
geneLists <-
  posTbls %>%
  arrange(spc,seqname,start) %>%
  group_by(spc,seqname) %>%
  do( geneList = paste0(.$geneID,.$strand)) %>%
  ungroup() %>%
  mutate( geneListFile = file.path(geneListsDir[spc],paste0(seqname,".lst")))

outPrepDir <- dirname(snakemake@output[["gene_family"]])
dir.create(outPrepDir,recursive = T)
species_tree <- read.tree(snakemake@input[["species_tree"]])
nodeId <- getMRCA(species_tree, species)
ancestralNode <- species_tree$node.label[nodeId - length(species_tree$tip.label)]
cat("Generating i-ADHoRE input files in",outPrepDir,"...\n")

geneFams <- read_tsv(snakemake@input[["og_tbl"]], col_types = cols(.default = "c")) %>%
  filter(spc %in% species, !is.na(ancestralNode)) %>%
  select(ancestralNode, geneID) %>%
  set_colnames(c("OG","geneID"))

outResDir <- file.path(gsub("Input$","Results",dirname(outPrepDir)), basename(outPrepDir))
generateIAdhoreInput(outPrepDir, outResDir, spcs = species, geneFams, posTbls,geneLists, lvl2only = F)
