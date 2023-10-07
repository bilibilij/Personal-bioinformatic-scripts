library(yaml)
source("scripts/getOrthoSubGroups.R")
treesRDSfile = snakemake@input[["rds_tree"]]
spcTreeFile = snakemake@input[["species_tree"]]
OGtblFileOut = snakemake@output[["og_tbl_file_out"]]

geneTrees <- readRDS(treesRDSfile)
spcTree <- read.tree(spcTreeFile)

splitOG <-
  imap_dfr(geneTrees, ~ splitOrthoGroups(spcTree, treedata = .x, OGid = .y)) %>%
  mutate(
    spc = sub(".*_([^_]+)$", "\\1", geneID),
    geneID = sub("_[^_]+$", "", geneID)
  ) %>%
  arrange(OG, spc) %>%
  write_tsv(path = OGtblFileOut)
