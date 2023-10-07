library(dplyr)
library(stringr)
library(treeio)
library(magrittr)
library(readr)
options(stringsAsFactors = F)
OGtbl <- read_tsv(snakemake@input[["og_tbl"]], col_types = cols(.default = "c")) %>%
  select(-OG,-spc) %>% 
  mutate(geneID = as.character(geneID))
nc <- colnames(OGtbl)
nc <- nc[nc!="geneID"]
genes <- read.table(snakemake@input[["genes"]], sep = "\t", header =T) %>%
  select(-remapped_coordinate,-remapped,-contains("tandem")) %>%
  mutate(id = as.character(id))
multiplicons <- read.table(snakemake@input[["multiplicons"]], sep = "\t", header = T) 
if (snakemake@params[["level2"]]){
  multiplicons <- filter(multiplicons, level == 2)
}
anchorpoints <- read.table(snakemake@input[["anchorpoints"]], sep = "\t", header =T) %>%
  filter(multiplicon %in% multiplicons$id) %>% 
  select(-id,-basecluster,-is_real_anchorpoint,-coord_x,-coord_y) %>%
  mutate(gene_x = as.character(gene_x), gene_y = as.character(gene_y)) %>%
  left_join(set_colnames(genes, paste0(colnames(genes),"_x")), 
                          by = c("gene_x" = "id_x")) %>%
  left_join(set_colnames(genes, paste0(colnames(genes),"_y")), 
            by = c("gene_y" = "id_y")) %>%
  left_join(set_colnames(OGtbl, paste0(colnames(OGtbl),"_x")), 
            by = c("gene_x" = "geneID_x")) %>%
  left_join(set_colnames(OGtbl, paste0(colnames(OGtbl),"_y")), 
            by = c("gene_y" = "geneID_y")) %>%
  mutate("speciesTreeNode" = "OG", "node" = NA)
for (n in nc){
  id <- which(anchorpoints[[paste0(n,"_x")]] == anchorpoints[[paste0(n,"_y")]])
  anchorpoints$speciesTreeNode[id] <- n
  anchorpoints$node[id] <- anchorpoints[[paste0(n,"_x")]][id]
  anchorpoints <- select(anchorpoints, -contains(n, ignore.case = F))
}
anchorpoints <- mutate(anchorpoints, "OG" = gsub("[.][0-9]+$","",node), "geneTreeNode" = as.integer(str_extract(node,"[0-9]+$"))) %>%
  select(-node)
trees <- readRDS(snakemake@input[["rds_tree"]])
anchorpoints <- lapply(split(anchorpoints, anchorpoints$OG), function(an_og){
  og <- an_og$OG[1]
  tree <- trees[[og]]
  if (!is.null(tree) && "B" %in% names(tree@data)){
    support = as.integer(tree@data$B[an_og$geneTreeNode])  
  } else {
    support = NA
  }
  mutate(an_og, "support" = support)
}) %>% bind_rows() %>%
  arrange(multiplicon, coordinate_x, coordinate_y) %>%
  select(-OG, OG)
write_tsv(anchorpoints, snakemake@output[["ohnolog_tbl"]])
