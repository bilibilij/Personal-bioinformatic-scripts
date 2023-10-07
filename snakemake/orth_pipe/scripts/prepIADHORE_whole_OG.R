source("scripts/generateIADHORE_input.R")
species <- basename(snakemake@output[["geneListDir"]])
# load gene positions
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
geneListsDir <- set_names(snakemake@output[["geneListDir"]],species)
walk(geneListsDir, dir.create)
geneLists <- posTbls %>%
  arrange(spc,seqname,start) %>%
  group_by(spc,seqname) %>%
  do( geneList = paste0(.$geneID,.$strand)) %>%
  ungroup() %>%
  mutate( geneListFile = file.path(geneListsDir[spc],paste0(seqname,".lst")))
with(geneLists, pwalk(.l = list(text = geneList, con = geneListFile), .f = writeLines))

outPrepDir <- dirname(snakemake@output[["gene_family"]])
cat("Generating i-ADHoRE input files in",outPrepDir,"...\n")
geneFams <- read_tsv(snakemake@input[["orthofinder_og_table"]]) %>%
  select(OG,species) %>%
  gather(key="spc",value="geneID",-OG) %>%
  drop_na() %>% 
  mutate(geneID=strsplit(geneID,split = ", ")) %>% 
  unnest() %>%
  select(-spc)

# generate all input files
outResDir <- file.path(gsub("Input","Results",dirname(outPrepDir)), basename(outPrepDir))
generateIAdhoreInput(outPrepDir, outResDir, spcs = species, geneFams, posTbls,geneLists, lvl2only = F)

