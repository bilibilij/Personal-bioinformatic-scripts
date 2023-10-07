# generate small test dataset
library(tidyverse)

OGfile <- "/mnt/SCRATCH/thut/pipeline_results_all_data/Ortho_pipeline/OrthoFinder/Orthogroups/Orthogroups.tsv"
allDataDir <- "/mnt/SCRATCH/thut/pipeline_results_all_data/PrepInput_pipeline"

# define and create output directories
outBaseDir <- "/mnt/SCRATCH/thut/testData/PrepInput_pipeline"
outDirs <- sapply(c("proteinFastas","cdsFastas","genePosTbls","geneIDtbls"), function(.x)(file.path(outBaseDir,.x)),simplify = F)
walk(outDirs, dir.create, recursive=T)

# function to extract protein and CDS sequences for a set of genes in a species
extractSubset <- function( geneIDs, spc){
  AAseqs <-
    file.path(allDataDir,"proteinFastas",paste0(spc,".faa")) %>% 
    seqinr::read.fasta(as.string = T,seqtype="AA")
  
  seqinr::write.fasta(sequences = AAseqs[geneIDs], names=geneIDs,
                      file.out = file.path(outDirs$proteinFastas,paste0(spc,".faa")) )
  
  CDSseqs <-
    file.path(allDataDir,"cdsFastas",paste0(spc,".fna")) %>% 
    seqinr::read.fasta(as.string = T,seqtype="DNA")
      
  seqinr::write.fasta(sequences = CDSseqs[geneIDs], names=geneIDs,
                      file.out = file.path(outDirs$cdsFastas,paste0(spc,".fna")) )
  read.table(file.path(allDataDir,"genePosTbls",paste0(spc,"_genePos.tsv")), header = TRUE) %>%
      mutate(geneID = as.character(geneID)) %>%
      filter(geneID %in% geneIDs) %>%
      write.table(file = file.path(outDirs$genePosTbls,paste0(spc,"_genePos.tsv")),quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
  
  read.table(file.path(allDataDir,"geneIDtbls",paste0(spc,"_IDtbl.tsv")), header = TRUE) %>%
    mutate(geneID = as.character(geneID)) %>%
    filter(geneID %in% geneIDs) %>%
    write.table(file = file.path(outDirs$geneIDtbls,paste0(spc,"_IDtbl.tsv")),quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
  
}



spcs <- c("Drer","Eluc","Ssal","Omyk")

OG <- 
  read_tsv(OGfile,col_types = cols(.default=col_character())) %>% 
  dplyr::rename(OG=Orthogroup) %>% 
  select(OG, spcs) %>% 
  gather(key="spc",value="geneID",-OG) %>% 
  drop_na() %>% 
  mutate(geneID=strsplit(geneID,split = ", ")) %>% 
  unnest()

Eluc_genePos <- read_tsv(file.path(allDataDir,"genePosTbls/Eluc_genePos.tsv"),col_types = "ciiccc")

Eluc_genePos %>% 
    # select a single chromosome in Eluc
    filter(seqname == "NC_025969.3", end <= 1e7) %>% 
    select(geneID) %>% 
    # get corresponding OG
    inner_join(OG,by="geneID") %>% 
    # add OG size and remove large OGs
    left_join( dplyr::count(OG,OG), by="OG") %>% 
    filter(n>12) %>% 
    # get geneIDs for selected OGs
    select(OG) %>% 
    unique() %>%
    left_join(OG,by="OG") %>% 
    split(x = .$geneID,f = .$spc) %>%
    iwalk( ~ extractSubset( geneIDs=.x,spc=.y))

