library(tidyverse)
library(yaml)
library(ape)
library(treeio)

generateIAdhoreInput <- function(jobPrepOutDir, jobResOutDir, spcs, geneFams, posTbls,geneLists, lvl2only){
    # save geneFamily file after adding the genes with no family
    posTbls %>%
        filter(spc %in% spcs) %>%
        select( geneID ) %>%
        left_join( geneFams, by="geneID") %>%
        # genes with no OG will get a unique id based on the geneID
        mutate( OG = ifelse(is.na(OG),paste0("noOG_",geneID),OG)) %>%
        write_tsv(col_names = F, path = file.path(jobPrepOutDir,"gene_family.txt"))
    
    # generate settings.ini:
    writeLines( con = file.path(jobPrepOutDir,"settings.ini"),
                text = generateIAdhoreSettingsFile(geneLists = filter(geneLists, spc %in% spcs),
                                                   geneFamFile = file.path(jobPrepOutDir,"gene_family.txt"),
                                                   resDir = jobResOutDir,
                                                   lvl2only = lvl2only))
}

generateIAdhoreSettingsFile <- function(geneLists,geneFamFile,resDir,lvl2only=T,nThreads=1){
    paste(
        sep="\n",
        # add genelist file names for each genome
        geneLists %>% 
            group_by(spc) %>% 
            summarise( geneListsTxt = paste(seqname,geneListFile,collapse = "\n")) %>% 
            with( paste(paste0("genome=",spc),geneListsTxt,sep="\n",collapse = "\n\n")),
        "",
        paste0("blast_table=",geneFamFile),
        paste0("output_path=",resDir),
        "gap_size= 30",
        "cluster_gap= 50",
        "q_value=0.9",
        "prob_cutoff=0.01",
        "anchor_points=3",
        "alignment_method=gg2",
        paste0("level_2_only=",ifelse(lvl2only,"true","false")),
        "table_type=family",
        "multiple_hypothesis_correction=FDR",
        "visualizeGHM=false",
        "visualizeAlignment=false",
        paste0("number_of_threads=",nThreads))
}