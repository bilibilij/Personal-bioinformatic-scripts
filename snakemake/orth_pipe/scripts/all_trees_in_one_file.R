library(tidyverse)
library(treeio)

saveAllTreesInOneFile <- function(geneTreesDir,outFile){
    currentDir <- getwd()
    setwd(geneTreesDir)
    system(paste0("ls OG*nhx | xargs cat  >> ",outFile))
    allTree <- readLines(outFile)
    ends <- grepl(";$",allTree)
    newLine <- rep("",length(allTree))
    newLine[ends] <- "\n"
    newTree <- paste0(allTree,newLine,collapse = "")
    cat(newTree,file = outFile)
    newTree <- readLines(outFile)
    ogs <- gsub(".nhx","",list.files(pattern = ".nhx"))
    data = data.frame("og"=ogs,"cdsTree"=newTree)
    write_tsv(data,outFile)
    setwd(currentDir)
}

saveTheTrees <- function(treeDir, treesRDSfile){
    allTrees <- 
        dir(treeDir,full.names = T) %>% 
        set_names(.,sub(".nhx","",basename(.))) %>% 
        map(read.nhx)
    
    saveRDS(allTrees,file = treesRDSfile)
}

geneTrees <- sort(unlist(snakemake@input))
outFile <- snakemake@output[["nhx_tree"]]
file.create(outFile)
saveAllTreesInOneFile(dirname(geneTrees[1]),outFile)
saveTheTrees(treeDir = dirname(snakemake@input[[1]]), treesRDSfile = snakemake@output[["rds_tree"]])