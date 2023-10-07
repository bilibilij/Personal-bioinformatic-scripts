library(dplyr)
library(magrittr)
orthogroups <- readr::read_tsv(paste0(snakemake@input[["orthofinder_og_table"]])) 
names(orthogroups)[1] <- "OG"
orthogroups <- mutate(orthogroups, OG=paste0("OG",snakemake@params[["version"]],gsub("OG","",OG)))
write.table(orthogroups,snakemake@output[["OG_table_with_version"]], quote = FALSE,col.names = TRUE, row.names = FALSE, sep = "\t")
