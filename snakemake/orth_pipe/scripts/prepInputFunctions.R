
filterRedundantGenePos <- function(genePosTbl,redundantSaveFile){
  if(!any(duplicated(gffTbls$genePosTbl$geneID))){
    # Don't do anything if there are no redundant gene positions
    return(genePosTbl)
  }
  cat("Duplicated gene positions! Saving alternative positions to:",redundantSaveFile,"\n")
  
  # figure out which of the duplicated positions to use
  redundantGFF_geneIDs <-
    genePosTbl %>% 
    filter(geneID %in% geneID[duplicated(geneID)]) %>% 
    # prioritize by gene length
    mutate( length=end-start) %>% 
    # prioritize genes on chromosomes before scaffolds
    mutate( length=length+1e6*grepl("NC_",seqname)) %>%
    group_by(geneID) %>%
    filter( (1:n() != match(max(length),length) ) ) %>% 
    with( gff_geneID)
  
  # store the redundant gene positions just in case
  genePosTbl %>% 
    filter( gff_geneID %in% redundantGFF_geneIDs) %>% 
    write_tsv( redundantSaveFile )
  
  # Return the rest
  return( filter(genePosTbl, !(gff_geneID %in% redundantGFF_geneIDs) ) )
}

# add protein length to ID conversion table and mark longest protein per gene
findLongestProtein <- function(IDtbl, seqLengths){
  IDtbl %>% 
    mutate( proteinLength = seqLengths[proteinID]  ) %>% 
    group_by(geneID) %>% 
    # find the first longest protein per gene
    arrange(-proteinLength) %>% 
    mutate( is_longest = (1:n())==1) 
}

# Check if IDs are matching
checkIDmissing <- function(inFasta,inGFF){
  if( any(duplicated(inFasta)) ) cat("Duplicated IDs in fasta file!")
  cat( "  Sequences in fasta file:",length(inFasta),"\n")
  cat( "  Fasta seq IDs not found in GFF:",sum(!(inFasta %in% inGFF)),"\n")
  cat( "  GFF IDs missing in fasta file:",sum(!inGFF %in% inFasta),"\n")
}



# NCBI fasta headers:
#   Protein: ">XP_019904295.1 ..."
#   CDS: ">lcl|NC_025975.3_cds_XP_019904295.1_14772 ..."
convertIdFun_CDS_NCBI = function(seqID){sub("lcl|.*_cds_(.*)_.*","\\1",seqID)}

# remove the version number suffix
convertIdFun_ENSEMBL = function(seqID){sub("\\.[0-9]+$","",seqID)}
