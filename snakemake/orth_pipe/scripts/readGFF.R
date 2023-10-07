# generate geneID conversion table and gene position table from gff



### GFF format columns:
#
# seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
# source - name of the program that generated this feature, or the data source (database or project name)
# feature - feature type name, e.g. Gene, Variation, Similarity
# start - Start position of the feature, with sequence numbering starting at 1.
# end - End position of the feature, with sequence numbering starting at 1.
# score - A floating point value.
# strand - defined as + (forward) or - (reverse).
# frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
# attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

# Split and separate GFF attribute tag-value pairs into columns
gffSplitAttributes <- function(tbl){
  tbl %>% 
    # Create a column with unique values outside the attribute column
    # or else the spread will not work.
    mutate(gffSplitAttributes_tmpidx = 1:n()) %>% 
    mutate(attribute = strsplit(attribute,split = ";(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)",perl = T)) %>% 
    unnest(attribute) %>%
    mutate(key = str_extract(attribute,"^[^=]+")) %>%
    mutate(attribute = gsub("^[^=]+=","",attribute)) %>%
    spread(key = key,value = attribute) %>% 
    select(-gffSplitAttributes_tmpidx) 
}

convertGFF_ncbi <- function(gff){
  gffTbl <-
    read_tsv(gff,col_names = c("seqname","feature","start","end","strand","attribute"),comment="#", col_types="c-cii-c-c")
  
  # get transcriptID, geneID, gff_geneID, rnaID(to link with proteinID) from mRNA features
  mRNAs <-
    gffTbl %>% 
    filter(feature == "mRNA") %>% 
    select(attribute) %>%
    gffSplitAttributes() %>%
    select(Dbxref,ID,partial,Parent) %>% 
    # split and separate the Dbxref tag-value pairs
    mutate(Dbxref = strsplit(Dbxref,split = ",")) %>% 
    unnest() %>%
    separate(Dbxref,c("key","value"), sep=":") %>% 
    spread(key = key,value = value) %>% 
    dplyr::rename( rnaID = ID, transcriptID=Genbank, gff_geneID=Parent, geneID=GeneID)
  
  # get proteinID and parentID(to link with transcriptID or geneID) from CDS features
  CDSs <-
    gffTbl %>% 
    filter(feature == "CDS") %>% 
    select(attribute) %>% 
    # add temporary index
    mutate(idx = 1:n()) %>% 
    # split the attributes
    mutate(attribute = strsplit(attribute,split = ";")) %>% 
    unnest() %>% 
    separate(attribute,c("key","value"), sep="=") %>% 
    filter( key %in% c("Parent", "Name","Dbxref")) %>% 
    spread(key = key,value = value) %>%
    transmute( cdsParent=Parent, proteinID=Name) %>% 
    # Ther may be some CDSs with no proteinID. Don't include them.
    filter(!is.na(proteinID)) %>% 
    distinct()
  
  geneTbl <- 
    gffTbl %>% 
    filter(feature == "gene") %>% 
    gffSplitAttributes() %>%
    # filter(ID %in% IDtbl$gff_geneID) %>% 
    select(seqname,start,end,strand,ID,Dbxref) %>% 
    # split and separate the Dbxref tag-value pairs
    mutate(Dbxref = strsplit(Dbxref,split = ",")) %>% 
    unnest() %>%
    separate(Dbxref,c("key","value"), sep=":") %>% 
    spread(key = key,value = value) %>% 
    dplyr::rename( gff_geneID = ID, geneID=GeneID)
  
  # Merge protein and transcript ID tables
  IDtbl <- 
    bind_rows(
      inner_join(CDSs,mRNAs,by=c("cdsParent"="rnaID")), 
      # There may be some proteins without transcript IDs.. include them
      inner_join(CDSs,select(geneTbl,gff_geneID,geneID),by=c("cdsParent"="gff_geneID")) %>% 
        mutate(gff_geneID=cdsParent) 
    )
  # note that there are some CDSs that do not map to mRNAs (geneID=NA)
  # and some mRNAs that dont have CDSs (not protein coding, not included)
  
  
  return( 
    list(
      genePosTbl=filter(geneTbl,gff_geneID %in% IDtbl$gff_geneID),
      IDtbl=IDtbl %>% select( geneID, proteinID, transcriptID) %>% distinct()
    )
  )  
}

convertGFF_ensembl <- function(gff){
  gffTbl <-
    read_tsv(gff,col_names = c("seqname","feature","start","end","strand","attribute"),comment="#", col_types="c-cii-c-c")
  
  # get transcriptID and geneID from mRNA features
  mRNAs <-
    gffTbl %>% 
    filter(grepl("^ID=transcript:",attribute)) %>% 
    select(attribute) %>%
    gffSplitAttributes() %>%
    transmute( geneID=sub("gene:","",Parent), transcriptID=transcript_id)
  
  # get proteinID and parentID(to link with transcriptID or geneID) from CDS features
  CDSs <-
    gffTbl %>% 
    filter(grepl("^ID=CDS:",attribute)) %>% 
    select(attribute) %>%
    gffSplitAttributes() %>% 
    transmute( transcriptID=sub("transcript:","",Parent), proteinID=protein_id) %>% 
    distinct()
  
  
  geneTbl <- 
    gffTbl %>% 
    filter(grepl("^ID=gene:",attribute)) %>% 
    gffSplitAttributes() %>% 
    dplyr::select(seqname,start,end,strand,gene_id) %>% 
    dplyr::rename(geneID=gene_id)
  
  
  # Merge protein and transcript ID tables
  IDtbl <- left_join(CDSs,mRNAs,by="transcriptID")
  # note that there are some CDSs that do not map to mRNAs (geneID=NA)
  
  return( 
    list(
      genePosTbl=filter(geneTbl,geneID %in% IDtbl$geneID),
      IDtbl=IDtbl %>% select( geneID, proteinID, transcriptID)
    )
  )  
}

convertGFF_local <- function(gff){
  gffTbl <-
    read_tsv(gff,col_names = c("seqname","feature","start","end","strand","attribute"),comment="#", col_types="c-cii-c-c")
  
  # get transcriptID and geneID from mRNA features
  mRNAs <-
    gffTbl %>% 
    filter(feature == "mRNA") %>% 
    select(attribute) %>%
    gffSplitAttributes() %>%
    transmute( geneID=Parent, transcriptID=ID)
  
  geneTbl <- 
    gffTbl %>% 
    filter(feature == "gene") %>% 
    gffSplitAttributes() %>%
    select(seqname,start,end,strand,ID) %>% 
    dplyr::rename(geneID=ID)
  
  
  # protein IDs are same as and transcript IDs
  IDtbl <- mRNAs %>% mutate(proteinID=gsub(".mRNA$","",transcriptID))
  
  return( 
    list(
      genePosTbl=filter(geneTbl,geneID %in% IDtbl$geneID),
      IDtbl=IDtbl %>% select( geneID, proteinID, transcriptID)
    )
  )  
}

convertGFF_mRNA <- function(gff){
  gffTbl <-
    read_tsv(gff,col_names = c("seqname","feature","start","end","strand","attribute"),comment="#", col_types="c-cii-c-c")
  
  geneTbl <- 
    gffTbl %>% 
    filter(feature == "mRNA") %>% 
    gffSplitAttributes() %>%
    select(seqname,start,end,strand,ID) %>% 
    dplyr::rename(geneID=ID) 
  
  
  # protein IDs are same as and transcript IDs
  IDtbl <- geneTbl %>%
    mutate(proteinID = geneID, transcriptID = geneID)
  
  return( 
    list(
      genePosTbl = geneTbl,
      IDtbl=IDtbl %>% select( geneID, proteinID, transcriptID)
    )
  )  
}


convertGFF_accession <- function(gff){
  gffTbl <-
    read_tsv(gff,col_names = c("seqname","feature","start","end","strand","attribute"),comment="#", col_types="c-cii-c-c")
  
  mRNAs <-
    gffTbl %>% 
    filter(feature == "mRNA") %>% 
    select(attribute) %>%
    gffSplitAttributes() %>%
    transmute( geneID=Parent_Accession, transcriptID=Accession)
  
  CDSs <-
    gffTbl %>% 
    filter(feature == "CDS") %>% 
    select(attribute) %>%
    gffSplitAttributes() %>% 
    transmute( proteinID=Protein_Accession, transcriptID=Parent_Accession) %>% 
    distinct()
  
  geneTbl <- 
    gffTbl %>% 
    filter(feature == "gene") %>% 
    gffSplitAttributes() %>%
    select(seqname,start,end,strand,Accession) %>% 
    dplyr::rename(geneID=Accession)
  
  IDtbl <- left_join(CDSs,mRNAs,by="transcriptID")
  
  return( 
    list(
      genePosTbl=filter(geneTbl,geneID %in% IDtbl$geneID),
      IDtbl=IDtbl %>% select( geneID, proteinID, transcriptID)
    )
  )  
}
