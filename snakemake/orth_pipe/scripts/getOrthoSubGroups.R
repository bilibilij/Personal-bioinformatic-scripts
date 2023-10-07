library(tidyverse)
library(ape)
library(treeio)


getRootNode <- function(phylo){
    # get the node that has no edges pointing to it
    setdiff(phylo$edge[ ,1], phylo$edge[ ,2])
}

# get all ancestors including given node back to root
getAncestorNodes <- function(phylo,node){
    # recursive call to get ancestral nodes
    getAncestors <- function(curNode,ancestorNodes){
        ancNode = phylo$edge[phylo$edge[,2]==curNode,1]
        if( length(ancNode) == 1)
            ancestorNodes <- Recall(curNode=ancNode,ancestorNodes)
        return( c(curNode,ancestorNodes))
    }
    
    getAncestors(curNode=node,ancestorNodes=integer(0))  
}


#' Get orthogroups by traversing the tree and spliting it at ancestral duplication nodes
#'
#' @param treedata treedata class gene tree
#' @param ancNodeLabel vector of names of ancestral nodes
#'
#' @return Table of genes with sub-orthogroup ids given by ancestral post-duplication node numbers
splitOrthoGroups <- function(spcTree, treedata, OGid){
    nbTips <- length(spcTree$tip.label)
    ancNodeLabel <- spcTree$node.label
    # prepare data for the recursion call
    
    fromNode <- treedata@phylo$edge[ ,1]
    toNode <- treedata@phylo$edge[ ,2]
 

    treedata@data <- treedata@data[order(treedata@data$node),]
    # sanity check
    if( !identical(treedata@data$node,1:(Ntip(treedata@phylo)+Nnode(treedata@phylo))) )
        stop("treedata table not ordered")
    
    # if there are DD entries then make DD not dup tip
    isDD <- rep(NA,nrow(treedata@data))
    if("DD" %in% names(treedata@data)){
        isDD <- treedata@data$DD
    }
    isD <- rep(NA,nrow(treedata@data))
    if("D" %in% names(treedata@data)){
      isD <- treedata@data$D
    }
    OGid <- c(OGid,rep(NA_character_,length(ancNodeLabel)+1))
    names(OGid) <- c("OG",ancNodeLabel,"Tip")
    retVal <- tibble(OG=NA_character_,geneID=treedata@phylo$tip.label) 
    for (anc in ancNodeLabel){
        retVal[[anc]] <- NA_character_
    }
    retVal[["Tip"]] = NA_character_
    isFix <- OGid
    # define recursion function
    recurseFUN <- function(curNode, OGid, isFix){
        childNodes <- toNode[fromNode==curNode]
        curSpecies <- treedata@data$S[curNode]
        newOGid = OGid
        if (curSpecies %in% ancNodeLabel){
            if (is.na(isFix[curSpecies]) && isD[curNode]!="Y") newOGid[curSpecies] <- curNode
        } else {
            if (is.na(newOGid["Tip"])) newOGid["Tip"] <- curNode
        }
        ##consider the case DD=Y
        newIsFix = isFix
        if (!is.na(isDD[curNode])){#if the current node is DD
            newIsFix[curSpecies] <- TRUE
        }
        ##consider the lost event case
        if (curSpecies %in% ancNodeLabel){
            ancOfcurSpecies <- spcTree$node.label[getAncestorNodes(spcTree,which(spcTree$node.label==curSpecies)+nbTips)-nbTips][-1]    
        } else{
            ancOfcurSpecies <- spcTree$node.label[getAncestorNodes(spcTree,which(spcTree$tip.label==curSpecies))[-1]-nbTips]
        }
        if (length(ancOfcurSpecies)>0 && is.na(newOGid[ancOfcurSpecies[1]])){
            notNA <- which(is.na(newOGid[ancOfcurSpecies])==FALSE)
            if (length(notNA)>0){
                fillUptoNode <- notNA[1]-1
            } else{
                fillUptoNode <- length(ancOfcurSpecies)
            }
            lostAnc <- ancOfcurSpecies[1:fillUptoNode]
            newOGid[lostAnc] <- curNode
        }
        ##call the children recursively
        if (length(childNodes) == 0){
            retVal[curNode,names(newOGid)] <<- as.list(newOGid)
        } else {
            for (childNode in childNodes){
                Recall(curNode = childNode, OGid = newOGid, isFix = newIsFix)    
            }
        }
    }
    
    # call recursive function starting from root
    # return( recurseFUN( curNode = getRootNode(treedata@phylo), OGid=OGid) )
    recurseFUN( curNode = getRootNode(treedata@phylo), OGid=OGid, isFix = isFix) 
    cols <- setdiff(colnames(retVal),c("OG","geneID"))
    retVal[,cols] <- apply(retVal[,cols],2,function(col){
        gsub(".+.NA$",NA,paste0(retVal$OG,".",col))
    }
    )
    return(retVal)
}

getOrthoGroupsByMRCA <- function(trees, MRCAnodeLabel, spcTree){
  MRCAnode <- match(MRCAnodeLabel,spcTree$node.label)+Ntip(spcTree)
  spcs <- extract.clade(spcTree,node = MRCAnode)$tip.label
  
  # get ancestral nodes
  ancNodes <- getAncestorNodes(spcTree,MRCAnode)
  ancNodeLabel <- spcTree$node.label[ancNodes - Ntip(spcTree)]
  
  # iterarate through the trees
  imap_dfr( trees, ~ splitOrthoGroups(treedata=.x,ancNodeLabel,OGid = .y)) %>% 
    mutate( spc = sub(".*_([^_]+)$","\\1",geneID),
            geneID = sub("_[^_]+$","",geneID)) %>% 
    filter( spc %in% spcs)
}

getOrthoGroupsBySpcs <- function(trees, spcs, spcTree){
  MRCAnode <- getMRCA(spcTree,tip = spcs)
  
  # get ancestral nodes
  ancNodes <- getAncestorNodes(spcTree,MRCAnode)
  ancNodeLabel <- spcTree$node.label[ancNodes - Ntip(spcTree)]
  
  # iterarate through the trees
  imap_dfr( trees, ~ splitOrthoGroups(treedata=.x,ancNodeLabel,OGid = .y)) %>% 
    mutate( spc = sub(".*_([^_]+)$","\\1",geneID),
            geneID = sub("_[^_]+$","",geneID)) %>% 
    filter( spc %in% spcs)
}


