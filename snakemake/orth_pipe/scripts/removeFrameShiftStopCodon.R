library(seqinr)
library(stringr)
args <- commandArgs(trailingOnly=TRUE)
nt <- read.fasta(args[1],as.string = T)
aa <- read.fasta(args[2],as.string = T)
posFrameShift <- str_locate_all(nt,pattern = '!')
posStopCodon <- str_locate_all(aa,pattern = '\\*')
for (i in which(sapply(posFrameShift,nrow) > 0 | sapply(posStopCodon,nrow) > 0)){
  posFS <- posFrameShift[[i]][,1]
  startFS <- 3*unique(floor((posFS-1)/3))+1
  endFS <- startFS+2
  posSC <- posStopCodon[[i]][,1]
  startSC <- 3*(posSC-1) + 1
  endSC <- startSC + 2
  start <- c(startFS, startSC)
  end <- c(endFS, endSC)
  s <- as.character(nt[[i]])
  for (j in seq_along(start)){
    substr(s, start[j], end[j]) <- "---"  
  }
  nt[[i]] <- s
}
write.fasta(sequences = nt,names = names(nt),file.out = args[3])
