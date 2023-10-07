setwd("F:/工作/3.诸葛菜/00.diploid/27.TE_ACR")
library("tidyverse")


dup<-read_tsv("input-data/0.TE_ACR_filter/dup_gene.tsv")

remove<-(dup %>% filter(class!="wgd") %>% mutate(wgd=paste(dup1,dup2, sep="_")))$wgd


TE_fam_bed<-read_tsv("input-data/2.TE_ACR_subfamily/OV_ACR_merge_5tis_TEfam.bed", col_names = F) %>%
  mutate_at( .vars = vars(X35), function(x){ gsub( "Seed", "seed", gsub( "^S$" ,"stem",
                    gsub( "R", "root",  gsub( "L", "leaf",  gsub( "F", "flower", x ) ) )) )      } )


TE_fam_bed<-TE_fam_bed[which(!(TE_fam_bed$X23 %in% remove)),]


TE_tis<-read_tsv("input-data/core_table/1.TE-ACR-tis.tsv")

TE_tis<-TE_tis[which(!(TE_tis$X4 %in% remove)),]



TE_fam_bed<-left_join(TE_fam_bed, TE_tis %>% select(X6, X21, tis_num, flower, leaf, root, seed, stem, tis_spe,X4 ), by=c("X35"="X21", "X25"="X6", "X23"="X4" )    )

write_tsv(TE_fam_bed,  "input-data/core_table/2.TE_fam_ACR_tis.tsv", col_names = T  )




#non-merge
TE_fam_bed_nonmerge<-read_tsv("input-data/2.TE_ACR_subfamily/OV_ACR_5tis_TEfam.bed", col_names = F) %>%
  mutate_at( .vars = vars(X32), function(x){ gsub( "Seed", "seed", gsub( "^S$" ,"stem",
                                                                         gsub( "R", "root",  gsub( "L", "leaf",  gsub( "F", "flower", x ) ) )) )      } )

TE_fam_bed_nonmerge<-TE_fam_bed_nonmerge[which(!(TE_fam_bed_nonmerge$X20 %in% remove)),]

TE_fam_bed_nonmerge<-left_join(TE_fam_bed_nonmerge, TE_tis %>% select(X6, X21, tis_num, flower, leaf, root, seed, stem, tis_spe,X4 ), 
                               by=c("X32"="X21", "X22"="X6", "X20"="X4" )    )

write_tsv(TE_fam_bed_nonmerge,  "input-data/core_table/2.TE_fam_ACR_tis_nonMerge.tsv", col_names = T  )










#old analysis test



TE_bed<-data.frame()
for (i in dir("input-data/0.TE_ACR_filter") %>% grep("spTE.bed$", . ,value = T) ){
  
  TE_bed  <-rbind(TE_bed,
                  read_tsv( paste("input-data/0.TE_ACR_filter", i , sep="/"), col_names = F ) %>%
                    mutate(class= gsub( ".geneDistance.closest.ACR.ACR.conserved.bed.sort.ts.all.spTE.bed"  ,"",
                                        gsub("OV_ATAC_", "",i) )  ) %>% 
                    mutate_at( .vars = vars(class), function(x){ gsub( "Seed", "seed", gsub( "^S$" ,"stem",  gsub( "R", "root",  gsub( "L", "leaf",  gsub( "F", "flower", x ) ) )) )      } )
  )
}


dup<-read_tsv("input-data/0.TE_ACR_filter/dup_gene.tsv")

remove<-(dup %>% filter(class!="wgd") %>% mutate(wgd=paste(dup1,dup2, sep="_")))$wgd

TE_bed<-TE_bed[which(!(TE_bed$X4 %in% remove)),]



TE_fam<-read_tsv("input-data/2.TE_ACR_subfamily/TE_subfamily_stat.tsv", col_names = F)

c<-as.numeric((TE_bed %>%
  filter(class=="flower") %>%
  select(X1,X2,X3)  %>%
  unique(.) %>%
  mutate(sum=sum(X3-X2)))[1,4])



left_join(TE_bed,TE_fam,by=c("X16"="X1")) %>% filter(class=="flower") %>%
  mutate(genome_len=1198493200, ACR_region=c) %>% select(X17, X2.y,X3.y,genome_len,ACR_region,X20)%>%
  filter(X17!=".") %>%
  group_by(X17) %>%
  mutate(TEfam_num=n(), OpenTE_len=sum(X20) ) %>% select(-X20) %>% unique(.) %>%
  mutate(LER=log2(TEfam_num/X3.y/OpenTE_len*genome_len  )       ) %>%
  arrange(desc(LER) ) %>% filter(X3.y > 10, OpenTE_len>50)
  



cc<-
left_join(TE_bed,TE_fam,by=c("X16"="X1")) %>% filter(class=="flower") %>%
  mutate(genome_len=1198493200, ACR_region=c, AllTE_len=875517021) %>% select(X16, X2.y,X3.y,genome_len,ACR_region,X20,AllTE_len,X15)%>%
  filter(X16!=".") %>%
  group_by(X16) %>%
  mutate(TEfam_num=n(), OpenTE_len=sum(X20) )  %>% select(-X20) %>% unique(.) %>%
  mutate(ER=OpenTE_len/ACR_region/X2.y*AllTE_len ) %>% filter(X3.y>10) %>%
  arrange( desc(TEfam_num), desc(ER)  )



