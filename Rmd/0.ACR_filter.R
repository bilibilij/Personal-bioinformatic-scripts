setwd("F:/工作/3.诸葛菜/00.diploid/27.TE_ACR")
library("tidyverse")



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



TE_bed_new<-TE_bed %>% group_by(class,X6,X4) %>%
  filter( !grepl("repeat_region",X17,  ) ) %>%
  #mutate(a= seq(0,1, length.out=n()) )  #this add a acol with 0, 0.33, 0.66, 1
  #mutate(TE_num= seq(1,length.out=n())  )  %>% ungroup() %>%    # this add a col with 1,2,3
  #mutate(TE_num= n(), TE_len=X14-X13,   ) %>%  # this add a col with 3,3,3, which is the suitable choice for filtering multi-TE-ACR (preference)
  #filter(TE_num>1)  %>%
  add_count() #add_count could be equal to n()


TE_bed_new_multi<-
  TE_bed_new %>% filter(n>1) %>% select(-n) %>%
  mutate(sp=ifelse(X19=="Ov_specific", 1, 0 ) , stat=sum(sp)) %>%  #first filter with Ov-specific preference
  filter(!(sp==0&stat>=1 )) %>% select(-sp, -stat) %>%
  mutate( unknown=ifelse(X15=="Unknown", 0, 1 ) , stat=sum(unknown) ) %>%
  filter(!(unknown==0&stat>=1)) %>% select(-unknown, -stat) %>%  #if unknown overlap with the other TE, we prefer to use known TE
  arrange(desc(X20)) %>%
  mutate( dis_rank= seq(1,length.out=n() ), max=max(X20) ) %>%
  filter(!(dis_rank>1&X20<max )) %>% select(-dis_rank, -max) %>% #filter length between TE and ACR, coverage much better is OK
  group_by(class,X6,X4) %>%
  arrange(X14-X13 ) %>%
  mutate( len_rank= seq(1,length.out=n() ), max=max(len_rank), TE_len=X14-X13 ) %>% # TE length  with smaller one
  filter( len_rank==1 ) %>%
  select(-len_rank,-max,-TE_len)
          

TE_ACR_final<-rbind(TE_bed_new%>% filter(n==1) %>% select(-n), TE_bed_new_multi  )  


dir.create("input-data/core_table/")
write_tsv(TE_ACR_final,"input-data/core_table/0.TE-ACR.tsv", col_names = F )


#(c %>% filter(X19=="Ov_specific") )$X6 
#dim(c)
#c[ %in% c$X6),'a']=1

#which((c %>% filter(X19=="Ov_specific"))$X6) 

#(c %>% filter(X19=="Ov_specific"))$X6
  
  
#seq(0,round(c,2), length.out=5)
