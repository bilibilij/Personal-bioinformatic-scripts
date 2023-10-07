
setwd("F:/工作/3.诸葛菜/00.diploid/27.TE_ACR")
library("tidyverse")

TE_fam_bed<-read_tsv("input-data/core_table/2.TE_fam_ACR_tis.tsv")


TE_fam_bed %>%
  mutate(merge_ACR=11597683, genome_len=1198493200,AllTE_len=875517021 ) %>%
  select(X4,)


cc<-
cbind(TE_fam_bed[,1:18],TE_fam_bed[,37] ) %>% unique(.) %>%
  mutate(merge_ACR=11597683, genome_len=1198493200,AllTE_len=875517021 ) %>%
  group_by(X4) %>%
  mutate( TE_subfam_len=sum(X37), ES =(TE_subfam_len/merge_ACR)/(X14/AllTE_len )) %>%
  mutate(TE_subfam_count=n(), LER=log2(TE_subfam_count/(X15*TE_subfam_len/genome_len))   )  %>%
  #select(-X1,-X2,-X3) %>%
  arrange(desc(ES)) %>% filter( X15>=10, X14>=1198493 )

#1198493200/1000
ccc<- cc %>% select(X4, X5, TE_subfam_len, ES, TE_subfam_count, LER ) %>% unique(.) %>% arrange(desc(ES))



TE_fam_bed_nonmerge<-read_tsv("input-data/core_table/2.TE_fam_ACR_tis_nonMerge.tsv")



d<-
TE_fam_bed_nonmerge %>%
  group_by(X32) %>%
  select(X17,X18,X19)  %>%
  unique(.) %>%
  mutate(ACR_len=sum(X19-X18)) %>%
  select(X32,ACR_len)%>%
  unique(.)
  

dd<-
TE_fam_bed_nonmerge %>%
  #group_by(X32) %>%
  mutate( genome_len=1198493200,AllTE_len=875517021 ) %>%
  left_join(d,by=c("X32") ) %>%
  mutate(length=X19-X18 ) %>%
  select(-X9,-X10,-X11,-X16,-X17,-X18,-X19,-X24,-X25) %>%
  group_by(X4,X32,X26) %>%  #group X26 : private conserved shared
  mutate( TE_subfam_len=sum(X33), ES =(TE_subfam_len/ACR_len)/(X14/AllTE_len )) %>%
  mutate(TE_subfam_count=n(), LER=log2(TE_subfam_count/(X15*TE_subfam_len/genome_len))   )  %>%
  #select(-X1,-X2,-X3) %>%
  arrange(desc(ES)) %>% filter( X15>=10, X14>=1198493 ) %>%
  select(-X2,-X3,-X20,-X27,-X28,-X29,-X30,-X31, -genome_len, -AllTE_len, -flower, -leaf,-root,-seed,-stem, -tis_spe   )  %>%
  mutate(subfam_len=sum(length)  ) 



library(ggrepel)

dd %>%
  mutate(fraction=subfam_len/ACR_len  ) %>%
  ungroup() %>%
  select(X4, X5,X26,X32,ES,TE_subfam_count, TE_subfam_len,fraction) %>%
  unique(.) %>% 
  group_by(X26,X32) %>%
  arrange( desc(ES) ) %>%
  group_by(X32) %>%
  mutate(rank= seq(1,length.out = n())) %>%
  ggplot() +
  geom_point( aes(x=fraction,y=ES, size=TE_subfam_len, color=X26)) +
  geom_text_repel(aes(x=fraction,y=ES,label=ifelse( rank<=5,paste(X5,"\n",X4), "") ) ) +
  facet_wrap(~X32, nrow = 1)+
  theme_bw() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 











e<-
  TE_fam_bed_nonmerge %>%
  group_by(X32,X26) %>%
  select(X17,X18,X19,X26)  %>%
  unique(.) %>%
  mutate(ACR_len=sum(X19-X18)) %>%
  select(X32,ACR_len,X26)%>%
  unique(.)


ee<-
  TE_fam_bed_nonmerge %>%
  #group_by(X32) %>%
  mutate( genome_len=1198493200,AllTE_len=875517021 ) %>%
  left_join(e,by=c("X32","X26") ) %>% 
  mutate(length=X19-X18 ) %>%
  select(-X9,-X10,-X11,-X16,-X17,-X18,-X19,-X24,-X25) %>%
  group_by(X4,X32,X26) %>%  #group X26 : private conserved shared
  mutate( TE_subfam_len=sum(X33), ES =(TE_subfam_len/ACR_len)/(X14/AllTE_len )) %>%
  mutate(TE_subfam_count=n(), LER=log2(TE_subfam_count/(X15*TE_subfam_len/genome_len))   )  %>%
  #select(-X1,-X2,-X3) %>%
  arrange(desc(ES)) %>%
  filter( X15>=10, X14>=1198493 ) %>%
  select(-X2,-X3,-X20,-X27,-X28,-X29,-X30,-X31, -genome_len, -AllTE_len, -flower, -leaf,-root,-seed,-stem, -tis_spe   )  %>%
  mutate(subfam_len=sum(length)  ) 
  
ee %>%
  mutate(fraction=subfam_len/ACR_len  ) %>%
  ungroup() %>%
  select(X4, X5,X26,X32,ES,TE_subfam_count, TE_subfam_len,fraction) %>%
  unique(.) %>% 
  group_by(X26,X32) %>%
  arrange( desc(ES) ) %>%
  group_by(X32) %>%
  mutate(rank= seq(1,length.out = n())) %>%
  ggplot() +
  geom_point( aes(x=fraction,y=ES, size=TE_subfam_len, color=X26)) +
  geom_text_repel(aes(x=fraction,y=ES,label=ifelse( rank<=5,paste(X5,"\n",X4), "") ) ) +
  facet_wrap(~X32, nrow = 1)+
  theme_bw() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 


TE_fam_CN<-read_tsv("input-data/2.TE_ACR_subfamily/TE_copy_number.tsv")

head(ccc,n = 10) %>%
  left_join(TE_fam_CN, by=c("X4"="TE_famliy") ) %>%
  mutate(TE=paste(X4,":",X5)) %>%
  select(-X4,-X5, -TE_subfam_len, -TE_subfam_count, -TE )  %>%
  gather( key = sp, value= copy_number, 4:22 ) %>%
  ggplot( aes(x=sp,y=copy_number,color=X4, group=1 ) ) +
  geom_line()


head(ccc,n = 10) %>%
  left_join(TE_fam_CN, by=c("X4"="TE_famliy") ) %>%
  mutate(TE=paste(X4,":",X5)) %>%
  select(-X4,-X5, -TE_subfam_len, -TE_subfam_count, -TE )  %>%
  gather( key = sp, value= copy_number, 4:22 ) %>%
  ggplot( aes(x=sp,y=copy_number,fill=X4, group=1 ) ) +
  geom_col()+
  theme_bw() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 

head(ccc,n = 10) %>%
  left_join(TE_fam_CN, by=c("X4"="TE_famliy") )

TE_fam_CN[which(TE_fam_CN$TE_famliy %in%lev),]



#private conserved shared 


a<-ee %>%
 #- dd%>%
  mutate(fraction=subfam_len/ACR_len  ) %>%
  ungroup() %>%
  select(X4, X5,X26,X32,ES,TE_subfam_count, TE_subfam_len,fraction) %>%
  unique(.) %>% 
  group_by(X26,X32) %>%
  arrange( desc(ES) ) %>%
  group_by(X32) %>%
  mutate(rank= seq(1,length.out = n())) 

ap<-a %>% filter(X26=="private")
ac<-a %>% filter(X26=="conserved")
as<-a %>% filter(X26=="shared")

aa<-inner_join (inner_join(ap, ac,by=c("X4", "X32") ),as ,by=c("X4", "X32")   )


#

ap<-a %>% filter(X26=="private")
ac<-a %>% filter(X26=="conserved")
as<-a %>% filter(X26=="shared")

aa<-full_join (full_join(ap, ac,by=c("X4", "X32") ),as ,by=c("X4", "X32")  )


p<-
a %>%
  group_by(X32,X26) %>%
  arrange(desc(ES) ) %>%
  mutate(rank=seq(1,length.out = n())) %>%
  filter(rank<=5) %>%
  group_by(X32) %>%arrange(desc(ES) ) %>%  mutate(rank2=seq(1,length.out = n()))  %>% ungroup() %>%
  select(X4,X32,rank,rank2) %>% unique(.) %>% 
  arrange(desc(rank2))  %>%
  filter(rank2<=10) %>%
  select(-rank,-rank2) 

lev<-(p %>% select(X4) %>% unique(.)  %>% left_join( ccc,by=c("X4") ) %>% arrange(desc(ES)) %>% select(X4,ES))$X4

rbind(rbind(left_join(p,ap, by=c("X4","X32") )%>% mutate(class="private") ,
            left_join(p,as, by=c("X4","X32") ) %>% mutate(class="shared")) , 
      left_join(p,ac, by=c("X4","X32") )  %>% mutate(class="conserved"))  %>% #arrange(desc(fraction))
  #group_by(X4) %>% add_count() %>% arrange(desc(n))
  #mutate_all(.funs = function(x){ ifelse(is.na(x), 0 , x) }) %>%  #arrange(desc(ES))
  ggplot() +
  geom_point(aes(y=factor(X4, levels = rev(lev)) ,x=X32,fill=log2(ES),size=fraction, group=1), shape=21) +
  facet_wrap(~factor(class,levels = c("shared","conserved", "private") ))  +
  scale_fill_gradient( low="white", high="#2E62AD")+
  theme_bw() +
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 













  left_join(ap, by=c("X4","X32") ) %>%
  mutate_all(.funs = function(x){ ifelse(is.na(x), 0 , x) }) %>%
  ggplot(aes(y=X4,x=X32,color=ES.x,size=fraction.x, group=1)) +
  geom_point() +
  scale_color_continuous()

a %>%
  group_by(X32,X26) %>%
  arrange(desc(ES) ) %>%
  mutate(rank=seq(1,length.out = n())) %>%
  filter(rank<=5) %>%
  group_by(X32) %>%arrange(desc(ES) ) %>%  mutate(rank2=seq(1,length.out = n()))  %>% ungroup() %>%
  select(X4,X32,rank,rank2) %>% unique(.) %>% 
  arrange(desc(rank2))  %>%
  filter(rank2<=10) %>%
  select(-rank,-rank2) %>%
  left_join(aa, by=c("X4","X32") ) %>%
  mutate_all(.funs = function(x){ ifelse(is.na(x), 0 , x) }) %>%
  ggplot(aes(y=X4,x=X32,color=ES.y,size=fraction.y, group=1)) +
  geom_point()


a %>%
  group_by(X32,X26) %>%
  arrange(desc(ES) ) %>%
  mutate(rank=seq(1,length.out = n())) %>%
  filter(rank<=5) %>%
  group_by(X32) %>%arrange(desc(ES) ) %>%  mutate(rank2=seq(1,length.out = n()))  %>% ungroup() %>%
  select(X4,X32,rank,rank2) %>% unique(.) %>% 
  arrange(desc(rank2))  %>%
  filter(rank2<=10) %>%
  select(-rank,-rank2) %>%
  left_join(aa, by=c("X4","X32") ) %>%
  mutate_all(.funs = function(x){ ifelse(is.na(x), 0 , x) }) %>%
  ggplot(aes(y=X4,x=X32,color=ES,size=fraction, group=1)) +
  geom_point()
  #select(-ES,-TE_subfam_count,-fraction) %>%
