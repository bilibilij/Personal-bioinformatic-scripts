setwd("F:/工作/3.诸葛菜/00.diploid/27.TE_ACR")
library("tidyverse")

TE_ACR_final<-read_tsv("input-data/core_table/0.TE-ACR.tsv",col_names = F  )


# pie plot of TE,non-TE, dupACR

TE_ACR_final %>%
  select(X10, X21, X19) %>%
  mutate(class=ifelse( X19==".", "Non-TE", "TE" ) ) %>%
  select(-X19) %>%
  #stat count of Non-TE, ACR class
  group_by(X10,X21,class) %>%
  add_count() %>%
  ungroup()  %>% 
  group_by(X21) %>%
  add_count() %>% 
  group_by(X21,class) %>% 
  add_count()%>%
  group_by(X21) %>% 
  #add ymin and ymax of each ACRclass and w/wo TE
  arrange(X21,factor(class,levels = c("TE","Non-TE")), factor(X10,levels=c("private", "conserved", "shared") )) %>%
  mutate(prop=seq(0,1,length.out=n())) %>% 
  group_by(X10,X21,class)%>% 
  mutate(ymin=min(prop),ymax=max(prop)) %>% 
  group_by(X21,class) %>%
  mutate(yTEmin=min(prop),yTEmax=max(prop)) %>% 
  select(-prop) %>% 
  unique(.) %>%
  #add label of each class 
  mutate(pie1=n/nn,pie2=nnn/nn,  label1=paste(round((n/nn)*100,2), "%", sep=""  ),label2=paste(round((nnn/nn)*100,2), "%", sep="" )) %>%
  ggplot() +
  geom_rect(aes(fill=X10,ymax=ymax,ymin=ymin, color=X10, xmax=4, xmin=3,))+
  geom_rect(aes(fill=class,ymax=ymax,ymin=ymin, color=class, xmax=0, xmin=2.95 ))+
  geom_text(aes(4.7, as.numeric(  (ymax+ymin)/2), label=label1), color="black")+
  geom_text(aes(2, as.numeric(yTEmax+yTEmin )/2, label=label2 ), color="white" )+
  facet_wrap(~X21, nrow = 1)+ xlim(c(0,5)) +
  #theme(aspect.ratio=1) + 
  coord_polar(theta="y") +
  theme_void() +
  theme(legend.position = 'top') + # 隐
  scale_fill_manual(values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD", "TE"= "#bb2f28", "Non-TE"= "#494544" ))+
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD", "TE"= "#bb2f28", "Non-TE"= "#494544"   ))

ggsave("1.multi_pie.pdf")



#Ov-shared, Ov-specific TE count plot 

TE_ACR_final %>%
  select( X21, X19) %>%
  filter(X19!=".") %>%
  arrange(X21,X19) %>%
  group_by(X21) %>%
  add_count() %>%
  group_by(X21,X19) %>%
  add_count() %>%
  mutate( prop=nn/n ) %>%
  unique(.) %>%
  mutate(label= paste(round((prop)*100,2), "%", sep = "" ) ) %>%
  ggplot(aes(x=X21,y=prop,fill=X19)) +
  geom_col() +
  geom_text(aes(x=X21, y=prop, label=label), color="white"  ) +
  theme_classic() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 



TE_ACR_final %>%
  select(X10,X19,X21)%>%
  filter(X19!=".") %>%
  group_by(X10,X19,X21) %>%
  add_count() %>%
  group_by(X10,X21) %>%
  add_count() %>%
  mutate( prop=n/nn ) %>%
  unique(.) %>%
  mutate(label= paste(round((prop)*100,2), "%", sep = "" ) ) %>%
  ggplot(aes(x=factor(X10, levels=c("shared","conserved","private")), y=prop, fill= X19,) )+
  geom_col() +
  geom_text(aes(x=factor(X10, levels=c("shared","conserved","private")), y=prop, label=label), color="white"  )+
  facet_wrap(~X21, nrow = 1) +
  theme_classic() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 




#TE superfamily 

TE_ACR_final %>%
  filter(X19!=".") %>%
  select(X10,X15,X21) %>%
  mutate_at(.vars=vars(X15), .funs = function(x){ gsub("/.*","", x)  }) %>%
  group_by(X10,X15,X21) %>%
  add_count() %>%
  group_by(X10,X21) %>%
  add_count() %>%
  mutate( prop=n/nn ) %>%
  unique(.) %>%
  mutate(label= paste(round((prop)*100,2), "%", sep = "" ) ) %>% 
  filter(X15!="pararetrovirus", X15!="Simple_repeat") %>%
  arrange(desc(prop)) %>%
  ggplot(aes(x=factor(X10, levels=c("shared","conserved","private")), y=prop, fill= factor(X15,levels=rev(c("LINE", "Unknown","MITE","DNA", "LTR")) )  ,) )+
  geom_col() +
  #geom_text(aes(x=factor(X10, levels=c("shared","conserved","private")), y=prop, label=label), color="white"  )+
  facet_wrap(~X21, nrow = 1) +
  theme_classic() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 
  


TE_ACR_final %>%
  filter(X19!=".") %>%
  select(X10,X15,X21) %>%
  mutate_at(.vars=vars(X15), .funs = function(x){ gsub("/.*","", x)  }) %>%
  group_by(X10,X15,X21) %>%
  add_count() %>%
  group_by(X10,X21) %>%
  add_count() %>%
  mutate( prop=n/nn ) %>%
  unique(.) %>%
  mutate(label= paste(round((prop)*100,2), "%", sep = "" ) ) %>% 
  filter(X15!="pararetrovirus", X15!="Simple_repeat") %>%
  arrange(desc(prop)) %>%
  ggplot(aes(x=factor(X10, levels=c("shared","conserved","private")), y=prop, fill= factor(X15,levels=c("LINE", "Unknown","MITE","DNA", "LTR") )  ,) )+
  geom_col() +
  #geom_text(aes(x=factor(X10, levels=c("shared","conserved","private")), y=prop, label=label), color="white"  )+
  facet_wrap(~X21, nrow = 1) +
  theme_classic() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 

### LTR butaihao 

TE_ACR_final %>%
  filter(X19!=".") %>%
  select(X10,X15,X21) %>%
  #filter(grepl("LTR",X15)) %>%
  #mutate_at(.vars=vars(X15), .funs = function(x){ gsub("/.*","", x)  }) %>%
  mutate( class= gsub("/.*","", X15)  ) %>%
  #rename(  X155=class,X15=class  ) %>%
  group_by(X10,X15,X21) %>%
  add_count() %>%
  group_by(X10,X21) %>%
  add_count() %>%
  mutate( prop=n/nn ) %>%
  unique(.) %>%
  mutate(label= paste(round((prop)*100,2), "%", sep = "" ) ) %>% 
  filter(X15!="pararetrovirus", X15!="Simple_repeat") %>%
  filter(grepl("LTR",X15)) %>%
  ggplot(aes(x=factor(X10, levels=c("shared","conserved","private")), y=prop, fill=X15 ) )+
  geom_col(position="dodge") +
  #geom_text(aes(x=factor(X10, levels=c("shared","conserved","private")), y=prop, label=label), color="white"  )+
  facet_wrap(~X21, nrow = 1) +
  theme_classic() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 




#sub family

TE_ACR_final %>%
  filter(X19!=".") %>%
  select(X10,X16,X21) 
  
  
  
