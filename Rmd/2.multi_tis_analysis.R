setwd("F:/工作/3.诸葛菜/00.diploid/27.TE_ACR")
library("tidyverse")


TE_ACR_filter_tis<-read_tsv("input-data/core_table/1.TE-ACR-tis.tsv",col_names = T  )

# the overview of 1,2,3,4,5tissues overlap ACR proportion
P1<-
TE_ACR_filter_tis %>%
  select(X6,X10,X11,X21,tis_num) %>%
  mutate_at(.vars = vars(tis_num), .funs = function(x){ paste(x,"tissues", sep = " ") } ) %>%
  group_by(X21,tis_num) %>%
  add_count() %>%
  group_by(X21)%>%
  add_count() %>%
  mutate(prop=n/nn) %>%
  select(-X6,-X10,-X11) %>%
  unique(.)%>%
  ggplot(aes(x=X21,y=prop,fill=tis_num)) +
  geom_col() +
  theme_classic() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) +
  xlab("") +
  ylab("Proportion of overall duplicates ACR")

P1
ggsave("results/2.tis_spe/0.overallTisSpe.pdf")
ggsave("results/2.tis_spe/0.overallTisSpe.PNG")

#
P2<-
TE_ACR_filter_tis %>%
  select(X6,X10,X11,X21,tis_num) %>%
  mutate_at(.vars = vars(tis_num), .funs = function(x){ paste(x,"tissues", sep = " ") } ) %>%
  group_by(X21,tis_num) %>%
  add_count() %>%
  group_by(X21)%>%
  add_count() %>%
  mutate(prop=n/nn) %>% #above line has calculated 1,2,3,4,5 tissues shared ACR proprotion in each tissue
  group_by(X21,tis_num,X10) %>%
  add_count() %>%
  mutate(ACR_class_prop=nnn/n ) %>%
  select(-X6,-X11) %>%
  unique(.) %>%
  ggplot(aes(x=tis_num,y=ACR_class_prop, fill=factor(X10, levels=c("shared", "private", "conserved")) )) +
  geom_col(position = "stack") +
  facet_wrap(~X21, nrow = 1) +
  theme_classic() +
  scale_fill_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" ))+
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) +
  xlab("")+ylab("Proportion of ACR")

P2
ggsave("results/2.tis_spe/1.dupACR_tis_spe.pdf", width = 8, height = 5)
ggsave("results/2.tis_spe/1.dupACR_tis_spe.PNG", width = 8, height = 5)


# relationship between tissue specificity and TE within ACR
P3<-
TE_ACR_filter_tis %>%
  mutate(class=ifelse( X19==".", "Non-TE", "TE" ) ) %>%
  select(X6,X10,X11,X21,tis_num,class) %>%
  mutate_at(.vars = vars(tis_num), .funs = function(x){ paste(x,"tissues", sep = " ") } ) %>%
  group_by(X21,tis_num) %>%
  add_count()%>%
  group_by(X21,tis_num,X10,class)%>%
  add_count() %>%
  mutate(prop=nn/n,a=sum(prop)) %>%
  select(-X6) %>% unique(.) %>%
  ggplot()+
  geom_segment(aes(y=0, yend= prop, x=tis_num, xend=tis_num ), size=2, alpha=0.2, color="grey") +
  geom_point(aes(x=tis_num, y=prop, color=X10) , size = 5) +
  #geom_line(aes(x=tis_num , y=prop, group=1, color=X10 ))+
  facet_grid(cols = vars(X21), rows = vars(class), scales = "free" ) +
  theme_classic() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" )) +
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) +
  xlab("") + ylab("Proportion of TE-ACR")
P3
ggsave("results/2.tis_spe/2.1.TE_dupACR_tis_spe.pdf", width = 8, height = 5)
ggsave("results/2.tis_spe/2.1.TE_dupACR_tis_spe.PNG", width = 8, height = 5)

  
P4<-
TE_ACR_filter_tis %>%
  mutate(class=ifelse( X19==".", "Non-TE", "TE" ) ) %>%
  select(X6,X10,X11,X21,tis_num,class) %>%
  mutate_at(.vars = vars(tis_num), .funs = function(x){ paste(x,"tissues", sep = " ") } ) %>%
  group_by(X21,tis_num) %>%
  add_count()%>%
  group_by(X21,tis_num,X10,class)%>%
  add_count() %>%
  mutate(prop=nn/n,a=sum(prop)) %>%
  select(-X6) %>% unique(.) %>%
  ggplot()+
  geom_segment(aes(y=0, yend= prop, x=tis_num, xend=tis_num ), size=2, alpha=0.2, color="grey") +
  geom_point(aes(x=tis_num, y=prop, color=class) , size = 5) +
  #geom_line(aes(x=tis_num,y=prop))+
  facet_grid(cols = vars(X21), rows = vars(X10), scales = "free" ) +
  theme_classic() +
  #scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" )) +
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) +
  xlab("") + ylab("Proportion of TE-ACR")
P4
ggsave("results/2.tis_spe/2.2.TE_dupACR_tis_spe.pdf", width = 9, height = 7)
ggsave("results/2.tis_spe/2.2.TE_dupACR_tis_spe.PNG", width = 9, height = 7)




P5<-
TE_ACR_filter_tis %>%
  mutate(class=ifelse( X19==".", "Non-TE", "TE" ) ) %>%
  select(X6,X10,X11,X21,tis_num,class) %>%
  mutate_at(.vars = vars(tis_num), .funs = function(x){ paste(x,"tissues", sep = " ") } ) %>%
  group_by(X21,tis_num,X10) %>%
  add_count()%>%
  group_by(X21,tis_num,X10,class)%>%
  add_count() %>%
  mutate(prop=nn/n,a=sum(prop)) %>%
  select(-X6) %>% unique(.) %>%
  #filter(class=="TE") %>%
  ggplot()+
  geom_segment(aes(y=0, yend= prop, x=tis_num, xend=tis_num ), size=2, alpha=0.2, color="grey") +
  geom_point(aes(x=tis_num, y=prop, color=X10, shape=class, alpha=class) , size = 5) +
  #geom_line(aes(x=tis_num , y=prop, group=1, color=X10 ))+
  #facet_grid(cols = vars(X21), rows = vars(class), scales = "free" ) +
  facet_wrap(~X21,nrow = 1)+
  theme_classic() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" )) +
  scale_alpha_manual(values=c("Non-TE"=0.5, "TE"=1), ) +
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5))
P5
ggsave("results/2.tis_spe/2.3.TE_dupACR_tis_spe.pdf", width = 8, height = 6)
ggsave("results/2.tis_spe/2.3.TE_dupACR_tis_spe.PNG", width = 8, height = 6)


P6<-
TE_ACR_filter_tis %>%
  mutate(class=ifelse( X19==".", "Non-TE", "TE" ) ) %>%
  select(X6,X10,X11,X21,tis_num,class) %>%
  mutate_at(.vars = vars(tis_num), .funs = function(x){ paste(x,"tissues", sep = " ") } ) %>%
  group_by(X21,tis_num,X10) %>%
  add_count()%>%
  group_by(X21,tis_num,X10,class)%>%
  add_count() %>%
  mutate(prop=nn/n) %>%
  select(-X6,-X11,-n,-nn,) %>% unique(.) %>%
  spread(key = class, value=prop) %>%
  mutate(enr=TE/`Non-TE`) %>%
  #filter(class=="TE") %>%
  ggplot()+
  geom_segment(aes(y=0, yend= enr, x=tis_num, xend=tis_num ), size=2, alpha=0.2, color="grey") +
  geom_point(aes(x=tis_num, y=enr, color=X10) , size = 5) +
  #geom_line(aes(x=tis_num , y=prop, group=1, color=X10 ))+
  #facet_grid(cols = vars(X21), rows = vars(class), scales = "free" ) +
  facet_wrap(~X21,nrow = 1)+
  theme_classic() +
  scale_color_manual( values = c("shared"="#EEE12D", "conserved" = "#AB5785", "private" = "#2E62AD" )) +
  theme(
    axis.title.x=element_text(size=13,color="black",hjust=0.5, face="italic", vjust = -1),                                                                           axis.title.y=element_text(size=13,color="black", vjust=2.5),
    axis.text.x=element_text(size=13,color="black", angle = 90),
    axis.text.y=element_text(size=13,color="black"),
    legend.text=element_text(size=13,color="black"),
    legend.title=element_blank(),
    plot.title=element_text(size=13,color="black",face="italic",hjust=0.5)) 
P6
ggsave("results/2.tis_spe/2.4.TE_dupACR_tis_spe.pdf", width = 8, height = 6)
ggsave("results/2.tis_spe/2.4.TE_dupACR_tis_spe.PNG", width = 8, height = 6)




pheat_plot<-list()
library(pheatmap)


for(i in c("shared", "conserved", "private")){
  c<-
    TE_ACR_filter_tis %>%
    mutate(class=ifelse( X19==".", "Non-TE", "TE" ) ) %>%
    select(X6,X10,X11,X21,tis_num,class) %>%
    mutate_at(.vars = vars(tis_num), .funs = function(x){ paste(x,"tissues", sep = " ") } ) %>%
    group_by(X21,tis_num,X10) %>%
    add_count()%>%
    group_by(X21,tis_num,X10,class)%>%
    add_count() %>%
    ungroup() %>%
    mutate(prop=nn/n) %>%
    filter(class=="TE",X10==i) %>%
    select(-X6,-n,-nn, -X11, -class, -X10) %>% unique(.) %>%
    spread(key=tis_num,value = prop) %>%
    as.data.frame()
  rownames(c)<-c$X21
  c<-c[,-1]

  pheat_plot[[i]]<-pheatmap(as.matrix(c), scale = "row", cluster_cols = F, cluster_rows = F, 
                            show_rownames = T, main = i,legend_breaks = -1.5:1.5, legend = F,fontsize = 15,
                            color=colorRampPalette(c("#3d5488","white","#c01d20"))(50))
  #pheatmap(as.matrix(c), scale = "row", 
  #        cluster_cols = F, cluster_rows = F, show_rownames = T, main = i, breaks=c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01)),
  #         filename = paste(i,"pheat.pdf", sep="_"),
  #         legend_breaks = c(-2,-1.5, -1, -0.5, 0, 0.5, 1, 1.5,2),  legend_labels = c("-2", "-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5", "2")  )
}


P7<-gridExtra::grid.arrange(grobs= lapply( pheat_plot, function(x){ x$gtable } ), ncol= 3, labels=LETTERS[1:3] )
P7

pheat_plot<-list()

for(i in c("shared", "conserved", "private")){
  c<-
    TE_ACR_filter_tis %>%
    mutate(class=ifelse( X19==".", "Non-TE", "TE" ) ) %>%
    select(X6,X10,X11,X21,tis_num,class) %>%
    mutate_at(.vars = vars(tis_num), .funs = function(x){ paste(x,"tissues", sep = " ") } ) %>%
    group_by(X21,tis_num) %>%
    add_count()%>%
    group_by(X21,tis_num,X10,class)%>%
    add_count() %>%
    ungroup() %>%
    mutate(prop=nn/n) %>%
    filter(class=="TE",X10==i) %>%
    select(-X6,-n,-nn, -X11, -class, -X10) %>% unique(.) %>%
    spread(key=tis_num,value = prop) %>%
    as.data.frame()
  rownames(c)<-c$X21
  c<-c[,-1]
  
  pheat_plot[[i]]<-pheatmap(as.matrix(c), scale = "row", cluster_cols = F, cluster_rows = F, 
                            show_rownames = T, main = i,legend_breaks = -1.5:1.5, legend = F,fontsize = 15,
                            color=colorRampPalette(c("#3d5488","white","#c01d20"))(50))
  #pheatmap(as.matrix(c), scale = "row", 
  #        cluster_cols = F, cluster_rows = F, show_rownames = T, main = i, breaks=c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01)),
  #         filename = paste(i,"pheat.pdf", sep="_"),
  #         legend_breaks = c(-2,-1.5, -1, -0.5, 0, 0.5, 1, 1.5,2),  legend_labels = c("-2", "-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5", "2")  )
}


P8<-gridExtra::grid.arrange(grobs= lapply( pheat_plot, function(x){ x$gtable } ), ncol= 3, labels=LETTERS[1:3] )
P8



TE_ACR_final %>%
  #filter(X21=="flower") %>%
  filter(X10=="private") %>%
  filter(X12!=".") %>%
  select(X16) %>%
  group_by(X16) %>%
  add_count() %>% 
  arrange(desc(n)) %>%unique(.)


(TE_ACR_final %>% filter(X16=="TE_00018605_INT") )[,10:21 ] %>% filter(X10=="private") %>%
  select(X17) %>% unique(.)
# color followed by genome research of diploid wheat
#colorRampPalette(c("#3d5488","white","#c01d20"))(50)
#c01d20, #3d5488



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
  
  
  
