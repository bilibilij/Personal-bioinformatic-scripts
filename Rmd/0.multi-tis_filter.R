setwd("F:/工作/3.诸葛菜/00.diploid/27.TE_ACR")

library(tidyverse)

tis_bed<-data.frame()
for (i in dir("input-data/1.tis_ACR_filter") %>% grep("ine.ts$", . ,value = T) ){
  
  tis_bed  <-rbind(tis_bed,
                  read_tsv( paste("input-data/1.tis_ACR_filter", i , sep="/"), col_names = F ) %>%
                    mutate(class= gsub( ".geneDistance.closest.ACR.*"  ,"",
                                        gsub("OV_ATAC_", "",i) )  ) %>% 
                    mutate_at( .vars = vars(class), function(x){ gsub( "Seed", "seed", gsub( "^S$" ,"stem",  gsub( "R", "root",  gsub( "L", "leaf",  gsub( "F", "flower", x ) ) )) )      } ) %>%
                    mutate_at( .vars = vars(X33), function(x){ gsub( "Seed", "seed", gsub( "^S$" ,"stem",  gsub( "R", "root",  gsub( "L", "leaf",  gsub( "F", "flower", x ) ) )) )      } )
  )
}



tis_out<-
tis_bed %>%
  select(X6,class,X33) %>% unique(.) %>%
  mutate_at(.vars = vars(X33), .funs = function(x){ifelse(x==".", "tis_spe", x)}) %>%
  mutate( values=ifelse(X33=="tis_spe",0,1) ) %>%
  group_by(X6,class) %>%
  mutate( tis_num= sum(values)+1 ) %>%
  spread(key = X33, value = values) %>%
  mutate_at(.vars=vars(tis_spe), .funs = function(x){ gsub("0","1",x) } ) %>%
  ungroup() 

TE_ACR_final<-read_tsv("input-data/core_table/0.TE-ACR.tsv",col_names = F  )

TE_ACR_filter_tis<-left_join( TE_ACR_final, tis_out, by=c("X21"="class", "X6") )
  
write_tsv(TE_ACR_filter_tis, "input-data/core_table/1.TE-ACR-tis.tsv",col_names = T,quote = "none" )

