library(openxlsx)
library(clusterProfiler)
library(tidyr)
library(tidyverse)
mas_all_metadata=as.data.frame(read.csv("./6803_KO_MAA_8F_20260115.csv"))
key=unique(mas_all_metadata$metadata)
key1=key[1:28]

result=NULL
for (meta in unique(mas_all_metadata$metadata)){
  KO_up <- mas_all_metadata %>% filter(metadata==meta, coef>0, qval<0.1) %>% pull(feature)
  KO_down <- mas_all_metadata %>% filter(metadata==meta, coef<0, qval<0.1) %>% pull(feature)
  
  ko_up <- try(
    enrichKEGG(gene = KO_up,
               organism = "ko",
               keyType = "kegg",
               pAdjustMethod = "BH",
               pvalueCutoff = 1) %>%
      .@result %>% as_tibble() %>%
      mutate(metadata=meta, KO_list="up_regulated",Factor_category=ifelse(meta %in%key1,"Lifestyle","HealthDisease"),.before=1)
  )
  
  ko_down <- try(
    enrichKEGG(gene = KO_down,
               organism = "ko",
               keyType = "kegg",
               pAdjustMethod = "BH",
               pvalueCutoff = 1) %>%
      .@result %>% as_tibble() %>%
      mutate(metadata=meta, KO_list="down_regulated",Factor_category=ifelse(meta %in%key1,"Lifestyle","HealthDisease"), .before=1)
  )
  result=rbind(result,ko_up,ko_down)
}
write.csv(result,"KEGGID_enrich_20260115.csv")
