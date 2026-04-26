library(openxlsx)
library(parallel)

metadata=as.data.frame(read.csv("meta.atc.filter.csv"))
metakey=colnames(metadata)[6:33]

r=NULL
for (i in 1:(length(metakey)-1)){
  for (j in (i+1):length(metakey)){
    if (i ==j){next}
    if(i ==length(metakey)){next}
    temp=na.omit(metadata[,c(metakey[i],metakey[j])])
    n=nrow(temp)
    temp=cor.test(temp[,1],temp[,2],method = "spearman",exact = T)
    tempdata=data.frame(
      feature1=metakey[i],
      feature2=metakey[j],
      rho=temp$estimate,
      p=temp$p.value,
      n=n
      
    )
    r=rbind(r,tempdata)
    
    
    
    
  }
}
write.csv(r,"metadata_spearman_20260120.csv",row.names = F)
