library(vegan)
library(openxlsx)
library(Maaslin2)
library(doParallel)

metadata=as.data.frame(read.csv("meta.atc.filter.csv"))

fixkey <- colnames(metadata)[c(3:5,51:55)]
diskey <- colnames(metadata)[6:50]

rownames(metadata)=metadata$METAF

fdata=as.data.frame(read.csv("mo_7562samples_20250428.txt",sep="\t",check.names=F,row.names=1))
colnames(fdata)=c(colnames(fdata)[2:ncol(fdata)],"test")
fdata=fdata[,-ncol(fdata)]
fdata=as.data.frame(t(fdata))

modify_rownames <- function(df) {
old_names <- rownames(df)
new_names <- sapply(old_names, function(name) {
num_name <- suppressWarnings(as.numeric(name))
if (!is.na(num_name)) {
return(as.character(num_name))
} else {
return(name)
}
})
rownames(df) <- new_names
return(df)
}

fdata=modify_rownames(fdata)
metadata=modify_rownames(metadata)

metadata=metadata[rownames(metadata) %in% rownames(fdata),]
fdata=fdata[rownames(metadata),]
fdata=fdata[,(colSums(fdata>0)/nrow(fdata))>0.025]

write.csv(fdata,"6803_MO_20260117.csv")

data=fdata
spname=data.frame(feature=paste0("X",1:ncol(data)),name=colnames(data))
colnames(data)=spname$feature
rownames(spname)=spname$feature

pr=NULL
ncores <- 10
cl <- makeCluster(ncores)
registerDoParallel(cl)

parallel_results <- foreach(i = diskey,
.packages = c("Maaslin2","dplyr")) %dopar% {
set.seed(123)

tempmeta=metadata
tempmeta=tempmeta[,c(i,fixkey)]
tempmeta=na.omit(tempmeta)
colnames(tempmeta)[1]="test"

binkey <- sapply(colnames(tempmeta), function(k) {
col_vals <- tempmeta[[k]]
all(col_vals %in% c(0,1)) && any(col_vals==0) && any(col_vals==1)
})

binkey=binkey[binkey==TRUE]
for (j in names(binkey)){
tempmeta[,j]=factor(tempmeta[,j],levels=c(0,1))
}

temp <- Maaslin2(
data[rownames(tempmeta),],
tempmeta,
standardize=TRUE,
output="output",
fixed_effects=c("test",fixkey),
min_abundance=0,
min_prevalence=0,
normalization="NONE",
plot_heatmap=FALSE,
plot_scatter=FALSE
)

tempp=temp$results
tempp=tempp[tempp$metadata %in% c("test"),]
tempp$metadata=gsub("test",i,tempp$metadata)

tempp=merge(spname,tempp,by="feature")
tempp=tempp[,-1]
tempp$fix=paste(fixkey,collapse=",")

colnames(tempp)[1]="feature"
tempp
}

pr <- do.call(rbind, parallel_results)
stopCluster(cl)

pr=pr[,c(-3,-5,-7)]
write.csv(pr,"./6803_MO_MAA_8F_20260115.csv",row.names=F)


new_old_name=as.data.frame(read.csv("SP_new_old.csv",row.names=1))
rownames(new_old_name)=new_old_name$new

fdata=as.data.frame(read.xlsx("RAdata_motu3_7562samples_newsp_20250629.xlsx",rowNames=T,check.names=F,sep.names=" "))
fdata=fdata[rownames(metadata),]
fdata=fdata[,(colSums(fdata>0)/nrow(fdata))>0.025]

write.csv(fdata,"6803_motusp_20260117.csv")

data=fdata
data=data[,-1]

spname=data.frame(feature=paste0("X",1:ncol(data)),name=colnames(data))
colnames(data)=spname$feature
rownames(spname)=spname$feature

pr=NULL
ncores <- 10
cl <- makeCluster(ncores)
registerDoParallel(cl)

parallel_results <- foreach(i = diskey,
.packages=c("Maaslin2","dplyr")) %dopar% {
set.seed(123)

tempmeta=metadata
tempmeta=tempmeta[,c(i,fixkey)]
tempmeta=na.omit(tempmeta)
colnames(tempmeta)[1]="test"

binkey <- sapply(colnames(tempmeta), function(k){
col_vals <- tempmeta[[k]]
all(col_vals %in% c(0,1)) && any(col_vals==0) && any(col_vals==1)
})

binkey=binkey[binkey==TRUE]
for (j in names(binkey)){
tempmeta[,j]=factor(tempmeta[,j],levels=c(0,1))
}

temp <- Maaslin2(
data[rownames(tempmeta),],
tempmeta,
standardize=TRUE,
output="output",
fixed_effects=c("test",fixkey),
min_abundance=0,
min_prevalence=0,
normalization="NONE",
plot_heatmap=FALSE,
plot_scatter=FALSE
)

tempp=temp$results
tempp=tempp[tempp$metadata=="test",]
tempp$metadata=i

tempp=merge(spname,tempp,by="feature")
tempp=tempp[,-1]
tempp$fix=paste(fixkey,collapse=",")

colnames(tempp)[1]="feature"
tempp$old=new_old_name[tempp$feature,]$old
tempp
}

stopCluster(cl)
pr <- do.call(rbind, parallel_results)
pr=pr[,c(-3,-5,-7)]
pr=unique(pr)
write.csv(pr,"./6803_motusp_MAA_8F_20260115.csv",row.names=F)



fdata=as.data.frame(read.xlsx("ko_7562samples_20250428.xlsx",rowNames=T,check.names=F,sep.names=" "))
fdata=as.data.frame(t(fdata))

fdata=fdata[rownames(metadata),]
fdata=fdata[,(colSums(fdata>0)/nrow(fdata))>0.025]

write.csv(fdata,"6803_KO_20260117.csv")
saveRDS(fdata,"6803_KO_20260117.rds")

data=fdata
spname=data.frame(feature=paste0("X",1:ncol(data)),name=colnames(data))
colnames(data)=spname$feature
rownames(spname)=spname$feature

pr=NULL
ncores <- 10
cl <- makeCluster(ncores)
registerDoParallel(cl)

parallel_results <- foreach(i = diskey,
.packages=c("Maaslin2","dplyr")) %dopar% {
set.seed(123)

tempmeta=metadata
tempmeta=tempmeta[,c(i,fixkey)]
tempmeta=na.omit(tempmeta)
colnames(tempmeta)[1]="test"

binkey <- sapply(colnames(tempmeta), function(k){
col_vals <- tempmeta[[k]]
all(col_vals %in% c(0,1)) && any(col_vals==0) && any(col_vals==1)
})

binkey=binkey[binkey==TRUE]
for (j in names(binkey)){
tempmeta[,j]=factor(tempmeta[,j],levels=c(0,1))
}

temp <- Maaslin2(
data[rownames(tempmeta),],
tempmeta,
standardize=TRUE,
output="output",
fixed_effects=c("test",fixkey),
min_abundance=0,
min_prevalence=0,
normalization="NONE",
plot_heatmap=FALSE,
plot_scatter=FALSE
)

tempp=temp$results
tempp=tempp[tempp$metadata %in% c("test"),]
tempp$metadata=gsub("test",i,tempp$metadata)

tempp=merge(spname,tempp,by="feature")
tempp=tempp[,-1]
tempp$fix=paste(fixkey,collapse=",")

colnames(tempp)[1]="feature"
tempp
}

stopCluster(cl)
pr <- do.call(rbind, parallel_results)
pr=pr[,c(-3,-5,-7)]
pr=unique(pr)
write.csv(pr,"./6803_KO_MAA_8F_20260115.csv",row.names=F)
