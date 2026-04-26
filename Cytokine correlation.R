
library("openxlsx")
library(doParallel)
library(foreach)
library(ppcor)

cytokine=as.data.frame(read.xlsx("Cytokine.xlsx",rowNames = T,check.names = F,sep.names = " ",sheet = 1))
rownames(cytokine)=cytokine$METAF

modify_rownames <- function(df) {
new_names <- sapply(rownames(df), function(x) {
v <- suppressWarnings(as.numeric(x))
if (!is.na(v)) return(as.character(v))
x
})
rownames(df) <- new_names
df
}
cytokine <- modify_rownames(cytokine)

cytokine=cytokine[,c(-1,-2,-3)]

metadata=as.data.frame(read.csv("meta.atc.filter.csv"))
rownames(metadata)=metadata$METAF

newtoold=as.data.frame(read.csv("motu3/SP_new_old.csv",row.names = 1))
rownames(newtoold)=newtoold$new

fdata=as.data.frame(read.xlsx("6803/RAdata_motu3_7562samples_newsp_20250629.xlsx",rowNames = T,check.names = F,sep.names = " "))
fdata =fdata[, colnames(fdata) != "Unassigned species", drop = FALSE]
fdata=fdata[rownames(metadata),]
fdata=fdata[,(colSums(fdata>0)/nrow(fdata))>0.025]

metadata=metadata[rownames(metadata) %in% rownames(cytokine),]
cytokine=cytokine[rownames(metadata),]

fdata=fdata[rownames(metadata),]
data=fdata

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

signif_table <- foreach(i = 1:(ncol(data)), .combine = rbind) %dopar% {
local_results <- data.frame()
for (j in 1:ncol(cytokine)) {
temp <- tryCatch({
tempdata=data.frame(cytokine=cytokine[,j],meta=data[,i])
result <- cor.test(tempdata[, 1], tempdata[, 2], method = "spearman",exact = T)
data.frame(n = result$estimate, p = result$p.value)
}, error = function(e) {
data.frame(n = NA, p = NA)
})

temp <- data.frame(
  Old = newtoold[colnames(data)[i],]$old,
  New = colnames(data)[i],
  Cytokine = colnames(cytokine)[j],
  Rho = temp$n,
  P = temp$p,
  n = nrow(tempdata)
)

local_results <- rbind(local_results, temp)

}
local_results
}

stopCluster(cl)
write.csv(signif_table,"cytokine_vs_sp_spearman_1740_20260210.csv",row.names = F)