library(vegan)
library(openxlsx)
library(Maaslin2)
library(ppcor)

metadata=as.data.frame(read.csv("meta.atc.csv"))
metadata=metadata[,c(-3,-4)]
metadata=metadata[!is.na(metadata$METAF),]
metadata=metadata[metadata$METAF %in% rownames(fdata),]

fixkey <- colnames(metadata)[c(3:5,51:55)]
diskey <- colnames(metadata)[6:50]
rownames(metadata)=metadata$METAF

fdata=as.data.frame(read.xlsx("RAdata_motu3_7562samples_newsp_20250629.xlsx",rowNames = T,check.names = F,sep.names = " "))
fdata=fdata[rownames(metadata),]
fdata =fdata[, colnames(fdata) != "Unassigned species", drop = FALSE]
fdata <- fdata / rowSums(fdata)
fdata=fdata[,(colSums(fdata>0)/nrow(fdata))>0.025]

ml=read.table("Japanese4D_feces_load_7562samples_20250911.tsv",sep="\t",header = 1,row.names = 1)
ml=modify_rownames(ml)

shannon <- diversity(fdata, index = "shannon")
richness <- rowSums(fdata > 0)
ml=ml[rownames(fdata),,drop=F]
alpha_df <- data.frame(
SampleID = rownames(fdata),
Shannon = shannon,
Richness = richness,
Microbialload=ml
)

habit_cols <- colnames(metadata)[6:33]
disease_cols <- colnames(metadata)[34:50]
test_cols <- c(habit_cols, disease_cols)

results <- list()

for (var in test_cols) {
temp_meta <- metadata[, c(var,fixkey), drop = FALSE]
temp_meta$SampleID <- rownames(temp_meta)

merged <- merge(alpha_df, temp_meta, by = "SampleID")
merged <- merged[!is.na(merged[[var]]), ]

cor1 <- pcor.test(merged$Shannon, merged[[var]],merged[,c(fixkey)],method = "spearman")
cor2 <- pcor.test(merged$Richness, merged[[var]],merged[,c(fixkey)],method = "spearman")
cor3 <- pcor.test(merged$load, merged[[var]],merged[,c(fixkey)],method = "spearman")

results[[var]] <- data.frame(
Variable = var,
Spearman_Shannon = cor1$estimate,
P_Shannon = cor1$p.value,
Spearman_Richness = cor2$estimate,
P_Richness = cor2$p.value,
Spearman_ML = cor2$estimate,
P_ML = cor2$p.value,
N = nrow(merged)
)
}

res_df <- do.call(rbind, results)
write.csv(res_df,"./6803_motusp_8fix_alpha.csv",row.names = F)