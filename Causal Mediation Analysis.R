library(mediation)
library(readxl)
library(CMAverse)
library(openxlsx)
library(tidyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

key_idx <- as.numeric(args[1])
dis_idx <- as.numeric(args[2])

output_file <- file.path("result", paste0("6803_causal_motu_", key_idx, "_", dis_idx, ".csv"))

meta = as.data.frame(read.csv("meta.atc.filter.csv", check.names = FALSE))
colnames(meta)[c(3:5,51:55)] = paste0("ADJ",1:8)

dis_all = colnames(meta)[c(34:41,43:50)]
key_all = colnames(meta)[6:33]

key = key_all[key_idx]
dis = dis_all[dis_idx]

colnames(meta)[2] = "METAF"
meta = meta[!is.na(meta$METAF),]
rownames(meta) = meta$METAF

fdata = as.data.frame(read.csv("6803_motusp_20260117.csv",
check.names = FALSE,
row.names = 1))

fdata[fdata == 0] = min(fdata[fdata != 0]) / 2
fdata = as.data.frame(log10(fdata))

maa = read.csv("6803_motusp_MAA_8F_20260115.csv")



food_maaslin = maa[maa$metadata == key,]
food_maaslin = food_maaslin[food_maaslin$qval < 0.1,]

dis_maaslin = maa[maa$metadata == dis,]
dis_maaslin = dis_maaslin[dis_maaslin$qval < 0.1,]

data = fdata[, colnames(fdata) %in%
intersect(food_maaslin$feature, dis_maaslin$feature),
drop = FALSE]

sp_name = colnames(data)
colnames(data) = 1:ncol(data)

cat("Feature number:", ncol(data), "\n")

result_list <- list()

for (num in 1:ncol(data)) {

sp_analyzed = sp_name[num]

temp_data = cbind(data[,num],
meta[,c(key, dis, paste0("ADJ",1:8))])

temp_data = as.data.frame(sapply(temp_data, as.numeric))
temp_data = na.omit(temp_data)

colnames(temp_data)[1:2] = c("medi","treat")
colnames(temp_data)[3] = "outcome"

set.seed(123)
model.m <- lm(medi ~ treat + ADJ1 + ADJ2 + ADJ3 + ADJ4 +
ADJ5 + ADJ6 + ADJ7 + ADJ8,
data = temp_data)

model.y <- glm(outcome ~ treat + medi + ADJ1 + ADJ2 + ADJ3 +
ADJ4 + ADJ5 + ADJ6 + ADJ7 + ADJ8,
data = temp_data,
family = binomial(link="logit"))

cont1 <- try(mediate(model.m, model.y,
sims = 1000,
treat = "treat",
mediator = "medi"))

set.seed(123)
model.m.inv <- glm(outcome ~ treat + ADJ1 + ADJ2 + ADJ3 +
ADJ4 + ADJ5 + ADJ6 +ADJ7 + ADJ8,
data = temp_data,
family = binomial(link="logit"))

model.y.inv <- lm(medi ~ outcome + treat + ADJ1 + ADJ2 +
ADJ3 + ADJ4 + ADJ5 + ADJ6 + ADJ7 +ADJ8,
data = temp_data)

cont2 <- try(mediate(model.m.inv, model.y.inv,
sims = 1000,
treat = "treat",
mediator = "outcome"))

result_temp <- data.frame(
Treat = key,
Disease = dis,
Mediater = sp_analyzed,
N = nrow(temp_data),

Avg_causal_mediation_effect = ifelse(class(cont1)!="try-error", cont1$d.avg, NA),
Avg_causal_mediation_pval   = ifelse(class(cont1)!="try-error", cont1$d.avg.p, NA),
Avg_causal_direct_effect    = ifelse(class(cont1)!="try-error", cont1$z.avg, NA),
Avg_causal_direct_pval      = ifelse(class(cont1)!="try-error", cont1$z.avg.p, NA),
Prop_mediated               = ifelse(class(cont1)!="try-error", cont1$n.avg, NA),
Prop_mediated_pval          = ifelse(class(cont1)!="try-error", cont1$n.avg.p, NA),
Avg_total_effect            = ifelse(class(cont1)!="try-error", cont1$tau.coef, NA),
Avg_total_pval              = ifelse(class(cont1)!="try-error", cont1$tau.p, NA),

Rev_Avg_causal_mediation_effect = ifelse(class(cont2)!="try-error", cont2$d.avg, NA),
Rev_Avg_causal_mediation_pval   = ifelse(class(cont2)!="try-error", cont2$d.avg.p, NA),
Rev_Avg_causal_direct_effect    = ifelse(class(cont2)!="try-error", cont2$z.avg, NA),
Rev_Avg_causal_direct_pval      = ifelse(class(cont2)!="try-error", cont2$z.avg.p, NA),
Rev_Prop_mediated               = ifelse(class(cont2)!="try-error", cont2$n.avg, NA),
Rev_Prop_mediated_pval          = ifelse(class(cont2)!="try-error", cont2$n.avg.p, NA),
Rev_Avg_total_effect            = ifelse(class(cont2)!="try-error", cont2$tau.coef, NA),
Rev_Avg_total_pval              = ifelse(class(cont2)!="try-error", cont2$tau.p, NA),

EXP_coef = food_maaslin[food_maaslin$feature==sp_analyzed,"coef"],
EXP_pval = food_maaslin[food_maaslin$feature==sp_analyzed,"pval"],
EXP_qval = food_maaslin[food_maaslin$feature==sp_analyzed,"qval"],

Dis_coef = dis_maaslin[dis_maaslin$feature==sp_analyzed,"coef"],
Dis_pval = dis_maaslin[dis_maaslin$feature==sp_analyzed,"pval"],
Dis_qval = dis_maaslin[dis_maaslin$feature==sp_analyzed,"qval"]

)

result_list[[num]] <- result_temp
}

result <- do.call(rbind, result_list)

result$Avg_causal_mediation_qval =
p.adjust(result$Avg_causal_mediation_pval, "fdr")

write.csv(result, output_file, row.names = FALSE)
