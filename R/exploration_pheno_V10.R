library(dplyr)
library(naniar)
library(FactoMineR)
library(ggplot2)

summary(dta_pheno_V10)
dta_pheno_V10$variety=as.factor(dta_pheno_V10$variety)
colSums(is.na(dta_pheno_V10)) #somme des NA par colonnes
vis_miss(dta_pheno_V10)

data_SN  <- dta_pheno_V10 %>% select(contains("SN"))
data_SY  <- dta_pheno_V10 %>% select(contains("SY"))
data_TSW <- dta_pheno_V10 %>% select(contains("TSW"))

data_SN_E <- dta_pheno_V10 %>% select(matches("variety|SN_E1|SN_E2|SN_E3|SN_E4"))
data_SY_E <- dta_pheno_V10 %>% select(matches("variety|SY_E1|SY_E2|SY_E3|SY_E4"))
data_TSW_E <- dta_pheno_V10 %>% select(matches("variety|TSW_E1|TSW_E2|TSW_E3|TSW_E4"))

SN_cols <- setdiff(names(data_SN), names(data_SN_E))
SY_cols <- setdiff(names(data_SY), names(data_SY_E))
TSW_cols <- setdiff(names(data_TSW), names(data_TSW_E))

data_SN_P <- dta_pheno_V10 %>% select(all_of(SN_cols), variety)
data_SY_P <- dta_pheno_V10 %>% select(all_of(SY_cols), variety)
data_TSW_P <- dta_pheno_V10 %>% select(all_of(TSW_cols), variety)
data_SY_P <- data_SY_P[,0:23]

names(data_SN_P) <- gsub("^SN_", "", names(data_SN_P))
names(data_SY_P) <- gsub("^SY_", "", names(data_SY_P))
names(data_TSW_P) <- gsub("^TSW_", "", names(data_TSW_P))

names_SN_P <- names(data_SN_P)
names_SY_P <- names(data_SY_P)
names_TSW_P <- names(data_TSW_P)

common_all <- Reduce(intersect, list(names_SN_P, names_SY_P, names_TSW_P))
common_all
#Les parcelles en communs dans les 3 jeux de données 

#Graphique pour MT et SN / SY / TSW
data_SN_P %>% ggplot() + aes(x=variety, y=MT) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_SY_P %>% ggplot() + aes(x=variety, y=MT) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
data_TSW_P %>% ggplot() + aes(x=variety, y=MT) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
