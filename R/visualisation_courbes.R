# Visualiser les 20 courbes temporelles de chaque exp en fonction des variables de ind_dev (les variables par période climatiques)
# On utilise le jeu de donnée dt du fichier r "crea_dta_clean"

library(tidyverse)
library(ggplot2)

# formatage DATE
dt_courbe <- dt %>%
  mutate(DATE = as_date(DATE, origin = "1899-12-30"))
dt_courbe$ANNEE <- year(dt_courbe$DATE)
dt_courbe$JOUR_MOIS <- format(dt_courbe$DATE, "%d-%m") 

summary(dt_courbe)

dt_courbe$exp=as.factor(dt_courbe$exp)
dt_courbe$JOUR_MOIS <- as.Date(paste0("2025/", dt_courbe$JOUR_MOIS), format="%Y/%d/%m")
 
