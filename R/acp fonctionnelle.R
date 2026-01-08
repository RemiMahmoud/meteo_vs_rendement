library(tidyverse)
library(fda)
library(RColorBrewer)
library(patchwork)
library(FactoMineR)
require(splines)
#Pour la régression fonctionnelle:
#https://cran.r-project.org/web/packages/robflreg/robflreg.pdf
library(robflreg)
library(tidyr)
library(dplyr)
library(ggplot2)

#meteo.qmd --> dt_meteo

dt_acpf <- dt_meteo  |> 
  mutate(
    m = month(DATE),
    d = day(DATE),
    fake_year = if_else(
      (m > 8) | (m == 8 & d >= 15), 1999, 2000),
    date_no_year = as.Date(sprintf("%04d-%02d-%02d", fake_year, m, d)))

#date du premier semis
premier_semis <- dt_acpf %>%
  group_by(exp) %>%
  summarise(premiere_date_semis = min(date_no_year, na.rm = TRUE))

date_ref <- as.Date(min(premier_semis$premiere_date_semis))
dt_acpf$compteur_jours <- as.numeric(dt_acpf$date_no_year - date_ref) + 1

summary(dt_acpf)

data_TEMP_T_Q <- dt_acpf %>%
  dplyr :: select(exp, compteur_jours, TEMP_T_Q) %>%
  tidyr :: pivot_wider(names_from = compteur_jours, values_from = TEMP_T_Q)

exp_col <- data_TEMP_T_Q["exp"]
num_cols <- data_TEMP_T_Q[, sapply(data_TEMP_T_Q, is.numeric)]
num_cols <- num_cols[, order(as.numeric(names(num_cols)))]

#Intervale avec aucun NAs
#Erreur dans eigen(Cmat) : valeurs infinies ou manquantes dans 'x'
#vient du fait que la matrice de covariance calculée à partir de ton objet fonctionnel (fdobj) contient des NA
good_cols <- apply(num_cols, 2, function(x) all(!is.na(x)))
num_cols_clean <- num_cols[, good_cols]
data_TEMP_T_Q <- cbind(exp_col, num_cols_clean)

matplot(
  t(data_TEMP_T_Q[, -1]),      # on transpose pour avoir les jours en x
  type = "b",                  # 'b' pour points et lignes
  pch = 19,                     # type de point
  col = 1:20,                  # couleurs différentes pour 20 parcelles
  xlab = "Jour",
  ylab = "Valeur",
  main = "Courbes des 20 parcelles"
)

jours <- as.numeric(names(num_cols_clean))
rangeval <- c(min(jours), max(jours)) 
nbasis <- 30 # Nombre de splines utilisées, superposer les approximations et nombre de plines  
basisobj <- create.bspline.basis(rangeval = rangeval, nbasis = nbasis) # création d'une base d'objet

# convertir les données discrètes en fd (estimation par régression sur la base)
fdsmooth <- smooth.basis(argvals = jours, y = t(as.matrix(num_cols_clean)), fdParobj = basisobj)
# pour y les lignes doivent être les valeurs discrètes et les colonnes les individus

fdobj <- fdsmooth$fd

# ACP fonctionnelle avec lissage medium, on peut avoir au maximum 20 variables
res.fpca <- pca.fd(fdobj, nharm = 5, harmfdPar = fdPar(fdobj = fdobj, lambda = 10))
plot.pca.fd(res.fpca, pointplot = FALSE)

res.fpca$harmonics$coefs
#Les harmoniques sont les dimensions. Une harmonique = une fonction
# Espace fonctionnel de dimension 50 
# PC1 et PC2 sont les 2 harmoniques, et les lignes sont les coefficients des 50 splines qui servent à recontruire les fonctions 

res.fpca$meanfd$coefs
#coefficients de la base (ici splines) pour la courbe moyenne

# Affichage du nuage des individus
plot(res.fpca$scores, type = "p")

nuage_individus <- data.frame(
  CP1 = res.fpca$scores[,1],
  CP2 = res.fpca$scores[,2],
  exp = data_TEMP_T_Q$exp)

rendement <- nuage_individus %>% left_join(data_rendement, by = "exp")

ggplot(rendement %>% filter(type == "SY", variety == "BRA_008_Aviso"),
  aes(CP1, CP2, label = exp, size = value)) +
  geom_point(alpha = 0.7) +
  geom_text(size = 5, vjust = -0.5) +
  scale_size_continuous(range = c(1, 10)) +
  labs(title = "Nuage des individus – ACP fonctionnelle", size = "Rendement SY") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.4)


# On essaye de lisser les courbes pour lisser les tendances
# Avec lissage fort 
res.fpca <- pca.fd(fdobj, harmfdPar = fdPar(fdobj = basisobj,lambda = 1e6))
p1 <- eval.fd(jours, res.fpca$harmonics) %>% 
  as_tibble() %>% 
  mutate(jours = jours) %>% 
  pivot_longer(-jours, names_to = "component") %>% 
  ggplot(aes(x =  jours, y=  value)) +
  geom_line(aes(group = component, color = component), linewidth = 1) +
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Avec lissage")+
  theme(plot.title  = element_text(hjust = 0.5))
#transforme des fonctions abstraites en valeurs numériques exploitables 

#Sans lissage
res.fpca <- pca.fd(fdobj)
p2 <- eval.fd(jours, res.fpca$harmonics) %>% 
  as_tibble() %>% 
  mutate(jours = jours) %>% 
  pivot_longer(-jours, names_to = "component") %>% 
  ggplot(aes(x =  jours, y=  value)) +
  geom_line(aes(group = component, color = component), linewidth = 1) +
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Sans lissage")+
  theme(plot.title = element_text(hjust = 0.5))

p1 / p2


# Affichage des 2 premières dimensions AVEC ou SANS lissages
eval.fd(jours, res.fpca$harmonics) %>% 
  as_tibble() %>% 
  mutate(jours = jours) %>% 
  pivot_longer(-jours, names_to = "component") %>% 
  ggplot(aes(x =  jours, y=  value)) +
  geom_line(aes(group = component, color = component), linewidth = 1) +
  scale_color_brewer(palette = "Dark2")

# Affichage des 2 premières dimensions avec la moyenne
eval.fd(jours, res.fpca$harmonics) %>% 
  as_tibble() %>% 
  mutate(jours = wavelength, mean = eval.fd(wavelength, res.fpca$meanfd)) %>% 
  pivot_longer(-wavelength, names_to = "component") %>% 
  ggplot(aes(x =  wavelength, y=  value)) +
  geom_line(aes(group = component, color = component), linewidth = 1) +
  scale_color_brewer(palette = "Dark2")

# Affichage de tous les spectres
plot(fdobj)
lines(res.fpca$meanfd, lwd = 1, col = "red")
plot(res.fpca$meanfd, lwd = 2, col = "red")
lines(eval.fd(wavelength, res.fpca$harmonics))

