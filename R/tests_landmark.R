# Jours de suivi pour deux parcelles (en jours après semis)
jours_parcelle1 <- 1:120
jours_parcelle2 <- 1:120

# NDVI simulé (croissance puis plateau puis déclin)
ndvi_parcelle1 <- c(
  seq(0.1, 0.9, length.out=40),  # croissance
  rep(0.9, 40),                  # plateau
  seq(0.9, 0.2, length.out=40)   # sénescence
)
ndvi_parcelle2 <- c(
  seq(0.1, 0.9, length.out=50),  # croissance plus lente
  rep(0.9, 30),                  # plateau plus court
  seq(0.9, 0.2, length.out=40)   # sénescence
)

# Dates des stades pour chaque parcelle
stades <- c("levée", "floraison", "maturité")
stades_p1 <- c(10, 60, 110)   
stades_p2 <- c(15, 70, 100)   




library(ggplot2)

# Timeline de référence : 0 (levée), 0.5 (floraison), 1 (maturité)
timeline_ref <- c(0, 0.5, 1)

align_with_landmarks <- function(jours, valeurs, stades_jours, timeline_ref, n_points=200) {
  # Crée la timeline recalée
  timeline_aligned <- seq(0, 1, length.out = n_points)
  
  # Fonction de mapping : des dates des stades vers timeline_ref
  map_fun <- approxfun(stades_jours, timeline_ref, rule=2)
  
  # Pour chaque point original, on calcule son temps recalé
  recal_time <- map_fun(jours)
  
  # print("ok")
  # Interpoler la variable sur la timeline recalée
  approx(x = recal_time, y = valeurs, xout = timeline_aligned, rule=2)$y
}

ndvi_p1_aligned <- align_with_landmarks(jours_parcelle1, ndvi_parcelle1, stades_p1, timeline_ref)
ndvi_p2_aligned <- align_with_landmarks(jours_parcelle2, ndvi_parcelle2, stades_p2, timeline_ref)




df_plot <- data.frame(
  Timeline = rep(seq(0,1,length.out=200), 2),
  NDVI = c(ndvi_p1_aligned, ndvi_p2_aligned),
  Parcelle = rep(c("Parcelle 1", "Parcelle 2"), each=200)
)

ggplot(df_plot, aes(Timeline, NDVI, color=Parcelle)) +
  geom_line(linewidth=1.2) +
  scale_x_continuous(breaks=timeline_ref, labels=stades) +
  labs(title="Alignement par landmark registration",
       x="Timeline phénologique recalée",
       y="NDVI") +
  theme_minimal()



df_plot <- data.frame(
  Timeline = c(jours_parcelle1, jours_parcelle2),
  NDVI = c(ndvi_parcelle1, ndvi_parcelle2),
  Parcelle = rep(c("Parcelle 1", "Parcelle 2"), each=120)
)

ggplot(df_plot, aes(Timeline, NDVI, color=Parcelle)) +
  geom_line(linewidth=1.2) +
  # scale_x_continuous(breaks=timeline_ref, labels=stades) +
  labs(title="Alignement par landmark registration",
       x="Timeline phénologique recalée",
       y="NDVI") +
  theme_minimal()
