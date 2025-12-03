# Jours de suivi pour deux parcelles (en jours après semis)
jours_parcelle1 <- 1:120
jours_parcelle2 <- 1:120

# NDVI simulé
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
  timeline_aligned <- seq(0, 1, length.out = n_points)
  
  map_fun <- approxfun(stades_jours, timeline_ref, rule=2)
  
  # Pour chaque point original, on calcule son temps recalé
  recal_time <- map_fun(jours)
  
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








library(SimilarityMeasures)
y = 1
x <- seq(0, 3, 0.1)
s1 <- sin(x)
s2 <- sin(x + y)
path1 <- matrix(s1, 31)
path2 <- matrix(s2, 31)
DTW(path1, path2, 10)

plot(x, s1, type="l", col="blue", lwd=2,
     main="Sinus et Sinus décalé",
     ylab="Valeur", xlab="x")
lines(x, s2, col="red", lwd=2)
legend("topright", legend=c("s1 = sin(x)", "s2 = sin(x+0.3)"),
       col=c("blue","red"), lwd=2)







library(dtw)
demo(dtw)

?dtw
?plot.dtw

n = 20

idx<-seq(0,6.28,len=n);
query<-sin(idx)+runif(n)/10;

plot(x = 1:n, y = query)

## A cosine is for template; sin and cos are offset by 25 samples
template<-cos(idx)

## Find the best match with the canonical recursion formula
library(dtw);
alignment<-dtw(query,template,keep=TRUE);

## Display the warping curve, i.e. the alignment curve
plot(alignment,type="threeway")

## Align and plot with the Rabiner-Juang type VI-c unsmoothed recursion
plot(
  dtw(query,template,keep=TRUE,
      step=rabinerJuangStepPattern(6,"c")),
  type="twoway",offset=-2);

## See the recursion relation, as formula and diagram
rabinerJuangStepPattern(6,"c")
plot(rabinerJuangStepPattern(6,"c"))

alignment$stepPattern


## And much more!

# ca ca fonctionne bien mais ca n'a pas de landmark, possibilité de faire ca 
# par morceau. Ou alors si on a des fonction qui se ressemble peut etre que ca marche


idx<-seq(0,6.28,len=100)
query<-matrix(c(sin(idx)+runif(100)/10, sin(idx)+runif(100)/10), 2)

plot(x = 1:100, y = query)

## A cosine is for template; sin and cos are offset by 25 samples
template<-matrix(c(cos(idx), cos(idx)), 2)

## Find the best match with the canonical recursion formula
library(dtw)
alignment<-dtw(query,template,keep=TRUE)

## Display the warping curve, i.e. the alignment curve
plot(alignment,type="threeway")

## Align and plot with the Rabiner-Juang type VI-c unsmoothed recursion
plot(
  dtw(query,template,keep=TRUE,
      step=rabinerJuangStepPattern(6,"c")),
  type="twoway",offset=-2)

## See the recursion relation, as formula and diagram
rabinerJuangStepPattern(6,"c")
plot(rabinerJuangStepPattern(6,"c"))

## And much more!

plot(x = alignment$index2, y = alignment$index1)






t  <- c(0, 1, 3, 5)
tp <- c(0, 2, 2.5, 5)

# fonction qui va de t à tp
f <- splinefun(x = t, y = tp, method = "monoH.FC")

# creation des points de la fonction
x <- seq(min(t), max(t), length.out = 100)
values <- f(x)

# plot
plot(t, tp, col = "red", pch = 19, main = "Spline monotone (Fritsch–Carlson)")
lines(x, values, col = "blue", lwd = 2)


# Fonction inverse

f_inv_num <- function(y) {
  solve_one <- function(y0) {
    uniroot(
      f = function(x) f(x) - y0,
      interval = c(min(t), max(t)),
      tol = .Machine$double.eps^0.5
    )$root
  }
  sapply(y, solve_one)
}

# Test
vals_inv <- f_inv_num(x)

plot(tp, t, col = "red", pch = 19, main = "Spline monotone (Fritsch–Carlson)")
lines(x, vals_inv, col = "blue", lwd = 2)

# vérification
xx1 <- f(vals_inv)
plot(x, xx1)

xx2 <- f_inv_num(values)
plot(x, xx2)





