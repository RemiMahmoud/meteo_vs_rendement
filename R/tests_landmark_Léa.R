# Installer si nécessaire
install.packages("fda")
library(fda)
library(ggplot2)

# Générer 100 points pour t de 0 à 1
t <- seq(0, 1, length.out = 100)

# Simuler 5 courbes sinusoïdales avec un décalage
set.seed(123)
curves <- sapply(1:5, function(i) sin(2 * pi * t + runif(1, -0.2, 0.2)))

# Visualiser les courbes
matplot(t, curves, type = "l", lty = 1, col = 1:5,
        xlab = "Temps", ylab = "Valeur",
        main = "Courbes originales")

# Extraire les landmarks (position du maximum)
landmarks_idx <- apply(curves, 2, which.max)
landmarks <- matrix(t[landmarks_idx], nrow = 5, ncol = 1)
landmarks

# x0marks = moyenne des landmarks
x0marks <- colMeans(landmarks)
x0marks

# X0lim pour landmarkreg
x0lim <- c(0,1)

# -------------------------------
# Créer fd_obj pour les courbes via smooth.morph()
# -------------------------------

fd_list <- vector("list", ncol(curves))
for(i in 1:ncol(curves)) {
  y <- curves[,i]
  
  # Base initiale pour smooth.morph
  nbasis <- 7
  basisW <- create.fourier.basis(rangeval = c(0,1), nbasis = nbasis)
  Wfd0 <- fd(coef = rep(0, nbasis), basisobj = basisW)
  WfdPar0 <- fdPar(Wfd0, Lfdobj = 2, lambda = 1e-2)
  
  # smooth.morph : renvoie un monfd
  mono <- smooth.morph(x = t, y = y, WfdPar = WfdPar0, ylim = range(y))
  
  # récupérer l'objet fd et mettre dans la liste
  fd_list[[i]] <- mono$fd
}

# Vérifier que fd_list n'est pas vide
length(fd_list)  # doit être 5

# Créer fd_obj combiné
coef_matrix <- sapply(fd_list, coef)
fd_obj <- fd(coef_matrix, fd_list[[1]]$basis)

# -------------------------------
# Créer WfdPar pour landmarkreg
# -------------------------------
# On peut réutiliser la base utilisée pour smooth.morph
WfdPar <- fdPar(fd_list[[1]]$basis, Lfdobj = 2, lambda = 1e-2)

# -------------------------------
# Landmark registration
# -------------------------------
result <- landmarkreg(
  unregfd = fd_obj,
  ximarks = landmarks,
  x0marks = x0marks,
  x0lim = x0lim,
  WfdPar = WfdPar,
  WfdPar0 = NULL,
  ylambda = 1e-10
)

# Visualiser les courbes alignées
matplot(t, eval.fd(t, result$fdreg),
        type = "l", col = 1:5, lty = 1,
        main = "Courbes alignées")

