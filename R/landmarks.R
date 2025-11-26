install.packages("ncdf4")
library(ncdf4)
library("raster")
dt = raster("data/rendement/soybean/yield_2009.nc4")
dt
plot(dt)


library(fda)

# 1. CRÉER DES DONNÉES SIMULÉES
# ==============================
set.seed(123)
t <- seq(0, 1, length.out = 100)  # Grille temporelle
n_curves <- 5  # Nombre de courbes

# Fonction de base : une courbe sinusoïdale avec un pic
base_curve <- sin(2 * pi * t) + 2 * exp(-((t - 0.5)^2) / 0.01)

# Créer plusieurs courbes avec des décalages temporels
curves <- matrix(NA, nrow = length(t), ncol = n_curves)
for (i in 1:n_curves) {
  shift <- (i - 3) * 0.05  # Décalage temporel
  t_shifted <- t - shift
  t_shifted <- pmax(0, pmin(1, t_shifted))  # Garder dans [0,1]
  curves[, i] <- approx(t, base_curve, xout = t_shifted, rule = 2)$y + rnorm(length(t), 0, 0.1)
}

# Visualisation des courbes originales
fda::matplot(t, curves, type = 'l', lty = 1, col = 1:n_curves, lwd = 2,
             main = "Courbes originales (désalignées)", xlab = "Temps", ylab = "Valeur")
legend("topright", legend = paste("Courbe", 1:n_curves), col = 1:n_curves, lty = 1, lwd = 2)


# 2. IDENTIFIER LES LANDMARKS
# ============================
# Trouver le maximum de chaque courbe (notre landmark)
landmarks <- apply(curves, 2, function(y) t[which.max(y)])
print("Landmarks (position du maximum) :")
print(landmarks)

# Choisir la position cible (par exemple, la médiane)
target_landmark <- median(landmarks)
print(paste("Landmark cible:", target_landmark))


# 3. CRÉER LES OBJETS FONCTIONNELS
# =================================
# Base B-spline pour les courbes
basis <- create.bspline.basis(rangeval = range(t), nbasis = 20, norder = 4)

# Convertir en objets fonctionnels
fdobj <- Data2fd(argvals = t, y = curves, basisobj = basis)


# 4. LANDMARK REGISTRATION
# =========================
# Ajouter les extrémités comme landmarks (nécessaire pour nbasis >= norder)
ximarks <- t(rbind(
  rep(min(t), n_curves),      # Début
  landmarks,                   # Maximum
  rep(max(t), n_curves)        # Fin
))

x0marks <- c(min(t), target_landmark, max(t))

# Base pour la fonction de warping
warp_basis <- create.bspline.basis(rangeval = range(t), nbasis = 5, norder = 4)
warp_fdPar <- fdPar(fdobj = warp_basis, Lfdobj = 2, lambda = 1e-6)

# Registration
reg <- landmarkreg(fdobj,
                   ximarks = ximarks,
                   x0marks = x0marks,
                   WfdPar = warp_fdPar)


# 5. VISUALISER LES RÉSULTATS
# ============================
# Extraire les courbes alignées
registered_curves <- eval.fd(t, reg$regfd)

# Graphique comparatif
par(mfrow = c(1, 2))

# Avant registration
matplot(t, curves, type = 'l', lty = 1, col = 1:n_curves, lwd = 2,
        main = "Avant registration", xlab = "Temps", ylab = "Valeur")
abline(v = landmarks, col = 1:n_curves, lty = 2)
abline(v = target_landmark, col = "red", lwd = 2, lty = 2)

# Après registration
matplot(t, registered_curves, type = 'l', lty = 1, col = 1:n_curves, lwd = 2,
        main = "Après registration", xlab = "Temps", ylab = "Valeur")
abline(v = target_landmark, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))


# 6. VISUALISER LES FONCTIONS DE WARPING
# =======================================
warp_functions <- eval.fd(t, reg$warpfd)
matplot(t, warp_functions, type = 'l', lty = 1, col = 1:n_curves, lwd = 2,
        main = "Fonctions de warping", xlab = "Temps original", ylab = "Temps aligné")
abline(0, 1, lty = 2, col = "gray")  # Ligne identité


# 7. VÉRIFIER L'ALIGNEMENT DES LANDMARKS
# =======================================
print("\nVérification - Position des maxima après registration:")
registered_landmarks <- apply(registered_curves, 2, function(y) t[which.max(y)])
print(registered_landmarks)
print(paste("Écart-type avant:", round(sd(landmarks), 4)))
print(paste("Écart-type après:", round(sd(registered_landmarks), 4)))











library(fda)

# 1. CRÉER DES DONNÉES SIMULÉES
# ==============================
set.seed(123)
t <- seq(0, 1, length.out = 100)  # Grille temporelle
n_curves <- 5  # Nombre de courbes

# Fonction de base : une courbe sinusoïdale avec un pic bien défini
base_curve <- sin(2 * pi * t) + 2 * exp(-((t - 0.5)^2) / 0.01)

# Créer plusieurs courbes avec des décalages temporels
curves <- matrix(NA, nrow = length(t), ncol = n_curves)
for (i in 1:n_curves) {
  shift <- (i - 3) * 0.05  # Décalage temporel
  t_shifted <- t - shift
  t_shifted <- pmax(0, pmin(1, t_shifted))  # Garder dans [0,1]
  curves[, i] <- approx(t, base_curve, xout = t_shifted, rule = 2)$y + rnorm(length(t), 0, 0.05)
}

# Visualisation des courbes originales
fda::matplot(t, curves, type = 'l', lty = 1, col = 1:n_curves, lwd = 2,
             main = "Courbes originales (désalignées)", xlab = "Temps", ylab = "Valeur")
legend("topright", legend = paste("Courbe", 1:n_curves), col = 1:n_curves, lty = 1, lwd = 2)


# 2. IDENTIFIER LES LANDMARKS (EN S'ASSURANT QU'ILS SONT DANS LA PLAGE)
# ======================================================================
# Trouver le maximum de chaque courbe, en excluant les 5 premiers et derniers points
landmarks <- apply(curves, 2, function(y) {
  idx <- which.max(y[6:(length(y)-5)]) + 5  # Chercher dans la zone centrale
  t[idx]
})

# S'assurer que les landmarks sont strictement dans la plage
landmarks <- pmax(t[2], pmin(t[length(t)-1], landmarks))

print("Landmarks (position du maximum) :")
print(landmarks)

# Choisir la position cible (par exemple, la médiane)
target_landmark <- median(landmarks)
print(paste("Landmark cible:", target_landmark))


# 3. CRÉER LES OBJETS FONCTIONNELS
# =================================
# Base B-spline pour les courbes
basis <- create.bspline.basis(rangeval = range(t), nbasis = 20, norder = 4)

# Convertir en objets fonctionnels
fdobj <- Data2fd(argvals = t, y = curves, basisobj = basis)


# 4. LANDMARK REGISTRATION
# =========================
# Structure : chaque LIGNE = un landmark, chaque COLONNE = une courbe
ximarks <- t(matrix(landmarks, nrow = 1, ncol = n_curves))

# Landmark cible
x0marks <- target_landmark

# Base pour la fonction de warping (simple avec peu de fonctions)
warp_basis <- create.bspline.basis(rangeval = range(t), nbasis = 4, norder = 4)
warp_fdPar <- fdPar(fdobj = warp_basis, Lfdobj = 2, lambda = 1e-4)

# Registration
reg <- landmarkreg(fdobj,
                   ximarks = ximarks,
                   x0marks = x0marks,
                   WfdPar = warp_fdPar)


# 5. VISUALISER LES RÉSULTATS
# ============================
# Extraire les courbes alignées
registered_curves <- eval.fd(t, reg$regfd)

# Graphique comparatif
par(mfrow = c(1, 2))

# Avant registration
matplot(t, curves, type = 'l', lty = 1, col = 1:n_curves, lwd = 2,
        main = "Avant registration", xlab = "Temps", ylab = "Valeur")
abline(v = landmarks, col = 1:n_curves, lty = 2, alpha = 0.5)
abline(v = target_landmark, col = "red", lwd = 3, lty = 2)
legend("topright", legend = "Cible", col = "red", lwd = 3, lty = 2)

# Après registration
matplot(t, registered_curves, type = 'l', lty = 1, col = 1:n_curves, lwd = 2,
        main = "Après registration", xlab = "Temps", ylab = "Valeur")
abline(v = target_landmark, col = "red", lwd = 3, lty = 2)

par(mfrow = c(1, 1))


# 6. VISUALISER LES FONCTIONS DE WARPING
# =======================================
warp_functions <- eval.fd(t, reg$warpfd)
matplot(t, warp_functions, type = 'l', lty = 1, col = 1:n_curves, lwd = 2,
        main = "Fonctions de warping", xlab = "Temps original", ylab = "Temps aligné")
abline(0, 1, lty = 2, col = "gray", lwd = 2)  # Ligne identité
legend("topleft", legend = "Identité", col = "gray", lwd = 2, lty = 2)


# 7. VÉRIFIER L'ALIGNEMENT DES LANDMARKS
# =======================================
print("\nVérification - Position des maxima après registration:")
registered_landmarks <- apply(registered_curves, 2, function(y) {
  idx <- which.max(y[6:(length(y)-5)]) + 5
  t[idx]
})
print(registered_landmarks)
print(paste("Écart-type avant:", round(sd(landmarks), 4)))
print(paste("Écart-type après:", round(sd(registered_landmarks), 4)))