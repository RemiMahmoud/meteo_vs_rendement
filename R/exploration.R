library(tidyverse)
library(openxlsx)


data_environments <- as_tibble(
  read.xlsx(xlsxFile = "data/climate_data/corlouer_supp_data.xlsx",
            fillMergedCells = TRUE, colNames = TRUE,startRow = 4,
            rows = c(4:51), na.strings = "-")) |> 
  rename("Year" = "Year(a)",
         "Location" = "Location(b)", 
         "Environment.name" = "Environment.name(c)",
         "Nitrogen" = "Nitrogen(d)",
         "Population" = "Population(e)",
         "MET" = "MET(f)",
         "Seed.yield" = "Seed.yield(g)",
         "NNI." = "NNI.(h)",         
         "Envirotype." = "Envirotype.(i)")



data_climate_aggregated <-  as_tibble(
  read.xlsx(xlsxFile = "data/climate_data/corlouer_supp_data.xlsx",
            fillMergedCells = TRUE, colNames = TRUE,startRow = 4,
            rows = c(4:83), sheet = "Supp Data 3"))



ggplot(data_environments) +
  geom_point(aes(x = `Year`, y = `Seed.yield`, colour = `NNI.`))

data <-  read.xlsx(xlsxFile = "data/climate_data/SuppData.xlsx", fillMergedCells = TRUE, colNames = TRUE)




plot_rend_par_var <- function(var){
  ggplot(data) + aes(x = {{var}}, y = `Seed.yield(g)`) + 
    geom_point() + geom_smooth(method = "lm")
}
plot_rend_par_var(TMIN_F)


plot_rend_par_var <- function(var){
  ggplot(data) + aes(x = .data[[var]], y = .data[["Seed.yield(g)"]]) + 
    geom_point() + geom_smooth(method = "lm")
}

vars <- as.list(colnames(data)[12:90])

lapply(vars, plot_rend_par_var)

######################  Ratio pente/sd ###############################


ratio_var <- function(var){
  
  x <- data[[var]]
  y <- data[["Seed.yield(g)"]]
  
  mod <- lm(y ~ x)
  slope <- coef(mod)[2]
  
  ratio <- slope * sd(x, na.rm = TRUE) / sd(y, na.rm = TRUE)
  
  return(ratio)
}
vars <- colnames(data)[12:90]



ratios <- sapply(vars, ratio_var)
df_ratios <- tibble(
  variable = vars,
  ratio = ratios,
  abs_ratio = abs(ratios)
)

ggplot(df_ratios, aes(x = reorder(variable, ratio), 
                      y = abs_ratio, 
                      fill = abs_ratio)) +
  geom_col() +
  scale_fill_gradient(low = "lightblue", high = "red") +
  coord_flip() +
  labs(
    x = "Variable climatique",
    y = "Ratio (pente *sd(x)/ sd(y))",
    fill = "|ratio|",
    title = "Importance relative des variables climatiques"
  ) +
  theme_minimal()

######################  R2  ###############################

calc_R2 <- function(var) {
  x <- data[[var]]
  y <- data[["Seed.yield(g)"]]
  
  mod <- lm(y ~ x)
  
  summary(mod)$r.squared
}


vars <- colnames(data)[12:90]

R2_values <- sapply(vars, calc_R2)

df_R2 <- tibble(
  variable = vars,
  R2 = R2_values
)

ggplot(df_R2, aes(x = reorder(variable, R2), y = R2, fill = R2)) +
  geom_col() +
  scale_fill_gradient(low = "lightblue", high = "red") +
  coord_flip() +
  labs(
    x = "Variable climatique",
    y = "R²",
    fill = "R²",
    title = "Variance expliquée par chaque variable climatique"
  ) +
  theme_minimal()



################################################################################


data_climate_aggregated$ind.name <- sub("_.*", "", data_climate_aggregated$Indicator.name)

data_climate_aggregated_pivot <- data_climate_aggregated |> 
  pivot_longer(cols = 5:51, names_to = "Environment.name")

table(data_climate_aggregated_pivot$Period, data_climate_aggregated_pivot$ind.name)  
# présente pour toutes les saisons (sauf NA)
# LSR, SSR, TMAX, TMIN, TMN, WSC, HT

data_climate_aggregated_pivot_propre <- data_climate_aggregated_pivot |> 
  extract(
    col = Environment.name, 
    into = c("lieux", "annee", "nitro"),
    regex = "^([A-Za-z]{2,3})([0-9]{2})_(N\\+|N-)$",
    remove = FALSE
  ) |> 
  left_join(data_environments)

var <- "TMIN"

data_climate_aggregated_pivot_propre |> 
  filter(ind.name == var) |> 
  filter(nitro == "N-") |> 
  mutate(Period = factor(Period, levels = unique(Period))) |> 
  ggplot() +
  aes(x = Period, y = value,
      color = Seed.yield, alpha = 0.8,
      shape = nitro,
      group = Environment.name
      ) +        
  geom_line() +                      
  geom_jitter(aes(size = Seed.yield), width = 0.1, height = 0, alpha = 0.8) +
  scale_color_gradient(low = "lightgreen", high = "darkgreen") +
  labs(y = var)




###################################   ACP

library(FactoMineR)
data_acp <- data_climate_aggregated[, -c(1, 2, 3, 4)] |> 
  column_to_rownames("ind.name") |> 
  t() |> 
  as.data.frame()

PCA(data_acp)

PCA(data_acp[ , c("VERN", "TMAX_FLO", "HT_P600", "SSR_P600")])



mod.pca = PCA(data_acp)


library(factoextra)
fviz_pca_var(mod.pca, select.var = list(cos2 = 5))


library(Factoshiny)

Factoshiny(data_acp)
# minin clustering des var et des pacerelles





#####   PLS

require(leaps)             # For variable selection
require(glmnet)            # For penalized regression procedures
require(fields)            # For image.plot
require(pls)               # For cvsegments, plsr, ...
require(MASS)              # For LDA
require(viridis)           # For nice colors

data_climate_aggregated_small <- data_climate_aggregated[, -c(1:4)]

str(data_climate_aggregated_small)

n = nrow(data_climate_aggregated_small)
cvpred = rep(0,n)

segs = cvsegments(n,10)
cvpred = rep(0,n)

lmp.pls = plsr(ind.name~.,data=data_climate_aggregated_small,ncomp=1,scale=TRUE)

for (k in 1:10) {
  cvmod = plsr(ind.name~.,data=data_climate_aggregated_small[-segs[[k]],],ncomp=30,scale=TRUE,validation="CV",segments=10)
  bestncomp = selectNcomp(lmp.pls,method="onesigma")
  cvpred[segs[[k]]] = predict(cvmod,newdata=data_climate_aggregated_small[segs[[k]],])[,,bestncomp]
  print(k)
}

PRESS.pls = sum((data_climate_aggregated_small$ind.name-cvpred)^2) ; PRESS.pls


