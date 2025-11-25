library(tidyverse)

data_environments <- as_tibble(read.xlsx(xlsxFile = "data/climate_data/corlouer_supp_data.xlsx", fillMergedCells = TRUE, colNames = TRUE,startRow = 4, rows = c(4:51), na.strings = "-"))



data_climate_aggregated <-  as_tibble(read.xlsx(xlsxFile = "data/climate_data/corlouer_supp_data.xlsx", fillMergedCells = TRUE, colNames = TRUE,startRow = 4, rows = c(4:83), sheet = "Supp Data 3"))



ggplot(data_environments) +
  geom_point(aes(x = `Year(a)`, y = `Seed.yield(g)`, colour = `NNI.(h)`))

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
    y = "Ratio (pente / sd)",
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













