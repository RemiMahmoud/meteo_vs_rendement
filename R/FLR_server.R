
# Tidyverse and plotting functions
library(dplyr)
library(stringr)
library(forcats)
library(readxl)
library(ggplot2)
library(readr)
library(tidyr)
library(lubridate)
library(refund)
library(purrr)
library(pROC)
library(MLmetrics)
library(PRROC)
library(scoring)
library(ResourceSelection)



running_platform = sessionInfo()$running


if(grepl("Ubuntu",running_platform)){
  
  # go to root for ubuntu arborescence
  setwd("..")
  
  # input_dir = "simulated_datasets"
  # output_dir = "results"
  # compression_file = "brotli"
  
  
} else if(grepl("Windows", running_platform))  {
  # 
  # input_dir = "R/simulations/simulated_datasets"
  # output_dir = "models/results"
  # compression_file = "zstd"
}





extract_scalar_effect <-  function(fit, scalar_cov, se_mult = 2) {
  
  stopifnot(inherits(fit, "pfr"))
  
  
  summary_coeff <- summary(fit)$p.table
  
  df <- data.frame(coef = rownames(summary_coeff),
                   beta_hat = summary_coeff[,"Estimate"],
                   lower = summary_coeff[,"Estimate"] - se_mult * summary_coeff[, "Std. Error"],
                   upper = summary_coeff[,"Estimate"] + se_mult * summary_coeff[, "Std. Error"]
  )
  
  
  if(length(scalar_cov) > 0) {df = df%>% 
    mutate(variable = str_extract(coef, paste0(c("Intercept", scalar_cov),  collapse ="|"))) } else {
      df = df %>% 
  mutate(variable = "Intercept") 
    }
  
  rownames(df) = NULL
  
  return(df)
}


extract_functional_effect <- function(fit, which_component = 1, se_mult = 2) {
  
  stopifnot(inherits(fit, "pfr"))
  
  
  # Récupération des données de plot via mgcv
  pdf(NULL)
  plt <- mgcv::plot.gam(fit,seWithMean = TRUE, rug = FALSE, pages = 0)
  dev.off()
  
  
  
  if (!is.list(plt) || length(plt) < which_component) {
    stop("Impossible d'extraire les données pour le terme sélectionné.")
  }
  
  # Récupérer les valeurs
  x_vals <- plt[[which_component]]$x
  fit_vals <- plt[[which_component]]$fit[,1]
  se_vals <- plt[[which_component]]$se
  
  # Créer un data.frame pour ggplot2
  df <- data.frame(
    t = x_vals,
    beta_hat = fit_vals,
    se = se_vals,
    lower = fit_vals - se_mult * se_vals,
    upper = fit_vals + se_mult * se_vals
  )
  
  return(df)
}



get_intervals <- function(v, base_diff_time_interval = NULL, tol  = 1e-8) {
  
  
  # Estime le pas si non fourni
  if (is.null(base_diff_time_interval)) {
    base_diff_time_interval <- median(diff(v))
  }
  
  diffs <- diff(v)
  is_gap <- abs(diffs - base_diff_time_interval) > tol
  
  split_indices <- c(0, which(is_gap), length(v))
  
  # split_indices <- c(0, which(diff(v) != base_diff_time_interval), length(v))
  intervals <- data.frame(start = v[split_indices[-length(split_indices)] + 1],
                          end = v[split_indices[-1]])
  return(intervals)
}


fit_pfr_multi_fun_covariates <- function(var_list,
                                         data_climate,
                                         ID_train,
                                         ID_test,
                                         response_var = "NIV_sup_thresh",
                                         scalar_covariates = NULL,
                                         fun_params = NULL,
                                         k_spline = 30,
                                         family = "binomial") {
  
  # helper : retourne les paramètres pour chaque variable (sinon valeurs par défaut)
  get_params <- function(var) {
    if (!is.null(fun_params) && var %in% fun_params$var) {
      row <- fun_params %>% filter(var == !!var)
      list(k = row$k, bs = row$bs, m = row$m)
    } else {
      list(k = k_spline, bs = "ps", m = NULL)
    }
  }
  
  # helper pour pivoting
  get_pivoted_data <- function(var, ID_data) {
    data_climate %>%
      inner_join(ID_data, by = "ID_ECH") %>%
      distinct(ID_ECH, .data[[var]], diff_flo_dd, .data[[response_var]], concentration, !!!syms(scalar_covariates)) %>%
      pivot_wider(values_from = all_of(var), names_from = "diff_flo_dd", names_prefix = "T")
  }
  
  # initialisation
  data_train_list <- list()
  data_test_list <- list()
  matrix_train_list <- list()
  matrix_test_list <- list()
  
  for (var in var_list) {
    data_train_var <- get_pivoted_data(var, ID_train)
    data_test_var  <- get_pivoted_data(var, ID_test)
    
    mat_train <- data_train_var %>% select(-c(ID_ECH, any_of(c(contains("concentration"), contains("sup_thresh"))), all_of(scalar_covariates))) %>% as.matrix()
    mat_test  <- data_test_var  %>% select(-c(ID_ECH, any_of(c(contains("concentration"), contains("sup_thresh"))), all_of(scalar_covariates))) %>% as.matrix()
    
    matrix_train_list[[var]] <- mat_train
    matrix_test_list[[var]] <- mat_test
    
    if (length(data_train_list) == 0) {
      data_train_list[["meta"]] <- data_train_var %>% select(ID_ECH, all_of(response_var), all_of(scalar_covariates))
      data_test_list[["meta"]]  <- data_test_var  %>% select(ID_ECH, all_of(response_var), all_of(scalar_covariates))
    }
  }
  
  # assemblage final
  data_train <- data_train_list[["meta"]]
  data_test  <- data_test_list[["meta"]]
  
  for (var in var_list) {
    data_train[[var]] <- I(matrix_train_list[[var]])
    data_test[[var]]  <- I(matrix_test_list[[var]])
  }
  
  # construire la formule
  lf_terms <- purrr::map_chr(var_list, function(var) {
    p <- get_params(var)
    m_str <- if (!is.null(p$m)) paste0(", m = ", p$m) else ""
    paste0("lf(", var, ", k = ", p$k, ", bs = '", p$bs, "'", m_str, ")")
  })
  
  rhs <- c(lf_terms, scalar_covariates) %>% paste(collapse = " + ")
  full_formula <- as.formula(paste(response_var, "~", rhs))
  
  # fit
  fit <- pfr(full_formula, data = data_train, family = "binomial")
  
  # prédiction
  data_test <- data_test %>%
    mutate(pred = predict(fit, newdata = data_test, type = "response"))
  
  return(list(
    model = fit,
    data_train = data_train,
    data_test = data_test,
    formula = full_formula
  ))
}





eval_model_metrics <- function(model, data_test, probas = NULL,
                               threshold = "F1", plot = TRUE, response_var = "NIV_sup_thresh", verbose = FALSE) {
  
  # 1. Prédictions si non fournies
  if (is.null(probas)) {
    probas <- predict(model, newdata = data_test, type = "response")
    
  }
  # which_na_proba = which(probas == 1000)
  which_na_proba <- which(is.na(probas))
  if(length(which_na_proba) > 0){
    
    
    if(verbose) {print(paste0("Some missing value(s), dropping observation(s) ", paste0(which_na_proba), collapse = " - "))}
    
    probas = probas[- which_na_proba]
    y_true <- as.numeric(data_test[[response_var]])
    y_true <- y_true[-which_na_proba]
  } else {
    y_true <- as.numeric(data_test[[response_var]])
  }
  
  
  # 2. Détermination du seuil optimal si besoin
  if (is.character(threshold)) {
    roc_obj <- roc(y_true, probas)
    thresholds <- roc_obj$thresholds
    sensitivities <- roc_obj$sensitivities
    specificities <- roc_obj$specificities
    
    if (threshold == "youden") {
      youden_index <- sensitivities + specificities - 1
      threshold <- thresholds[which.max(youden_index)]
    } else if (threshold == "F1") {
      f1_scores <- sapply(thresholds, function(t) {
        y_pred_tmp <- ifelse(probas > t, 1, 0)
        TP <- sum(y_true == 1 & y_pred_tmp == 1)
        FP <- sum(y_true == 0 & y_pred_tmp == 1)
        FN <- sum(y_true == 1 & y_pred_tmp == 0)
        precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
        recall <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
        if (precision + recall == 0) return(0)
        return(2 * precision * recall / (precision + recall))
      })
      threshold <- thresholds[which.max(f1_scores)]
    } else {
      stop("Seuil inconnu : utilisez un nombre ou 'F1' ou 'youden'")
    }
  }
  
  # 3. Prédiction binaire
  y_pred <- ifelse(probas > threshold, 1, 0)
  
  # 4. Courbes PR et ROC
  roc_obj <- roc(y_true, probas)
  auc_roc <- auc(roc_obj)
  
  pr_obj <- pr.curve(scores.class0 = probas,
                     weights.class0 = y_true,
                     curve = TRUE)
  auc_pr <- pr_obj$auc.integral
  
  # 5. Brier Score
  # brier <- BrierScore(probas, y_true)
  brier <- mean((probas - y_true)^2)
  
  # 6. Confusion matrix & métriques
  TP <- sum(y_true == 1 & y_pred == 1)
  TN <- sum(y_true == 0 & y_pred == 0)
  FP <- sum(y_true == 0 & y_pred == 1)
  FN <- sum(y_true == 1 & y_pred == 0)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  f1 <- ifelse((precision + sensitivity) == 0, 0, 2 * precision * sensitivity / (precision + sensitivity))
  
  # accordance between distrib and frequencies
  hoslem <- hoslem.test(y_true, probas, g = 10) 
  
  
  
  plot_calibration <- function(probas, y_true, n_bins = 10) {
    stopifnot(length(probas) == length(y_true))
    
    df <- data.frame(probas = probas, y_true = y_true) %>%
      mutate(bin = ntile(probas, n_bins)) %>% 
      group_by(bin) %>%
      summarise(
        pred_mean = mean(probas),
        obs_rate = mean(y_true),
        n = n(),
        .groups = "drop"
      )
    
    ggplot(df, aes(x = pred_mean, y = obs_rate)) +
      geom_point(size = 2) +
      geom_line() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      labs(
        x = "Probabilité prédite",
        y = "Fréquence observée",
        title = "Courbe de calibration"
      ) +
      theme_minimal()
  }
  # 7. Affichage des courbes
  if (plot) {
    par(mfrow = c(1, 3))
    plot(roc_obj, main = "ROC Curve", col = "blue")
    plot(pr_obj, main = "Precision-Recall Curve", col = "darkgreen")
    # plot(calibration.plot(y = y_true, p = probas), main = "Calibration Curve")
    plot_calibration(probas, y_true)
    par(mfrow = c(1, 1))
  }
  
  # 8. Résumé
  metrics <- list(
    AUC_ROC = as.numeric(auc_roc), # avoid weird auc type 
    AUC_PR  = auc_pr,
    Threshold = threshold,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    Brier = brier,
    hoslem_pval = hoslem$p.value, 
    AIC = AIC(model)
  )
  
  return(metrics)
}



data_itk <- readxl::read_xlsx("data/data_itk_coordinates.xlsx") %>% 
  mutate(all_itk = paste(ID_ESPECE, PRECEDENT2_IDTRAVAIL3, sep = "_")) 
# 
data_climate <-data.table::fread(file= "data/data_all_climate_degree_days_Pierre.csv") %>%
  mutate(ID_ECH = as.factor(ID_ECH)) %>%
  left_join(data_itk %>% distinct(ID_ECH, ID_ESPECE, NUM_ANNEE, ID_ESPECE, ID_TRAVAIL_SOL_3, PRECEDENT2)) 

# data_climate_all <- readRDS("data/data_climate_flr_all.rds")

table_sowing_harvest_date <- data_itk %>% distinct( ID_ECH, NUM_ANNEE, Semis.date, Recolte.date)


table_legal_thresholds <- tibble(toxin = rep(c("ZEA", "DON", "T2_HT2", "NIV","ENNIATINE"),2), legal_threshold = c(100,1500,100, 100,100, 100, 1000, 50, 100,100) , ID_ESPECE = rep(c("BLE DUR", "BLE TENDRE"), each = 5))

table_threshold_toxins <- read_xlsx("data/data_arvalis.xlsx", sheet = "metadata_mycotoxins_fusa") %>% 
  distinct(labo, toxin, threshold) %>% 
  rename(lab_threshold = threshold) %>% 
  filter(!is.na(toxin))  %>% 
  mutate(max_threshold = max(as.numeric(str_replace(lab_threshold, "ND", "0")) ), .by = toxin)  %>% distinct(labo, toxin, max_threshold, lab_threshold)%>% 
  left_join(table_legal_thresholds) %>% 
  mutate(lab_threshold = as.numeric(lab_threshold)) 




data_toxins <- read_xlsx("data/data_arvalis.xlsx") %>% 
  select(ID_ECH,ZEA_valeurs, NIV_valeurs, DON_valeurs)  %>% 
  pivot_longer(contains("valeurs"), names_to = "toxin", values_to = "concentration", names_transform =  function(x) str_remove(x, "_valeurs")) %>% 
  left_join(read_xlsx("data/data_arvalis.xlsx") %>% 
              select(ID_ECH,ZEA_labo, NIV_labo, DON_labo) %>% 
              pivot_longer(contains("labo"), names_to = "toxin", values_to = "labo", names_transform =  function(x) str_remove(x, "_labo"))) %>% 
  left_join(data_itk %>%
              distinct(ID_ECH, ID_ESPECE)) %>% 
  left_join(table_threshold_toxins %>%
              distinct(labo, toxin, max_threshold, lab_threshold, legal_threshold, ID_ESPECE)) %>% # deal with MYCSA ND later
  mutate(above_max_threshold = concentration  >
           max_threshold) %>% 
  mutate(above_lab_threshold = concentration  > 
           lab_threshold) %>% 
  mutate(above_legal_threshold = concentration  > legal_threshold)  %>% 
  mutate(above_red_leg_thresh = concentration  > legal_threshold / if_else(toxin == "NIV", 2,4))%>% 
  filter(!is.na(concentration)) %>% 
  left_join(data_itk %>% 
              select(ID_ECH, ID_ESPECE,  PRECEDENT2, ID_TRAVAIL_SOL_3, PRECEDENT2_IDTRAVAIL3))




# name_toxin = "NIV"
#name_toxin = "DON"
# name_toxin = "FG"


vec_toxins <- c("NIV", "DON", "FG")



fun_vars <- c("P_ETP", "Dm", "Rg", "Amplitude", "Geslin")
scalar_vars <- c("ID_ESPECE", "ID_TRAVAIL_SOL_3", "PRECEDENT2")
cols = c("AUC_ROC", "AUC_PR", "Threshold", "Accuracy", "Sensitivity", 
         "Specificity", "Precision", "F1", "Brier", "hoslem_pval", "AIC")

fun_combinations <- unlist(lapply(1:length(fun_vars), function(k) combn(fun_vars, k, simplify = FALSE)), recursive = FALSE)
scalar_combinations <- unlist(lapply(0:length(scalar_vars), function(k) combn(scalar_vars, k, simplify = FALSE)), recursive = FALSE)

n_fits <- 50


for(i in 1:length(vec_toxins)){
    
  name_toxin = vec_toxins[i]
  
  if(name_toxin == "FG"){
    
    
    
    table_threshold_grami <- tibble(fungi = "FG", lab_threshold = .5, impact_threshold = 3.5)
    
    data_grami <- read_xlsx("data/data_arvalis.xlsx") %>% 
      select(ID_ECH,FG_valeurs)  %>% 
      pivot_longer(contains("valeurs"), names_to = "fungi", values_to = "concentration", names_transform =  function(x) str_remove(x, "_valeurs")) %>% 
      left_join(read_xlsx("data/data_arvalis.xlsx") %>% 
                  select(ID_ECH,FG_labo) %>% 
                  pivot_longer(contains("labo"), names_to = "fungi", values_to = "labo", names_transform =  function(x) str_remove(x, "_labo"))) %>% 
      left_join(data_itk %>%
                  distinct(ID_ECH, ID_ESPECE)) %>% 
      left_join(table_threshold_grami) %>%
      mutate(grami_sup_thresh_impact = concentration  >
               impact_threshold) %>% 
      mutate(grami_sup_thresh_lab = concentration  > 
               lab_threshold) %>% 
      filter(!is.na(concentration)) %>% 
      distinct(ID_ECH, concentration, grami_sup_thresh_impact, grami_sup_thresh_lab)
    
  }
  
  
  data_toxin <- data_toxins %>% 
    filter(toxin == name_toxin) %>% 
    mutate(legal_threshold = legal_threshold / if_else(name_toxin == "NIV", 2,4), toxin_sup_thresh = concentration > legal_threshold) %>% 
    distinct(ID_ECH, concentration, toxin_sup_thresh)
  
  response_var <- "toxin_sup_thresh"
  
  
  if(name_toxin == "FG"){
    data_toxin <- data_grami
    response_var <- "grami_sup_thresh_lab"
   # response_var <- "grami_sup_thresh_impact"
  }
  
  
  
  
  all_IDs =  data_climate %>% 
    left_join(data_toxin) %>% 
    filter(!is.na(concentration))%>% 
    distinct(ID_ECH, concentration, .data[[response_var]] ) 
  
  
  results_all <- list()
  
  comb_index <- 1
  
  k_splines <- c(10,20,30)
  for (k in k_splines){
    for(fun_covs in fun_combinations){
      for(scalar_covs in scalar_combinations){
      
      cat(">> Modèle", comb_index, "/", length(fun_combinations) * length(scalar_combinations), "\n")
      
      cat(">> k =", k, "; toxin ", name_toxin, "\n")
      
      tibble_metrics <- tibble::tibble(
        !!!setNames(rep(list(rep(NA_real_, n_fits)), length(cols)), cols)
      ) %>% 
        mutate(ID_train = vector("list", n_fits),
               ID_test = vector("list", n_fits),
               data_plot = vector("list", n_fits),
               data_scalar_cov = vector("list", n_fits),
               id_rep = 1:n_fits,
               .before = 0L,
               fun_covs = paste(fun_covs, collapse = ","),
               scalar_covs = paste(scalar_covs, collapse = ","))
      
      for(i in 1:n_fits){
        ID_train <- all_IDs %>% slice_sample(prop = .7, by = {{response_var}})
        ID_test <- anti_join(all_IDs, ID_train)
        
        fit_pfr <- fit_pfr_multi_fun_covariates(
          var_list = fun_covs,
          data_climate = data_climate,
          ID_train = ID_train,
          ID_test = ID_test,
          scalar_covariates = scalar_covs,
          response_var = response_var,
          k_spline = k
        )
        
        tibble_metrics$ID_train[[i]] <- ID_train$ID_ECH
        tibble_metrics$ID_test[[i]] <- ID_test$ID_ECH
        tibble_metrics$data_plot[[i]] <- lapply(1:length(fun_covs), function(x) extract_functional_effect(fit_pfr$model, x))
        tibble_metrics$data_scalar_cov[[i]] <- extract_scalar_effect(fit_pfr$model, scalar_covs)
        
        metrics <- eval_model_metrics(fit_pfr$model, fit_pfr$data_test, plot = FALSE, response_var = response_var)
        tibble_metrics[i, cols] <- metrics
      }
      
      results_all[[comb_index]] <- tibble_metrics
      
      cat(">> Modèle", comb_index, "/", length(fun_combinations) * length(scalar_combinations), "fini ! \n")
      comb_index <- comb_index + 1
    }
  }
  
  # Combine en un seul data.frame géant si besoin
  results_df <- bind_rows(results_all, .id = "model_id")
  
  #saveRDS(results_df,paste0("models/results_FLR_", name_toxin, ".rds"))
  saveRDS(results_df,paste0("models/results_FLR_k", k, name_toxin,"_",  response_var, ".rds"))

  }
  
}