# Paquet xgboost
if (!require("xgboost")) install.packages("xgboost")
if (!require("progress")) install.packages("progress")
if (!require("reshape2")) install.packages("reshape2")

library(xgboost)
library(dplyr)
library(progress)
library(reshape2)

path_images <- "data/output/images/XGBoost/"

##-- MACHINE-LEARNING SUPERVISAT AMB REGRESSIÓ: Algorisme de XGBoost --##

# Variable objectiu--> grau histològic (categòrica ordinal)
# Variables independents --> conjunt de 10.600 columnes referents a gens (valors continus normalitzats entre 0 i 1)

cat("-------XGBOOST ALGORITHM----------\n")

# Crear variables binàries per a cada grau histològic
dataset$is_grade_1 <- ifelse(dataset$HISTOLOGICAL_GRADE == 1, 1, 0)
dataset$is_grade_2 <- ifelse(dataset$HISTOLOGICAL_GRADE == 2, 1, 0)
dataset$is_grade_3 <- ifelse(dataset$HISTOLOGICAL_GRADE == 3, 1, 0)


# Funció genèrica per a trobar els hiperpàrmatres òptims (aquells que maximitzen ) emprant la validació creuada
# i la graella d'opcions de param_grid
# Paràmetres:
# - X :conjunt de variables independents)
# - y :variable objectiu binària: is_grade_n)
# - param_grid: matriu inicial de possibles hiperparàmetres amb possibles valors
hyperparams_xgboost_search <- function(X, y, param_grid) {
  cat("Cerca d'hiperparàmetres òptims: \n")
  
  # Inicialitzar les variables per guardar el millor AUC i els seus hiperparàmetres
  best_auc <- -Inf
  best_params <- NULL
  best_nrounds <- NULL
  
  # Barra de progrés per a tenir una noció del temps
  pb <- progress_bar$new(
    format = "Avaluant hiperparàmetres :current/:total (:percent) \n",
    total = nrow(param_grid),
    clear = FALSE,
    width = 60
  )
  
  # Número màxim de rondes (es fixa alt) i early stopping
  max_nrounds <- 1000
  early_stopping_rounds <- 15
  
  # Bucle sobre la graella
  for (i in 1:nrow(param_grid)) {
    
    # Actualitzar la barra de progrés
    pb$tick()
    
    params <- list(
      objective = "binary:logistic",
      max_depth = param_grid$max_depth[i],
      eta = param_grid$eta[i],
      subsample = param_grid$subsample[i],
      colsample_bytree = param_grid$colsample_bytree[i]
    )
    
    # Imprimir els hiperparàmetres actuals
    cat("\n Avaluant hiperparàmetres en XGBoost ", y ," : ", 
        paste(names(params), unlist(params), sep = "=", collapse = ", "), "\n")
    
    # Validació creuada emprant xgb.cv amb early stopping
    cv_results <- xgb.cv(
      params = params,
      data = xgb.DMatrix(X, label = dataset[[y]]),
      nrounds = max_nrounds,
      nfold = 3,
      metrics = "auc", # corba error
      verbose = 0,
      early_stopping_rounds = early_stopping_rounds
    )
    
    # Obtenir el millor AUC per aquests hiperparàmetres
    current_auc <- max(cv_results$evaluation_log$test_auc_mean)
    
    # Si el AUC actual és millor que el millor fins ara, en el quedem. Així capturem el que maximitza best_auc
    if (current_auc > best_auc) {
      best_auc <- current_auc
      best_params <- params
      best_nrounds <- cv_results$best_iteration # afegim la iteració òptima trobada
    }
  }
  
  # Retornar els millors hiperparàmetres
  return(list(best_params = best_params, best_nrounds = best_nrounds, best_auc = best_auc))
}

# Capturem les variables independents que són les corresponents als gens o conjunts de gens. Aquestes van de la columna 7
# a la N (10647) però n'hem d'excloure les 3 últimes variables temporals afegides en el pas anterior d'indicadors de grau específic.
# Ho convertim en matriu per a que sigui compatible amb l'algorisme.
X <- as.matrix(dataset[, 7:(ncol(dataset) - 3)])

# Definim el grid de possibilitats d'hiperparàmetres per a cercar els òptims
param_grid_xgboost <- expand.grid(
  max_depth = c(3, 6, 9),
  eta = c(0.01, 0.1, 0.3), # taxa aprenentatge
  subsample = c(0.6, 0.8, 1), # mostres de l'arbre en cada iteració
  colsample_bytree = c(0.5, 0.7, 1), # columnes per entrenar arbres en cada iteració
  objective = "binary:logistic"  # Classificació binària (is_grade)
)

# Cridar a la funció per a cada un dels graus histològics
best_hyperparams_xgboost_grade_1 <- hyperparams_xgboost_search(X,"is_grade_1",param_grid_xgboost)
best_hyperparams_xgboost_grade_2 <- hyperparams_xgboost_search(X,"is_grade_2",param_grid_xgboost)
best_hyperparams_xgboost_grade_3 <- hyperparams_xgboost_search(X,"is_grade_3",param_grid_xgboost)

# Mostra els millors hiperparàmetres per a cada grau histològic
cat("Els millors hiperparàmetres trobats per al grau histològic 1 són:",
    paste(names(best_hyperparams_xgboost_grade_1$best_params), "=", unlist(best_hyperparams_xgboost_grade_1$best_params), collapse = ", "),
    ", amb nrounds =", best_hyperparams_xgboost_grade_1$best_nrounds, "\n")

cat("Els millors hiperparàmetres trobats per al grau histològic 2 són:",
    paste(names(best_hyperparams_xgboost_grade_2$best_params), "=", unlist(best_hyperparams_xgboost_grade_2$best_params), collapse = ", "),
    ", amb nrounds =", best_hyperparams_xgboost_grade_2$best_nrounds, "\n")

cat("Els millors hiperparàmetres trobats per al grau histològic 3 són:",
    paste(names(best_hyperparams_xgboost_grade_3$best_params), "=", unlist(best_hyperparams_xgboost_grade_3$best_params), collapse = ", "),
    ", amb nrounds =", best_hyperparams_xgboost_grade_3$best_nrounds, "\n")

# Funció genèrica per obtenir la rellevància dels diferents gens o conjunts de gens a partir de cada grau histològic diferent i
# el conjunt d'hiperparàmetres òptims trobats amb la corss validation grid Search
# La variable objectiu ve marcada pel paràmetre indicador de grau i el retorn es el conjunt de resultats d'importàncies 
# per al grau concret.
xgboost_importance_by_grade <- function(X, histological_grade, best_params, best_nrounds) {
  
  cat("Cerca d'importàncies ",histological_grade, " \n")
  
  # L'algorisme requereix d'una matriu en format DMatrix
  # Volem saber la importància de cada gen respecte cada tipus de variable objectiu, per tant, no caldrà subdidir entre train i test
  # ja que no es volen realitzar prediccions de classificació. Tot el conjunt 'X' serà d'entrenament.
  dtrain <- xgb.DMatrix(data = X, label = dataset[[histological_grade]])
  
  # Entrenem model amb la matriu i els paràmetres
  model_final <- xgb.train(params = best_params, data = dtrain, nrounds = best_nrounds)
  
  # Obtenir i retornar la importància de les característiques 
  importance <- xgb.importance(model = model_final)
  
  return(importance)
  
}



# A partir de la funció definida anteriorment, obtenim les importàncies de les característiques (gens) 
# per a cada classe (grau histològic: 1, 2 i 3). Ho executarem 20 vegades per assgurar-nos l'exactitud del càlcul
# i el resultat final serà la mitjana aritmètica de l'acumulació de resultats de les importàncies de cada iteració

# iteracions
num_times <- 20

# acumulats
acc_importance_xgboost_g1 <- NULL
acc_importance_xgboost_g2 <- NULL
acc_importance_xgboost_g3 <- NULL

# Iterar y acumular resultats
for (i in 1:num_times) {

  cat("---Iteració: ", i, "---\n")
  # importàncies GRAU HISTOLÒGIC 1
  importance_xgboost_grade_1 <- xgboost_importance_by_grade(X,"is_grade_1", best_hyperparams_xgboost_grade_1$best_params, best_hyperparams_xgboost_grade_1$best_nrounds)
  # importàncies GRAU HISTOLÒGIC 2
  importance_xgboost_grade_2 <- xgboost_importance_by_grade(X,"is_grade_2", best_hyperparams_xgboost_grade_2$best_params, best_hyperparams_xgboost_grade_2$best_nrounds)
  # importàncies GRAU HISTOLÒGIC 3
  importance_xgboost_grade_3 <- xgboost_importance_by_grade(X,"is_grade_3", best_hyperparams_xgboost_grade_3$best_params, best_hyperparams_xgboost_grade_3$best_nrounds)

  # Acumulació
  if (is.null(acc_importance_xgboost_g1) && is.null(acc_importance_xgboost_g2) && is.null(acc_importance_xgboost_g3)) {
    # assignació la primera vegada
    accum_importance_grade_1 <- importance_xgboost_grade_1
    accum_importance_grade_2 <- importance_xgboost_grade_2
    accum_importance_grade_3 <- importance_xgboost_grade_3
  } else {
    # suma resultats quan no es la primera iteració
    accum_importance_grade_1 <- accum_importance_grade_1 + importance_xgboost_grade_1
    accum_importance_grade_2 <- accum_importance_grade_2 + importance_xgboost_grade_2
    accum_importance_grade_3 <- accum_importance_grade_3 + importance_xgboost_grade_3
  }

}

# Resultats finals a partir de la mitjana aritmètica dels acumulats en les iteracions anteriors.
# Excloem la columna 'Features' en l'operació
importance_xgboost_grade_1 <- importance_xgboost_grade_1 %>% mutate(across(where(is.numeric), ~ . / num_times))
importance_xgboost_grade_2 <- importance_xgboost_grade_2 %>% mutate(across(where(is.numeric), ~ . / num_times))
importance_xgboost_grade_3 <- importance_xgboost_grade_3 %>% mutate(across(where(is.numeric), ~ . / num_times))


# Establim el nombre de mostres que es volen mostrar en el rànking dels millors
top_num = 10

# Mostra del TOP 10 en el rànking de millors resultats d'importància de gens o conjunts de gens per a cada grau histològic
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 1
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 1 --> Gens/s més importants:\n")
print(head(importance_xgboost_grade_1, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 2 --> Gens/s més importants:\n")
print(head(importance_xgboost_grade_2, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 3 --> Gens/s més importants:\n")
print(head(importance_xgboost_grade_3, top_num))


# Visualització de resultats emprant gràfiques verticals amb 'ggplot'
# Resultats segons 'Gain'
feature_importance1 <- importance_xgboost_grade_1 %>% arrange(desc(Gain))
feature_importance2 <- importance_xgboost_grade_2 %>% arrange(desc(Gain))
feature_importance3 <- importance_xgboost_grade_3 %>% arrange(desc(Gain))

# Crida a la funció genèrica que ploteja els diagrames de barres hortizontals per a les característiques importants del dataset
xlabel <- "Gens"
ylabel <- "Gain"
title1 <- paste0("XGBoost - ", top_num, " Importància Genètica (Grau Histològic 1) segons Gain")
color1 <- "blue"
plot1 <- create_horitzontal_barchart_plot(feature_importance1, "Feature", "Gain", color1, title1, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_XGBoost.png"), plot = plot1, width = 15, height = 10)

title2 <- paste0("XGBoost - ", top_num, " Importància Genètica (Grau Histològic 2) segons Gain")
color2 <- "green"
plot2 <- create_horitzontal_barchart_plot(feature_importance2, "Feature", "Gain", color2, title2, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_XGBoost.png"), plot = plot2, width = 15, height = 10)
  
title3 <- paste0("XGBoost - ", top_num, " Importància Genètica (Grau Histològic 3) segons Gain")
color3 <- "red"
plot3 <- create_horitzontal_barchart_plot(feature_importance3, "Feature", "Gain", color3, title3, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_XGBoost.png"), plot = plot3, width = 15, height = 10)


# Crida a la funció que ploteja els diagrames de punts per a les característiques importants del dataset
# Visualització de resultats emprant diagrames de punts amb 'ggplot'
xlabel <- "Importància segons Gain"
ylabel <- "Gens"
title1 <- paste0("XGBoost - ", top_num, " Importància Genètica (Grau Histològic 1) segons Gain")
color1 <- "blue"
plot1 <- create_point_chart_plot(feature_importance1, "Feature", "Gain", color1, title1, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_XGBoost_Points.png"), plot = plot1, width = 15, height = 10)

title2 <- paste0("XGBoost - ", top_num, " Importància Genètica (Grau Histològic 2) segons Gain")
color2 <- "green"
plot2 <- create_point_chart_plot(feature_importance2, "Feature", "Gain", color2, title2, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_XGBoost_Points.png"), plot = plot2, width = 15, height = 10)

title3 <- paste0("XGBoost - ", top_num, " Importància Genètica (Grau Histològic 3) segons Gain")
color3 <- "red"
plot3 <- create_point_chart_plot(feature_importance3, "Feature", "Gain", color3, title3, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_XGBoost_Points.png"), plot = plot3, width = 15, height = 10)
