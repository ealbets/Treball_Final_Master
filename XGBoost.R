# Paquet xgboost
if (!require("xgboost")) install.packages("xgboost")
if (!require("progress")) install.packages("progress")
if (!require("reshape2")) install.packages("reshape2")
if (!require("pROC")) install.packages("pROC")

library(xgboost)
library(dplyr)
library(progress)
library(reshape2)
library(pROC)

path_images <- "data/output/images/XGBoost/"

##-- MACHINE-LEARNING SUPERVISAT AMB REGRESSIÓ: Algorisme de XGBoost --##

# Variable objectiu--> grau histològic (categòrica ordinal)
# Variables independents --> conjunt de 10.600 columnes referents a gens (valors continus normalitzats entre 0 i 1)

cat("-------XGBOOST ALGORITHM----------\n")


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
  early_stopping_rounds <- 20
  
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

# Funció que realitza:
# - Divisió del conjunt de dades en entrenament (2/3 parts) i test (1/3)
# - Entrenament del conjunt de dades de train amb els hiperparàmetres òptims trobats
# - Prediccions dels resultats a partir del conjunt de dades de test
# - Obtenció de la matriu de confusió com a indicador de rendiment
# - Càlcul de les importàncies concretes per a cada una de es variables objectiu indicadores de grau histològic
# El retorn es el conjunt de resultats d'importàncies, de mètriques auc i la matriu de confusió 
xgboost_train_predict_importances_by_grade <- function(X, histological_grade, best_params, best_nrounds) {
  
  cat("--Generació de subconjunts de dades de train i test en" , histological_grade ,"--\n")
  
  #Divisió del conjunt de dades en subconjunts de train (entrenament, 70% per conveni 2/3) i test (proves, 30% 1/3)
  set.seed(123)  # Reproducibilitat
  
  datasets_results <- train_test_datasets_generation(X,dataset,histological_grade)
  
  # Capturem resultats
  X_train <- as.matrix(datasets_results$X_train)
  X_test <- as.matrix(datasets_results$X_test)
  y_train <- datasets_results$y_train
  y_test <- datasets_results$y_test
  
  
  cat("Entrenament per a ",histological_grade, " \n")
  
  # Creacions de DMatrix per a entrenament i test
  # L'algorisme requereix d'una matriu en format DMatrix
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest <- xgb.DMatrix(data = X_test, label = y_test)
  
  # Entrenament del model amb els hiperparàmetres òptims
  model <- xgb.train(
    params = best_params,
    data = dtrain,
    nrounds = best_nrounds
  )
  
  
  cat("Prediccions sobre test en ",histological_grade, " \n")
  # realització de les prediccions sobre el conjunt de dades de test
  y_pred <- predict(model, newdata = dtest)
  
  # Convertim les prediccions en resultats binaris fent ús del 0.5 com a llindar
  # Resultats superiors a 0.5 s'etiqueten com a 1 i inferiors, com a 0
  y_pred_labels <- ifelse(y_pred > 0.5, 1, 0)
  
  # Càlcul de la matriu de confusió segons els resultats obtinguts en la predicció i la realitat
  confusion_matrix <- table(Predicted = y_pred_labels, Actual = y_test)
  # obtenció de la mètrica de rendiment de l'àrea sota la corba
  roc <- roc(y_test, y_pred)
  
  cat("Cerca d'importàncies ",histological_grade, " \n")

  # Obtenir i retornar la importància de les característiques però amb tot el conjunt de dades sencer, sense la divisió entre
  # entrenament i test ja que ens interessa és tenir-ho el més genèric possible
  
  # Entrenament del model amb els hiperparàmetres òptims i el conjunt de dades sencer
  model_final <- xgb.train(
    params = best_params,
    data = xgb.DMatrix(data = X, label = dataset[[histological_grade]]),
    nrounds = best_nrounds
  )
  
  # obtenció importàncies
  importance <- xgb.importance(model = model_final)
  
  # retorn de les importàncies, la matriu de confusió i l'àrea sota la corba
  return(list(
    importance = importance,
    confusion_matrix = confusion_matrix,
    roc = roc
  ))
  
}


xgboost_importances <- function(X, histological_grade, best_params, best_nrounds) {
  
}


# A partir de la funció definida anteriorment, obtenim les importàncies de les característiques (gens) i les
# mètriques de rendiment per a cada classe (grau histològic: 1, 2 i 3). 

#num_times <- 20

# acumulats
#acc_importance_xgboost_g1 <- NULL
#acc_importance_xgboost_g2 <- NULL
#acc_importance_xgboost_g3 <- NULL

# Iterar y acumular resultats
#for (i in 1:num_times) {

  #cat("---Iteració: ", i, "---\n")
  # importàncies GRAU HISTOLÒGIC 1
  results_xgboost_grade_1 <- xgboost_train_predict_importances_by_grade(X,"is_grade_1", best_hyperparams_xgboost_grade_1$best_params, best_hyperparams_xgboost_grade_1$best_nrounds)
  # importàncies GRAU HISTOLÒGIC 2
  results_xgboost_grade_2 <- xgboost_train_predict_importances_by_grade(X,"is_grade_2", best_hyperparams_xgboost_grade_2$best_params, best_hyperparams_xgboost_grade_2$best_nrounds)
  # importàncies GRAU HISTOLÒGIC 3
  results_xgboost_grade_3 <- xgboost_train_predict_importances_by_grade(X,"is_grade_3", best_hyperparams_xgboost_grade_3$best_params, best_hyperparams_xgboost_grade_3$best_nrounds)

  # Acumulació
  #if (is.null(acc_importance_xgboost_g1) && is.null(acc_importance_xgboost_g2) && is.null(acc_importance_xgboost_g3)) {
    # assignació la primera vegada
    #accum_importance_grade_1 <- importance_xgboost_grade_1
    #accum_importance_grade_2 <- importance_xgboost_grade_2
    #accum_importance_grade_3 <- importance_xgboost_grade_3
  #} else {
    # suma resultats quan no es la primera iteració
    #accum_importance_grade_1 <- accum_importance_grade_1 + importance_xgboost_grade_1
    #accum_importance_grade_2 <- accum_importance_grade_2 + importance_xgboost_grade_2
    #accum_importance_grade_3 <- accum_importance_grade_3 + importance_xgboost_grade_3
  #}

#}

# Captura importàncies
importancies_xgboost_grade1 <- results_xgboost_grade_1$importance
importancies_xgboost_grade2 <- results_xgboost_grade_2$importance
importancies_xgboost_grade3 <- results_xgboost_grade_3$importance

# Captura matriu confusió
matriu_confusio_xgboost_grade1 <- results_xgboost_grade_1$confusion_matrix
matriu_confusio_xgboost_grade2 <- results_xgboost_grade_2$confusion_matrix
matriu_confusio_xgboost_grade3 <- results_xgboost_grade_3$confusion_matrix

# Captura AUC
auc_xgboost_grade1 <- results_xgboost_grade_1$roc$auc
auc_xgboost_grade2 <- results_xgboost_grade_2$roc$auc
auc_xgboost_grade3 <- results_xgboost_grade_3$roc$auc
  
# Resultats finals a partir de la mitjana aritmètica dels acumulats en les iteracions anteriors.
# Excloem la columna 'Features' en l'operació
#importance_xgboost_grade_1 <- importance_xgboost_grade_1 %>% mutate(across(where(is.numeric), ~ . / num_times))
#importance_xgboost_grade_2 <- importance_xgboost_grade_2 %>% mutate(across(where(is.numeric), ~ . / num_times))
#importance_xgboost_grade_3 <- importance_xgboost_grade_3 %>% mutate(across(where(is.numeric), ~ . / num_times))


# Establim el nombre de mostres que es volen mostrar en el rànking dels millors
top_num = 10

# Mostra del TOP 10 en el rànking de millors resultats d'importància de gens o conjunts de gens per a cada grau histològic
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 1
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 1 --> Gens/s més importants:\n")
print(head(importancies_xgboost_grade1, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 2 --> Gens/s més importants:\n")
print(head(importancies_xgboost_grade2, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 3 --> Gens/s més importants:\n")
print(head(importancies_xgboost_grade3, top_num))


# Mostra les matrius de confusió en cada grau histològic
cat("MATRIU DE CONFUSIÓ GRAU HISTOLÒGIC 1 : \n")
print(matriu_confusio_xgboost_grade1)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("MATRIU DE CONFUSIÓ GRAU HISTOLÒGIC 2 :\n")
print(matriu_confusio_xgboost_grade2)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("MATRIU DE CONFUSIÓ GRAU HISTOLÒGIC 3 :\n")
print(matriu_confusio_xgboost_grade3)


# Totals i percentatges respecte al total global
# Casos reals positius en cada un dels graus histològics
total_grade1 <- sum(matriu_confusio_xgboost_grade1[, 2])
total_grade2 <- sum(matriu_confusio_xgboost_grade2[, 2])
total_grade3 <- sum(matriu_confusio_xgboost_grade3[, 2])

# Total global de casos de la clase "1"
total_class1_train <- total_grade1 + total_grade2 + total_grade3

cat("Total elements de Grau Histològic 1 en train: ", total_grade1, "(", round((total_grade1 / total_class1_train) * 100, 2), "% del total)\n")
cat("Total elements de Grau Histològic 2 en train: ", total_grade2, "(", round((total_grade2 / total_class1_train) * 100, 2), "% del total)\n")
cat("Total elements de Grau Histològic 3 en train: ", total_grade3, "(", round((total_grade3 / total_class1_train) * 100, 2), "% del total)\n")

# Càlcul de mètriques de rendiment a partir de la matriu de confusió i la fòrmula genèrica 'calcula_metrics'
metrics_cm_grade1 <- calcula_metrics(matriu_confusio_xgboost_grade1)
metrics_cm_grade2 <- calcula_metrics(matriu_confusio_xgboost_grade2)
metrics_cm_grade3 <- calcula_metrics(matriu_confusio_xgboost_grade3)

# Mostra resultats
cat("Mètriques de Matriu de Confusió per a Grau Histològic 1:\n")
cat(paste(names(metrics_cm_grade1), "=", unlist(metrics_cm_grade1), collapse = ", "), "\n")
cat("Mètriques de Matriu de Confusió per a Grau Histològic 2:\n")
cat(paste(names(metrics_cm_grade2), "=", unlist(metrics_cm_grade2), collapse = ", "), "\n")
cat("Mètriques de Matriu de Confusió per a Grau Histològic 3:\n")
cat(paste(names(metrics_cm_grade3), "=", unlist(metrics_cm_grade3), collapse = ", "), "\n")

# Mostra les mètriques AUC per a cada grau histològic
cat("AUC en GRAU HISTOLÒGIC 1 : \n")
print(auc_xgboost_grade1)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("AUC en GRAU HISTOLÒGIC 2 :\n")
print(auc_xgboost_grade2)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("AUC en GRAU HISTOLÒGIC 3 :\n")
print(auc_xgboost_grade3)


# Visualització de resultats emprant gràfiques verticals amb 'ggplot'
# Resultats segons 'Gain'
feature_importance1 <- importancies_xgboost_grade1 %>% arrange(desc(Gain))
feature_importance2 <- importancies_xgboost_grade2 %>% arrange(desc(Gain))
feature_importance3 <- importancies_xgboost_grade3 %>% arrange(desc(Gain))

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



# Generar las corbes ROC

# Generació datasets amb especifictats i sensibilitats a partir de roc()
roc_data_grade1 <- data.frame(
  Specificity = results_xgboost_grade_1$roc$specificities,
  Sensitivity = results_xgboost_grade_1$roc$sensitivities
)

roc_data_grade2 <- data.frame(
  Specificity = results_xgboost_grade_2$roc$specificities,
  Sensitivity = results_xgboost_grade_2$roc$sensitivities
)

roc_data_grade3 <- data.frame(
  Specificity = results_xgboost_grade_3$roc$specificities,
  Sensitivity = results_xgboost_grade_3$roc$sensitivities
)

# mostra diagrames corbes ROC
plot_roc1 <- plot_roc(roc_data_grade1, "Corba ROC - Grau Histològic 1 ")
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_XGBoost_ROC.png"), plot = plot_roc1, width = 15, height = 10)
plot_roc2 <- plot_roc(roc_data_grade2, "Corba ROC - Grau Histològic 2 ")
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_XGBoost_ROC.png"), plot = plot_roc2, width = 15, height = 10)
plot_roc3 <- plot_roc(roc_data_grade3, "Corba ROC - Grau Histològic 3 ")
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_XGBoost_ROC.png"), plot = plot_roc3, width = 15, height = 10)
