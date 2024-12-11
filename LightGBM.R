# Paquet xgboost
if (!require("lightgbm")) install.packages("lightgbm")
if (!require("pROC")) install.packages("pROC")

library(xgboost)
library(dplyr)
library(ggplot2)
library(progress)
library(reshape2)
library(pROC)

path_images <- "data/output/images/LightGBM/"

##-- MACHINE-LEARNING SUPERVISAT AMB REGRESSIÓ: Algorisme de LightGBM --##

cat("-------LIGHTGBM ALGORITHM----------\n")

# Funció genèrica per a trobar els hiperpàrmatres òptims (aquells que maximitzen AUC ) emprant la validació creuada
# i la graella d'opcions de param_grid
hyperparams_search_lgbm <- function(X, y, param_grid) {
  cat("Cerca d'hiperparàmetres òptims: \n")
  
  # Inicialitzar les variables per guardar el millor AUC i els seus hiperparàmetres
  best_auc <- -Inf
  best_params <- NULL
  best_nrounds <- NULL
  
  # Barra de progrés per a tenir una noció del temps
  pb <- progress_bar$new(
    format = " Avaluant hiperparàmetres en LightGBM :current/:total (:percent) \n",
    total = nrow(param_grid),
    clear = FALSE,
    width = 60
  )
  
  # Bucle sobre la graella
  for (i in 1:nrow(param_grid)) {
    
    # Actualitzar la barra de progrés
    pb$tick()
    
    params <- list(
      objective = param_grid$objective[i],
      num_leaves = param_grid$num_leaves[i],
      learning_rate = param_grid$learning_rate[i],
      bagging_fraction = param_grid$bagging_fraction[i],
      feature_fraction = param_grid$feature_fraction[i],
      max_depth = param_grid$max_depth[i]
    )
    
    # Imprimir els hiperparàmetres actuals
    cat("\n Avaluant hiperparàmetres LightGBM en ", y ," : ", 
        paste(names(params), unlist(params), sep = "=", collapse = ", "), "\n")
    
    
    # Número màxim de rondes (es fixa alt) i early stopping
    max_nrounds <- 1000
    early_stopping_rounds <- 20
    
    
    # Validació creuada emprant .cv 
    cv_results <- lgb.cv(
      params = params,
      data = lgb.Dataset(data = X, label = dataset[[y]]),
      nrounds = max_nrounds, # iteracions
      nfold = 3, # particions o plecs del conjunt de dades, 1/3 per provar i 2/3 per entrenar
      eval = "auc", # mètrica de rendiment
      verbose = -1,
      early_stopping_rounds = early_stopping_rounds
    )
    
    
    # Obtenir el millor AUC per aquests hiperparàmetres
    current_auc <- max(unlist(cv_results$record_evals$valid$auc$eval))
    
    # Si el AUC actual és millor que el millor fins ara, en el quedem. Així capturem el que maximitza best_auc
    if (current_auc > best_auc) {
      best_auc <- current_auc
      best_params <- params
      best_nrounds <- cv_results$best_iter # afegim la iteració òptima trobada
    }
    
  }
  
  # Retornar els millors hiperparàmetres
  return(list(best_params = best_params, best_nrounds = best_nrounds, best_auc = best_auc))
  
}

# Degut a les limitacions i les restriccions de l'algorisme LightGBM amb els noms de les columnes que actuen com a característiques
# assignem un mapa de variables i canviem els noms per X0... Xn (on 'n' es el nombre total de característiques)
column_map_lgbm <- setNames(paste0("X", seq_along(colnames(X)) - 1), colnames(X))
# Realitzem la substitució
colnames(X) <- column_map_lgbm[colnames(X)]

# Definim el grid de possibilitats d'hiperparàmetres per a cercar els òptims
param_grid_lgbm <- expand.grid(
  num_leaves = c(15, 31, 63), # nombre de fulles
  learning_rate = c(0.01, 0.1, 0.3), # taxa d'aprenentatge
  max_depth = c(5, 10, -1), # profunditat màxima, -1 significa sense límits en profunditat màxima en LightGBM
  bagging_fraction = c(0.6, 0.8, 1), # fracció de dades per iteració
  feature_fraction = c(0.5, 0.7, 1), # fracció de característiques per iteració
  objective = "binary" # problema de classificació binària
)

# Cridar a la funció per a cada un dels graus histològics
best_hyperparams_lgbm_grade_1 <- hyperparams_search_lgbm(X,"is_grade_1",param_grid_lgbm)
best_hyperparams_lgbm_grade_2 <- hyperparams_search_lgbm(X,"is_grade_2",param_grid_lgbm)
best_hyperparams_lgbm_grade_3 <- hyperparams_search_lgbm(X,"is_grade_3",param_grid_lgbm)

# Mostra els millors hiperparàmetres per a cada grau histològic
cat("Els millors hiperparàmetres trobats a LightGBM per al grau histològic 1 són:",
    paste(names(best_hyperparams_lgbm_grade_1$best_params), "=", unlist(best_hyperparams_lgbm_grade_1$best_params), collapse = ", "),
    ", amb nrounds =", best_hyperparams_lgbm_grade_1$best_nrounds, "\n")

cat("Els millors hiperparàmetres trobats a LightGBM per al grau histològic 2 són:",
    paste(names(best_hyperparams_lgbm_grade_2$best_params), "=", unlist(best_hyperparams_lgbm_grade_2$best_params), collapse = ", "),
    ", amb nrounds =", best_hyperparams_lgbm_grade_2$best_nrounds, "\n")

cat("Els millors hiperparàmetres trobats a LightGBM per al grau histològic 3 són:",
    paste(names(best_hyperparams_lgbm_grade_3$best_params), "=", unlist(best_hyperparams_lgbm_grade_3$best_params), collapse = ", "),
    ", amb nrounds =", best_hyperparams_lgbm_grade_3$best_nrounds, "\n")


# Funció que realitza:
# - Divisió del conjunt de dades en entrenament (2/3 parts) i test (1/3)
# - Entrenament del conjunt de dades de train amb els hiperparàmetres òptims trobats
# - Prediccions dels resultats a partir del conjunt de dades de test
# - Obtenció de la matriu de confusió com a indicador de rendiment
# - Càlcul de les importàncies concretes per a cada una de es variables objectiu indicadores de grau histològic
# El retorn es el conjunt de resultats d'importàncies, de mètriques auc i la matriu de confusió 
lightgbm_train_predict_importances_by_grade <- function(X, histological_grade, best_params, best_nrounds) {
  
  #Divisió del conjunt de dades en subconjunts de train (entrenament, 70% per conveni 2/3) i test (proves, 30% 1/3)
  set.seed(123)  # Reproducibilitat
  
  cat("--Generació de subconjunts de dades de train i test en" , histological_grade ,"--\n")
  
  datasets_results <- train_test_datasets_generation(X,dataset,histological_grade)
  
  # Capturem resultats
  X_train <- as.matrix(datasets_results$X_train)
  X_test <- as.matrix(datasets_results$X_test)
  y_train <- datasets_results$y_train
  y_test <- datasets_results$y_test
  
  cat("Entrenament per a ",histological_grade, " \n")
  
  # L'algorisme requereix un Dataset especial per LightGBM
  # Els generem per al conjunt de train i el de test
  dtrain <- lgb.Dataset(data = X_train, label = y_train)
  dtest <- lgb.Dataset(data = X_test, label = y_test)
  
  # Entrenament del model amb els hiperparàmetres òptims
  model <- lgb.train(params = best_params, data = dtrain, nrounds = best_nrounds, verbose = -1)
  
  cat("Prediccions sobre test en ",histological_grade, " \n")
  # realització de les prediccions sobre el conjunt de dades de test
  y_pred <- predict(model, newdata = X_test)
  # Convertim les prediccions en resultats binaris fent ús del 0.5 com a llindar
  # Resultats superiors a 0.5 s'etiqueten com a 1 i inferiors, com a 0
  y_pred_labels <- ifelse(y_pred > 0.5, 1, 0)
  
  # Càlcul de la matriu de confusió segons els resultats obtinguts en la predicció i la realitat
  # Inicialització de la matriu de confusió amb zeros
  confusion_matrix <- matrix(0, nrow = 2, ncol = 2, 
                        dimnames = list(Predicted = c(0, 1), Actual = c(0, 1)))
  # Càlcul
  confusion_matrix <- table(Predicted = y_pred_labels, Actual = y_test)
  # obtenció de la mètrica de rendiment de l'àrea sota la corba
  roc <- roc(y_test, y_pred)
  
  cat("Cerca d'importàncies ",histological_grade, " \n")
  
  # Obtenir i retornar la importància de les característiques però amb tot el conjunt de dades sencer, sense la divisió entre
  # entrenament i test ja que ens interessa és tenir-ho el més genèric possible
  
  # Entrenament del model amb els hiperparàmetres òptims i el conjunt de dades sencer
  model_final <- lgb.train(
    params = best_params,
    data = lgb.Dataset(data = X, label = dataset[[histological_grade]]),
    nrounds = best_nrounds,
    verbose = -1
  )
  
  # Obtenir i retornar la importància de les característiques 
  importance <- lgb.importance(model = model_final)
  
  # retorn de les importàncies, la matriu de confusió i l'àrea sota la corba
  return(list(
    importance = importance,
    confusion_matrix = confusion_matrix,
    roc = roc
  ))
  
}

# A partir de la funció definida anteriorment, obtenim les importàncies de les característiques (gens) 
# per a cada classe (grau histològic: 1, 2 i 3). Ho executarem 20 vegades per assgurar-nos l'exactitud del càlcul
# i el resultat final serà la mitjana aritmètica de l'acumulació de resultats de les importàncies de cada iteració

# iteracions
#num_times <- 20

# acumulats
#acc_importance_lightgbm_g1 <- NULL
#acc_importance_lightgbm_g2 <- NULL
#acc_importance_lightgbm_g3 <- NULL

# Iterar y acumular resultats
#for (i in 1:num_times) {
  
  #cat("---Iteració: ", i, "---\n")
  # A partir de la funció definida anteriorment, obtenim les importàncies de les característiques (gens) 
  # per a cada classe (grau histològic: 1, 2 i 3)
  # importàncies GRAU HISTOLÒGIC 1
  results_lgbm_grade_1 <- lightgbm_train_predict_importances_by_grade(X, "is_grade_1", best_hyperparams_lgbm_grade_1$best_params, best_hyperparams_lgbm_grade_1$best_nrounds)
  # importàncies GRAU HISTOLÒGIC 2
  results_lgbm_grade_2 <- lightgbm_train_predict_importances_by_grade(X, "is_grade_2", best_hyperparams_lgbm_grade_2$best_params, best_hyperparams_lgbm_grade_2$best_nrounds)
  # importàncies GRAU HISTOLÒGIC 3
  results_lgbm_grade_3 <- lightgbm_train_predict_importances_by_grade(X, "is_grade_3", best_hyperparams_lgbm_grade_3$best_params, best_hyperparams_lgbm_grade_3$best_nrounds)
  


# Captura importàncies
importancies_lightgbm_grade1 <- results_lgbm_grade_1$importance
importancies_lightgbm_grade2 <- results_lgbm_grade_2$importance
importancies_lightgbm_grade3 <- results_lgbm_grade_3$importance
  
# Captura matriu confusió
matriu_confusio_lightgbm_grade1 <- results_lgbm_grade_1$confusion_matrix
matriu_confusio_lightgbm_grade2 <- results_lgbm_grade_2$confusion_matrix
matriu_confusio_lightgbm_grade3 <- results_lgbm_grade_3$confusion_matrix
  
# Captura AUC
auc_lightgbm_grade1 <- results_lgbm_grade_1$roc$auc
auc_lightgbm_grade2 <- results_lgbm_grade_2$roc$auc
auc_lightgbm_grade3 <- results_lgbm_grade_3$roc$auc
  
  
# Establim el nombre de mostres que es volen mostrar en el rànking dels millors
top_num = 10


# Invertir el mapeig inicial del gens / característiques: de noms generats a noms originals
reversed_map_lgbm <- setNames(names(column_map_lgbm), column_map_lgbm)


# Reemplaçar els noms a la columna 'Feature' amb els noms originals
importancies_lightgbm_grade1$Feature <- reversed_map_lgbm[importancies_lightgbm_grade1$Feature]
importancies_lightgbm_grade2$Feature <- reversed_map_lgbm[importancies_lightgbm_grade2$Feature]
importancies_lightgbm_grade3$Feature <- reversed_map_lgbm[importancies_lightgbm_grade3$Feature]

# Resultats finals a partir de la mitjana aritmètica dels acumulats en les iteracions anteriors.
# Excloem la columna 'Features' en l'operació
#importance_lgbm_grade_1 <- importance_lgbm_grade_1 %>% mutate(across(where(is.numeric), ~ . / num_times))
#importance_lgbm_grade_2 <- importance_lgbm_grade_2 %>% mutate(across(where(is.numeric), ~ . / num_times))
#importance_lgbm_grade_3 <- importance_lgbm_grade_3 %>% mutate(across(where(is.numeric), ~ . / num_times))

# Mostra del TOP 15 en el rànking de millors resultats d'importància de gens o conjunts de gens per a cada grau histològic
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 1
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 1 --> Gens/s més importants:\n")
print(head(importancies_lightgbm_grade1, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 2 --> Gens/s més importants:\n")
print(head(importancies_lightgbm_grade2, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 3 --> Gens/s més importants:\n")
print(head(importancies_lightgbm_grade3, top_num))

# Mostra les matrius de confusió en cada grau histològic
cat("MATRIU DE CONFUSIÓ GRAU HISTOLÒGIC 1 : \n")
print(matriu_confusio_lightgbm_grade1)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("MATRIU DE CONFUSIÓ GRAU HISTOLÒGIC 2 :\n")
print(matriu_confusio_lightgbm_grade2)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("MATRIU DE CONFUSIÓ GRAU HISTOLÒGIC 3 :\n")
print(matriu_confusio_lightgbm_grade3)


# Totals i percentatges respecte al total global
# Casos reals positius en cada un dels graus histològics
total_grade1 <- sum(matriu_confusio_lightgbm_grade1[, 2])
total_grade2 <- sum(matriu_confusio_lightgbm_grade2[, 2])
total_grade3 <- sum(matriu_confusio_lightgbm_grade3[, 2])

# Total global de casos de la clase "1"
total_class1_train <- total_grade1 + total_grade2 + total_grade3

cat("Total elements de Grau Histològic 1 en train: ", total_grade1, "(", round((total_grade1 / total_class1_train) * 100, 2), "% del total)\n")
cat("Total elements de Grau Histològic 2 en train: ", total_grade2, "(", round((total_grade2 / total_class1_train) * 100, 2), "% del total)\n")
cat("Total elements de Grau Histològic 3 en train: ", total_grade3, "(", round((total_grade3 / total_class1_train) * 100, 2), "% del total)\n")

# Càlcul de mètriques de rendiment a partir de la matriu de confusió i la fòrmula genèrica 'calcula_metrics'
metrics_cm_grade1 <- calcula_metrics(matriu_confusio_lightgbm_grade1)
metrics_cm_grade2 <- calcula_metrics(matriu_confusio_lightgbm_grade2)
metrics_cm_grade3 <- calcula_metrics(matriu_confusio_lightgbm_grade3)

# Mostra resultats
cat("Mètriques de Matriu de Confusió per a Grau Histològic 1:\n")
cat(paste(names(metrics_cm_grade1), "=", unlist(metrics_cm_grade1), collapse = ", "), "\n")
cat("Mètriques de Matriu de Confusió per a Grau Histològic 2:\n")
cat(paste(names(metrics_cm_grade2), "=", unlist(metrics_cm_grade2), collapse = ", "), "\n")
cat("Mètriques de Matriu de Confusió per a Grau Histològic 3:\n")
cat(paste(names(metrics_cm_grade3), "=", unlist(metrics_cm_grade3), collapse = ", "), "\n")

# Mostra les mètriques AUC per a cada grau histològic
cat("AUC en GRAU HISTOLÒGIC 1 : \n")
print(auc_lightgbm_grade1)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("AUC en GRAU HISTOLÒGIC 2 :\n")
print(auc_lightgbm_grade2)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("AUC en GRAU HISTOLÒGIC 3 :\n")
print(auc_lightgbm_grade3)


# Visualització de resultants emprant gràfiques verticals amb 'ggplot'
feature_lbgm_importance1 <- importancies_lightgbm_grade1 %>% arrange(desc(Gain))
feature_lbgm_importance2 <- importancies_lightgbm_grade2 %>% arrange(desc(Gain))
feature_lbgm_importance3 <- importancies_lightgbm_grade3 %>% arrange(desc(Gain))


# Crida a la funció que ploteja els diagrames de barres hortizontals per a les característiques importants del dataset
xlabel <- "Gens"
ylabel <- "Gain"
title1 <- paste0("LightGBM - ", top_num, " Importància Genètica (Grau Histològic 1) segons Gain")
color1 <- "blue"
plot1 <- create_horitzontal_barchart_plot(feature_lbgm_importance1, "Feature", "Gain", color1, title1, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_LightGBM.png"), plot = plot1, width = 15, height = 10)

title2 <- paste0("LightGBM - ", top_num, " Importància Genètica (Grau Histològic 2) segons Gain")
color2 <- "green"
plot2 <- create_horitzontal_barchart_plot(feature_lbgm_importance2, "Feature", "Gain", color2, title2, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_LightGBM.png"), plot = plot2, width = 15, height = 10)

title3 <- paste0("LightGBM - ", top_num, " Importància Genètica (Grau Histològic 3) segons Gain")
color3 <- "red"
plot3 <- create_horitzontal_barchart_plot(feature_lbgm_importance3, "Feature", "Gain", color3, title3, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_LightGBM.png"), plot = plot3, width = 15, height = 10)


# Crida a la funció que ploteja els diagrames de punts per a les característiques importants del dataset
# Visualització de resultats emprant diagrames de punts amb 'ggplot'
xlabel <- "Importància segons Gain"
ylabel <- "Gens"
title1 <- paste0("LightGBM - ", top_num, " Importància Genètica (Grau Histològic 1) segons Gain")
color1 <- "blue"
plot1 <- create_point_chart_plot(feature_lbgm_importance1, "Feature", "Gain", color1, title1, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_LightGBM_Points.png"), plot = plot1, width = 15, height = 10)

title2 <- paste0("LightGBM - ", top_num, " Importància Genètica (Grau Histològic 2) segons Gain")
color2 <- "green"
plot2 <- create_point_chart_plot(feature_lbgm_importance2, "Feature", "Gain", color2, title2, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_LightGBM_Points.png"), plot = plot2, width = 15, height = 10)

title3 <- paste0("LightGBM - ", top_num, " Importància Genètica (Grau Histològic 3) segons Gain")
color3 <- "red"
plot3 <- create_point_chart_plot(feature_lbgm_importance3, "Feature", "Gain", color3, title3, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_LightGBM_Points.png"), plot = plot3, width = 15, height = 10)


# Generar las corbes ROC

# Generació datasets amb especifictats i sensibilitats a partir de roc()
roc_data_grade1 <- data.frame(
  Specificity = results_lgbm_grade_1$roc$specificities,
  Sensitivity = results_lgbm_grade_1$roc$sensitivities
)

roc_data_grade2 <- data.frame(
  Specificity = results_lgbm_grade_2$roc$specificities,
  Sensitivity = results_lgbm_grade_2$roc$sensitivities
)

roc_data_grade3 <- data.frame(
  Specificity = results_lgbm_grade_3$roc$specificities,
  Sensitivity = results_lgbm_grade_3$roc$sensitivities
)

# mostra diagrames corbes ROC
plot_roc1 <- plot_roc(roc_data_grade1, "Corba ROC - Grau Histològic 1 ")
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_LightGBM_ROC.png"), plot = plot_roc1, width = 15, height = 10)
plot_roc2 <- plot_roc(roc_data_grade2, "Corba ROC - Grau Histològic 2 ")
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_LightGBM_ROC.png"), plot = plot_roc2, width = 15, height = 10)
plot_roc3 <- plot_roc(roc_data_grade3, "Corba ROC - Grau Histològic 3 ")
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_LightGBM_ROC.png"), plot = plot_roc3, width = 15, height = 10)
