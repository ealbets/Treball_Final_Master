# Paquet xgboost
if (!require("glmnet")) install.packages("glmnet")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("progress")) install.packages("progress")

library(glmnet)
library(dplyr)
library(ggplot2)
library(progress)

path_images <- "data/output/images/ElasticNet/"

##-- MACHINE-LEARNING SUPERVISAT AMB REGRESSIÓ: Algorisme d'Elastic Net --##

# Variable objectiu--> grau histològic (categòrica ordinal)
# Variables independents --> conjunt de 10.600 columnes referents a gens (valors continus normalitzats entre 0 i 1)

cat("-------ELASTIC NET ALGORITHM----------\n")

# Funció genèrica per a trobar els hiperpàrmatres òptims (aquells que maximitzen ) emprant la validació creuada
# i la graella d'opcions de param_grid
hyperparams_en_search <- function(X, y, alpha, lambdes) {
  cat("Cerca d'hiperparàmetres òptims amb Elastic Net: \n")
  
  # Inicialitzar les variables per guardar el millor AUC i els seus hiperparàmetres
  best_auc <- -Inf
  best_params <- NULL
  
  
  # Bucle sobre valors alpha
  for (i in 1:length(alpha)) {
    
    alpha_val <- alpha[i]
    
    # Imprimir els hiperparàmetres actuals
    cat("\n Avaluant hiperparàmetres en Elastic Net ", y ," : ", 
        "alpha =" , alpha_val, " lambda seq: ", lambdes, "\n")
    
    
    # Validació creuada emprant .cv, la variable objectiu es binominal
    cv_model <- cv.glmnet(
      X, dataset[[y]],
      alpha = alpha_val,
      lambda = lambdes,
      family = "binomial",
      type.measure = "auc",
      nfolds = 9 # plecs en la validació creuada, com més petit menys temps de càlcul però menys precisió
    )
    
    # Obtenir el millor AUC per aquests hiperparàmetres
    current_auc <- max(cv_model$cvm)
    # la millor lambda
    best_lambda_for_alpha <- cv_model$lambda.min
    
    # Si el AUC actual és millor que el millor fins ara, en el quedem. Així capturem el que maximitza best_auc
    if (current_auc > best_auc) {
      best_auc <- current_auc
      best_alpha <- alpha_val
      best_lambda <- best_lambda_for_alpha
    }
    
  }
  
  # Retornar els millors hiperparàmetres
  return(list(best_alpha = best_alpha, best_lambda = best_lambda, best_auc = best_auc))
  
}

# Definim les possibilitats d'hiperparàmetres per a cercar els òptims
# seqüència de valors alpha i valors lambda
alpha <- seq(0.1, 0.95, by = 0.05)
lambdes <- 10^seq(-4, 2, length.out = 50)  # 50 valors entre 0.0001 i 100


# Cridar a la funció per a cada un dels graus histològics
best_hyperparams_elasticnet_grade_1 <- hyperparams_en_search(X,"is_grade_1",alpha,lambdes)
best_hyperparams_elasticnet_grade_2 <- hyperparams_en_search(X,"is_grade_2",alpha,lambdes)
best_hyperparams_elasticnet_grade_3 <- hyperparams_en_search(X,"is_grade_3",alpha,lambdes)

# Mostra els millors hiperparàmetres per a cada grau histològic
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 1, són: [alpha] --> ", best_hyperparams_elasticnet_grade_1$best_alpha,
    " [lambda] -->", best_hyperparams_elasticnet_grade_1$best_lambda ,"\n")
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 2, són: [alpha] --> ", best_hyperparams_elasticnet_grade_2$best_alpha,
    " [lambda] -->", best_hyperparams_elasticnet_grade_2$best_lambda ,"\n")
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 3, són: [alpha] --> ", best_hyperparams_elasticnet_grade_3$best_alpha,
    " [lambda] -->", best_hyperparams_elasticnet_grade_3$best_lambda ,"\n")


# Funció que realitza tot el procés de:
# - Divisió del conjunt de dades en entrenament (2/3 parts) i test (1/3)
# - Entrenament del conjunt de dades de train amb els hiperparàmetres òptims trobats
# - Prediccions i assignacions dels resultats a partir del conjunt de dades de test
# - Obtenció de la matriu de confusió com a indicador de rendiment
# - Càlcul de les importàncies concretes per a cada una de es variables objectiu indicadores de grau histològic
# El retorn es el conjunt de resultats dels coeficients d'importància, les mètriques ROC i la matriu de confusió 
elasticnet_process_metrics_importances_by_grade <- function(X, histological_grade, best_alpha, best_lambda) {
  
  cat("Cerca d'importàncies ",histological_grade, " \n")
  
  #Divisió del conjunt de dades en subconjunts de train (entrenament, 70% per conveni 2/3) i test (proves, 30% 1/3)
  set.seed(123)  # Reproducibilitat
  
  cat("--Generació de subconjunts de dades de train i test en" , histological_grade ,"--\n")
  
  datasets_results_elasticnet <- train_test_datasets_generation(X,dataset,histological_grade)
  
  # Capturem resultats
  X_train <- as.matrix(datasets_results_elasticnet$X_train)
  X_test <- as.matrix(datasets_results_elasticnet$X_test)
  y_train <- datasets_results_elasticnet$y_train
  y_test <- datasets_results_elasticnet$y_test
  
  cat("Entrenament per a ",histological_grade, " \n")
  # Entrenament del model final amb els millors hiperparàmetres
  # Volem saber la importància de cada gen respecte cada tipus de variable objectiu, per tant, no caldrà subdidir entre train i test
  # ja que no es volen realitzar prediccions de classificació. Tot el conjunt 'X' serà d'entrenament
  en_model <- glmnet(
    X_train, y_train,
    alpha = best_alpha,
    lambda = best_lambda,
    family = "binomial"
  )
  
  
  # Prediccions sobre el conjunt de test
  cat("Prediccions sobre test en ",histological_grade, " \n")
  y_pred_prob <- predict(en_model, s = best_lambda, newx = X_test, type = "response")
  y_pred_prob <- as.vector(y_pred_prob)  # Assegurem que sigui un vector numèric
  # Convertim les prediccions en resultats binaris fent ús del 0.5 com a llindar
  # Resultats superiors a 0.5 s'etiqueten com a 1 i inferiors, com a 0
  y_pred_labels <- ifelse(y_pred_prob > 0.5, 1, 0)
  
  # Generació matriu de confusió
  confusion_matrix <- generate_confusion_matrix(y_pred_labels, y_test, classes = c(0, 1))
  
  # Càlcul de les mètriques ROC
  roc_curve <- roc(y_test, y_pred_prob)
  
  cat("Cerca d'importàncies ",histological_grade, " \n")
  
  # Obtenir i retornar la importància de les característiques però amb tot el conjunt de dades sencer, sense la divisió entre
  # entrenament i test ja que ens interessa és tenir-ho el més genèric possible
  
  # Entrenament del model amb els hiperparàmetres òptims i el conjunt de dades sencer
  en_model_final <- glmnet(
    X, dataset[[histological_grade]],
    alpha = best_alpha,
    lambda = best_lambda,
    family = "binomial"
  )
  
  # Obtenir les importàncies de les característiques a partir dels coeficients
  coeficients <- coef(en_model_final, s = best_lambda)
  feature_importance <- data.frame(
    Feature = rownames(as.matrix(coeficients)),
    Coef = as.vector(as.matrix(coeficients))
  )
  
  # Ordenem els coeficients per valor absolut (importància).
  # Els valors negatius representen importàncies positivies de "protecció" i els positius de "risc" de patiment.
  feature_importance <- feature_importance %>%
    filter(Feature != "(Intercept)") %>%  # Excloure l'intercept.
    mutate(AbsCoef = abs(Coef)) %>%       # Calcular el valor absolut dels coeficients.
    arrange(desc(AbsCoef))               # Ordenar pel valor absolut de major a menor.
  
  # retorn de les importàncies, la matriu de confusió i l'àrea sota la corba
  return(list(
    importance = feature_importance,
    confusion_matrix = confusion_matrix,
    roc = roc_curve
  ))
  
}

# Establim el nombre de mostres que es volen mostrar en el rànking dels millors i els pitjors
top_num = 10

# A partir de la funció definida anteriorment, obtenim les importàncies de les característiques (gens) 
# per a cada classe (grau histològic: 1, 2 i 3)
# importàncies GRAU HISTOLÒGIC 1
results_elasticnet_grade_1 <- elasticnet_process_metrics_importances_by_grade(X,"is_grade_1", best_hyperparams_elasticnet_grade_1$best_alpha,best_hyperparams_elasticnet_grade_1$best_lambda)
# importàncies GRAU HISTOLÒGIC 2
results_elasticnet_grade_2 <- elasticnet_process_metrics_importances_by_grade(X,"is_grade_2", best_hyperparams_elasticnet_grade_2$best_alpha,best_hyperparams_elasticnet_grade_2$best_lambda)
# importàncies GRAU HISTOLÒGIC 3
results_elasticnet_grade_3 <- elasticnet_process_metrics_importances_by_grade(X,"is_grade_3", best_hyperparams_elasticnet_grade_3$best_alpha,best_hyperparams_elasticnet_grade_3$best_lambda)

# Captura importàncies
importancies_elasticnet_grade1 <- results_elasticnet_grade_1$importance
importancies_elasticnet_grade2 <- results_elasticnet_grade_2$importance
importancies_elasticnet_grade3 <- results_elasticnet_grade_3$importance

# Captura matriu confusió
matriu_confusio_elasticnet_grade1 <- results_elasticnet_grade_1$confusion_matrix
matriu_confusio_elasticnet_grade2 <- results_elasticnet_grade_2$confusion_matrix
matriu_confusio_elasticnet_grade3 <- results_elasticnet_grade_3$confusion_matrix

# Captura AUC
auc_elasticnet_grade1 <- results_elasticnet_grade_1$roc$auc
auc_elasticnet_grade2 <- results_elasticnet_grade_2$roc$auc
auc_elasticnet_grade3 <- results_elasticnet_grade_3$roc$auc

# Mostra del TOP 10 en el rànking de millors resultats d'importància de gens o conjunts de gens per a cada grau histològic
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 1
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 1 --> Gens/s més importants:\n")
print(head(importancies_elasticnet_grade1, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 2 --> Gens/s més importants:\n")
print(head(importancies_elasticnet_grade2, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 3 --> Gens/s més importants:\n")
print(head(importancies_elasticnet_grade3, top_num))

# Mostra les matrius de confusió en cada grau histològic
cat("MATRIU DE CONFUSIÓ GRAU HISTOLÒGIC 1 : \n")
print(matriu_confusio_elasticnet_grade1)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("MATRIU DE CONFUSIÓ GRAU HISTOLÒGIC 2 :\n")
print(matriu_confusio_elasticnet_grade2)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("MATRIU DE CONFUSIÓ GRAU HISTOLÒGIC 3 :\n")
print(matriu_confusio_elasticnet_grade3)


# Totals i percentatges respecte al total global
# Casos reals positius en cada un dels graus histològics
total_grade1 <- sum(matriu_confusio_elasticnet_grade1[, 2])
total_grade2 <- sum(matriu_confusio_elasticnet_grade2[, 2])
total_grade3 <- sum(matriu_confusio_elasticnet_grade3[, 2])

# Total global de casos de la clase "1"
total_class1_train <- total_grade1 + total_grade2 + total_grade3

cat("Total elements de Grau Histològic 1 en train: ", total_grade1, "(", round((total_grade1 / total_class1_train) * 100, 2), "% del total)\n")
cat("Total elements de Grau Histològic 2 en train: ", total_grade2, "(", round((total_grade2 / total_class1_train) * 100, 2), "% del total)\n")
cat("Total elements de Grau Histològic 3 en train: ", total_grade3, "(", round((total_grade3 / total_class1_train) * 100, 2), "% del total)\n")

# Càlcul de mètriques de rendiment a partir de la matriu de confusió i la fòrmula genèrica 'calcula_metrics'
metrics_cm_grade1 <- calcula_metrics(matriu_confusio_elasticnet_grade1)
metrics_cm_grade2 <- calcula_metrics(matriu_confusio_elasticnet_grade2)
metrics_cm_grade3 <- calcula_metrics(matriu_confusio_elasticnet_grade3)

# Mostra resultats
cat("Mètriques de Matriu de Confusió per a Grau Histològic 1:\n")
cat(paste(names(metrics_cm_grade1), "=", unlist(metrics_cm_grade1), collapse = ", "), "\n")
cat("Mètriques de Matriu de Confusió per a Grau Histològic 2:\n")
cat(paste(names(metrics_cm_grade2), "=", unlist(metrics_cm_grade2), collapse = ", "), "\n")
cat("Mètriques de Matriu de Confusió per a Grau Histològic 3:\n")
cat(paste(names(metrics_cm_grade3), "=", unlist(metrics_cm_grade3), collapse = ", "), "\n")

# Mostra les mètriques AUC per a cada grau histològic
cat("AUC en GRAU HISTOLÒGIC 1 : \n")
print(auc_elasticnet_grade1)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("AUC en GRAU HISTOLÒGIC 2 :\n")
print(auc_elasticnet_grade2)
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("AUC en GRAU HISTOLÒGIC 3 :\n")
print(auc_elasticnet_grade3)

# Visualització de resultants emprant gràfiques verticals amb 'ggplot'

# Afegim una columna per categoritzar els coeficients com a positius o negatius
importance_elasticnet_grade_1 <- importancies_elasticnet_grade1 %>%
  mutate(Signe = ifelse(Coef > 0, "Positiu (Risc)", "Negatiu (Protecció)"))

# Gràfic amb diferenciació de colors per signe del coeficient.
# blau pels coeficients positius (sentit negatiu, ja que denota risc de patiment)
# vermell pels coeficients negatius (sentit positiu, ja que denota protecció)

#definició variables
xlabel <- "Gens"
ylabel <- "Coeficient"
color1 <- "blue"
color2 <- "red"
title1 <- paste0("Elastic Net - Top ", top_num, " Importància Genètica (Grau Histològic 1)")
subtitle <- "Coeficients positius (risc) i negatius (protecció)"

plot1 <- create_horitzontal_barchart_with_sign_plot(importance_elasticnet_grade_1, "Feature", "Coef", "Signe", color1, color2, title1, subtitle, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_ElasticNet.png"), plot = plot1, width = 15, height = 10)


# Afegim una columna per categoritzar els coeficients com a positius o negatius
importance_elasticnet_grade_2 <- importancies_elasticnet_grade2 %>%
  mutate(Signe = ifelse(Coef > 0, "Positiu (Risc)", "Negatiu (Protecció)"))

title2 <- paste0("Elastic Net - Top ", top_num, " Importància Genètica (Grau Histològic 2)")

plot2 <- create_horitzontal_barchart_with_sign_plot(importance_elasticnet_grade_2, "Feature", "Coef", "Signe", color1, color2, title2, subtitle, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_ElasticNet.png"), plot = plot2, width = 15, height = 10)


# Afegim una columna per categoritzar els coeficients com a positius o negatius
importance_elasticnet_grade_3 <- importancies_elasticnet_grade3 %>%
  mutate(Signe = ifelse(Coef > 0, "Positiu (Risc)", "Negatiu (Protecció)"))

title3 <- paste0("Elastic Net - Top ", top_num, " Importància Genètica (Grau Histològic 3)")

plot3 <- create_horitzontal_barchart_with_sign_plot(importance_elasticnet_grade_3, "Feature", "Coef", "Signe", color1, color2, title3, subtitle, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_ElasticNet.png"), plot = plot3, width = 15, height = 10)


# Generació datasets amb especifictats i sensibilitats a partir de roc()
roc_data_grade1 <- data.frame(
  Specificity = results_elasticnet_grade_1$roc$specificities,
  Sensitivity = results_elasticnet_grade_1$roc$sensitivities
)

roc_data_grade2 <- data.frame(
  Specificity = results_elasticnet_grade_2$roc$specificities,
  Sensitivity = results_elasticnet_grade_2$roc$sensitivities
)

roc_data_grade3 <- data.frame(
  Specificity = results_elasticnet_grade_3$roc$specificities,
  Sensitivity = results_elasticnet_grade_3$roc$sensitivities
)

# mostra diagrames corbes ROC
plot_roc1 <- plot_roc(roc_data_grade1, "Corba ROC - Grau Histològic 1 ")
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_ElasticNet_ROC.png"), plot = plot_roc1, width = 15, height = 10)
plot_roc2 <- plot_roc(roc_data_grade2, "Corba ROC - Grau Histològic 2 ")
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_ElasticNet_ROC.png"), plot = plot_roc2, width = 15, height = 10)
plot_roc3 <- plot_roc(roc_data_grade3, "Corba ROC - Grau Histològic 3 ")
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_ElasticNet_ROC.png"), plot = plot_roc3, width = 15, height = 10)