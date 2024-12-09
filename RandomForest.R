# Paquet Random Forest
if (!require("randomForest")) install.packages("randomForest")
if (!require("progress")) install.packages("progress")
if (!require("caret")) install.packages("caret")

library(randomForest)
library(caret)
library(progress)

path_images <- "data/output/images/RandomForest/"

##-- MACHINE-LEARNING SUPERVISAT AMB REGRESSIÓ: Algorisme de Random Forest --##

cat("-------RANDOM FOREST ALGORITHM----------\n")

# Variable objectiu--> grau histològic (categòrica ordinal)
# Variables independents --> conjunt de 10.600 columnes referents a gens (valors continus normalitzats entre 0 i 1)

# Crear variables binàries per a cada grau histològic
dataset$is_grade_1 <- ifelse(dataset$HISTOLOGICAL_GRADE == 1, 1, 0)
dataset$is_grade_2 <- ifelse(dataset$HISTOLOGICAL_GRADE == 2, 1, 0)
dataset$is_grade_3 <- ifelse(dataset$HISTOLOGICAL_GRADE == 3, 1, 0)


# Funció per a trobar els hiperpàrmatres òptims (aquells que maximitzen els millors valors del Random Forest) emprant la validació creuada
# i la graella d'opcions de param_grid
hyperparams_search_randomForest <- function(X, y, param_grid) {
  cat("Cerca d'hiperparàmetres òptims: \n")
  
  # Inicialitzar les variables per guardar el millor AUC i els seus hiperparàmetres
  best_auc <- -Inf
  best_params <- NULL
  
  # Barra de progrés per a tenir una noció del temps
  pb <- progress_bar$new(
    format = " Avaluant hiperparàmetres en Random Forest :current/:total (:percent) Temps restant: \n",
    total = nrow(param_grid),
    clear = FALSE,
    width = 60
  )
  
  # Bucle sobre la graella
  for (i in 1:nrow(param_grid)) {
    
    # Actualitzar la barra de progrés
    pb$tick()
    
    # Obtenir els hiperparàmetres actuals
    params <- list(
      mtry = param_grid$mtry[i],
      ntree = param_grid$ntree[i],
      nodesize = param_grid$nodesize[i]
    )
    
    # Imprimir els hiperparàmetres actuals
    cat("\n Avaluant hiperparàmetres Random Forest en ", y ," : ", 
        paste(names(params), unlist(params), sep = "=", collapse = ", "), "\n")
    
    # Crear funció de validació creuada ('cv')
    train_control <- trainControl(method = "cv", number = 3, classProbs = TRUE, summaryFunction = twoClassSummary)
    
    # Entrenar el model Random Forest amb validació creuada. No accepta etiquetes numèriques en la variable objectiu.
    # Per tant, les canviem a 'true', 'false'
    rf_model <- train(
      x = X,
      y = as.factor(factor(dataset[[y]],levels = c("0", "1"),labels = c("true", "false"))),
      method = "rf",
      metric = "ROC",
      trControl = train_control,
      tuneGrid = data.frame(mtry = params$mtry),
      ntree = params$ntree,
      nodesize = params$nodesize
    )
    
    # Obtenir el millor AUC per aquests hiperparàmetres
    current_auc <- max(rf_model$results$ROC)
    
    # Si el AUC actual és millor que el millor fins ara, en el quedem. Així capturem el que maximitza best_auc
    if (current_auc > best_auc) {
      best_auc <- current_auc
      best_params <- params
    }
    
  }
  
  return(list(best_auc = best_auc, best_params = best_params))
  
}

# Capturem les variables independents que són les corresponents als gens o conjunts de gens. Aquestes van de la columna 7
# a la N (10647) però n'hem d'excloure les 3 últimes variables temporals afegides en el pas anterior d'indicadors de grau específic.
# Ho convertim en matriu per a que sigui compatible amb l'algorisme.
X <- as.matrix(dataset[, 7:(ncol(dataset) - 3)])

# Degut a les limitacions i les restriccions de l'algorisme RandomForest amb els noms de les columnes que actuen com a característiques
# assignem un mapa de variables i canviem els noms per X0... Xn (on 'n' es el nombre total de característiques)
column_map <- setNames(paste0("X", seq_along(colnames(X)) - 1), colnames(X))
# Realitzem la substitució
colnames(X) <- column_map[colnames(X)]

# Graella de possibles hiperparàmetres per a trobar-ne els valors òptims
param_grid_rf <- expand.grid(
  mtry = c(100, 500, 1000),     # Número de característiques considerades
  ntree = c(500, 1000, 3000),   # Número d'arbres
  nodesize = c(1, 5, 10)      # Mínim tamany de nodes
)

# Cridar a la funció per a cada un dels graus histològics
best_hyperparams_rf_grade_1 <- hyperparams_search_randomForest(X,"is_grade_1",param_grid_rf)
best_hyperparams_rf_grade_2 <- hyperparams_search_randomForest(X,"is_grade_2",param_grid_rf)
best_hyperparams_rf_grade_3 <- hyperparams_search_randomForest(X,"is_grade_3",param_grid_rf)

# Mostra els millors hiperparàmetres per a cada grau histològic
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 1, són:",paste(names(best_hyperparams_rf_grade_1$best_params),":", unlist(best_hyperparams_rf_grade_1$best_params)),"\n")
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 2, són:",paste(names(best_hyperparams_rf_grade_2$best_params),":", unlist(best_hyperparams_rf_grade_2$best_params)),"\n")
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 3, són:",paste(names(best_hyperparams_rf_grade_3$best_params),":", unlist(best_hyperparams_rf_grade_3$best_params)),"\n")


# Funció genèrica per obtenir la rellevància dels diferents gens o conjunts de gens a partir de cada grau histològic diferent i
# el conjunt d'hiperparàmetres òptims trobats amb la corss validation
# La variable objectiu ve marcada pel paràmetre indicador de grau i el retorn es el conjunt de resultats d'importàncies 
# per al grau concret.
rf_importance_by_grade <- function(X,histological_grade, best_params) {
  
  cat("Cerca d'importàncies ",histological_grade, " \n")
  
  # Entrenament de Random Forest amb els millors paràmetres
  model_rf <- randomForest(
    x = X,
    y = factor(dataset[[histological_grade]]),
    ntree = best_params$ntree,
    mtry = best_params$mtry,
    nodesize = best_params$nodesize,
    importance = TRUE  # activem importàncies
  )
  
  # Obtenir i retornar la importància
  importance <- as.data.frame(importance(model_rf))
  # importància per característiques ordenades per la mesura 'MeanDecreaseGini'
  importance$Feature <- rownames(importance)
  importance <- importance[order(-importance$MeanDecreaseGini), ]
  
  return(importance)
}

# Establim el nombre de mostres que es volen mostrar en el rànking dels millors
top_num = 10

# A partir de la funció definida anteriorment, obtenim les importàncies de les característiques (gens) 
# per a cada classe (grau histològic: 1, 2 i 3)
# importàncies GRAU HISTOLÒGIC 1
importance_rf_grade_1 <- rf_importance_by_grade(X,"is_grade_1", best_hyperparams_rf_grade_1$best_params)
# importàncies GRAU HISTOLÒGIC 2
importance_rf_grade_2 <- rf_importance_by_grade(X,"is_grade_2", best_hyperparams_rf_grade_2$best_params)
# importàncies GRAU HISTOLÒGIC 3
importance_rf_grade_3 <- rf_importance_by_grade(X,"is_grade_3", best_hyperparams_rf_grade_3$best_params)

# Invertir el mapeig inicial del gens / característiques: de noms generats a noms originals
reversed_map <- setNames(names(column_map), column_map)
# Reemplaçar els noms a la columna 'Feature' amb els noms originals
importance_rf_grade_1$Feature <- reversed_map[importance_rf_grade_1$Feature]
importance_rf_grade_2$Feature <- reversed_map[importance_rf_grade_2$Feature]
importance_rf_grade_3$Feature <- reversed_map[importance_rf_grade_3$Feature]

# Mostra del TOP 15 en el rànking de millors resultats d'importància de gens o conjunts de gens per a cada grau histològic
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 1
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 1 --> Gens/s més importants:\n")
print(head(importance_rf_grade_1[, c("Feature", "MeanDecreaseGini", "MeanDecreaseAccuracy")], top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 2 --> Gens/s més importants:\n")
print(head(importance_rf_grade_2[, c("Feature", "MeanDecreaseGini", "MeanDecreaseAccuracy")], top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 3 --> Gens/s més importants:\n")
print(head(importance_rf_grade_3[, c("Feature", "MeanDecreaseGini", "MeanDecreaseAccuracy")], top_num))



# Crida a la funció que ploteja els diagrames de barres hortizontals per a les característiques importants del dataset
xlabel <- "Gens"
ylabel <- "Mean Decrease Gini"
title1 <- paste0("RANDOM FOREST - ", top_num, " Importància Genètica (Grau Histològic 1) segons Mean Decrease Gini")
color1 <- "blue"
plot1 <- create_horitzontal_barchart_plot(importance_rf_grade_1, "Feature", "MeanDecreaseGini", color1, title1, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_RandomForest.png"), plot = plot1, width = 15, height = 10)

title2 <- paste0("RANDOM FOREST - ", top_num, " Importància Genètica (Grau Histològic 2) segons Mean Decrease Gini")
color2 <- "green"
plot2 <- create_horitzontal_barchart_plot(importance_rf_grade_2, "Feature", "MeanDecreaseGini", color2, title2, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_RandomForest.png"), plot = plot2, width = 15, height = 10)

title3 <- paste0("RANDOM FOREST - ", top_num, " Importància Genètica (Grau Histològic 3) segons Mean Decrease Gini")
color3 <- "red"
plot3 <- create_horitzontal_barchart_plot(importance_rf_grade_3, "Feature", "MeanDecreaseGini", color3, title3, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_RandomForest.png"), plot = plot3, width = 15, height = 10)


# Visualització de resultats emprant diagrames de punts amb 'ggplot'
xlabel <- "Importància segons MeanDecreaseGini"
ylabel <- "Gens"
title1 <- paste0("RANDOM FOREST - ", top_num, " Importància Genètica (Grau Histològic 1) segons MeanDecreaseGini")
color1 <- "blue"
plot1 <- create_point_chart_plot(importance_rf_grade_1, "Feature", "MeanDecreaseGini", color1, title1, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_RandomForest_Points.png"), plot = plot1, width = 15, height = 10)

title2 <- paste0("RANDOM FOREST - ", top_num, " Importància Genètica (Grau Histològic 2) segons MeanDecreaseGini")
color2 <- "green"
plot2 <- create_point_chart_plot(importance_rf_grade_2, "Feature", "MeanDecreaseGini", color2, title2, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_RandomForest_Points.png"), plot = plot2, width = 15, height = 10)

title3 <- paste0("RANDOM FOREST - ", top_num, " Importància Genètica (Grau Histològic 3) segons MeanDecreaseGini")
color3 <- "red"
plot3 <- create_point_chart_plot(importance_rf_grade_3, "Feature", "MeanDecreaseGini", color3, title3, top_num, xlabel, ylabel)
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_RandomForest_Points.png"), plot = plot3, width = 15, height = 10)
