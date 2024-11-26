# Paquet xgboost
if (!require("xgboost")) install.packages("xgboost")
if (!require("progress")) install.packages("progress")
if (!require("reshape2")) install.packages("reshape2")

library(xgboost)
library(dplyr)
library(ggplot2)
library(progress)
library(reshape2)

path_images <- "images/XGBoost/"

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
hyperparams_xgboost_search <- function(X, y, param_grid) {
  cat("Cerca d'hiperparàmetres òptims: \n")
  
  # Inicialitzar les variables per guardar el millor AUC i els seus hiperparàmetres
  best_auc <- -Inf
  best_params <- NULL
  
  # Barra de progrés per a tenir una noció del temps
  pb <- progress_bar$new(
    format = "Avaluant hiperparàmetres :current/:total (:percent) Temps restant: \n",
    total = nrow(param_grid),
    clear = FALSE,
    width = 60
  )
  
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
    
  
    # Validació creuada emprant .cv 
    cv_results <- xgb.cv(
      params = params,
      data = xgb.DMatrix(X, label = dataset[[y]]),
      nrounds = 100, # iteracions
      nfold = 3,
      metrics = "auc", # corba error
      verbose = 0
    )

    # Obtenir el millor AUC per aquests hiperparàmetres
    current_auc <- max(cv_results$evaluation_log$test_auc_mean)
    
    # Si el AUC actual és millor que el millor fins ara, en el quedem. Així capturem el que maximitza best_auc
    if (current_auc > best_auc) {
      best_auc <- current_auc
      best_params <- params
    }
    
  }
  
  # Retornar els millors hiperparàmetres
  return(list(best_params = best_params, best_auc = best_auc))
  
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
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 1, són:",paste(names(best_hyperparams_xgboost_grade_1$best_params),":", unlist(best_hyperparams_xgboost_grade_1$best_params)),"\n")
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 2, són:",paste(names(best_hyperparams_xgboost_grade_2$best_params),":", unlist(best_hyperparams_xgboost_grade_2$best_params)),"\n")
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 3, són:",paste(names(best_hyperparams_xgboost_grade_3$best_params),":", unlist(best_hyperparams_xgboost_grade_3$best_params)),"\n")


# Funció genèrica per obtenir la rellevància dels diferents gens o conjunts de gens a partir de cada grau histològic diferent i
# el conjunt d'hiperparàmetres òptims trobats amb la corss validation grid Search
# La variable objectiu ve marcada pel paràmetre indicador de grau i el retorn es el conjunt de resultats d'importàncies 
# per al grau concret.
xgboost_importance_by_grade <- function(X, histological_grade, best_params) {
  
  cat("Cerca d'importàncies ",histological_grade, " \n")
  
  # L'algorisme requereix d'una matriu en format DMatrix
  # Volem saber la importància de cada gen respecte cada tipus de variable objectiu, per tant, no caldrà subdidir entre train i test
  # ja que no es volen realitzar prediccions de classificació. Tot el conjunt 'X' serà d'entrenament.
  dtrain <- xgb.DMatrix(data = X, label = dataset[[histological_grade]])
  
  # Entrenem model amb la matriu i els paràmetres
  model_final <- xgb.train(params = best_params, data = dtrain, nrounds = 100)
  
  # Obtenir i retornar la importància de les característiques 
  importance <- xgb.importance(model = model_final)
  
  return(importance)
  
}

# Establim el nombre de mostres que es volen mostrar en el rànking dels millors i els pitjors
top_num = 15

# A partir de la funció definida anteriorment, obtenim les importàncies de les característiques (gens) 
# per a cada classe (grau histològic: 1, 2 i 3)
# importàncies GRAU HISTOLÒGIC 1
importance_xgboost_grade_1 <- xgboost_importance_by_grade(X,"is_grade_1", best_hyperparams_xgboost_grade_1$best_params)
# importàncies GRAU HISTOLÒGIC 2
importance_xgboost_grade_2 <- xgboost_importance_by_grade(X,"is_grade_2", best_hyperparams_xgboost_grade_2$best_params)
# importàncies GRAU HISTOLÒGIC 3
importance_xgboost_grade_3 <- xgboost_importance_by_grade(X,"is_grade_3", best_hyperparams_xgboost_grade_3$best_params)


# Mostra del TOP 15 en el rànking de millors resultats d'importància de gens o conjunts de gens per a cada grau histològic
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 1
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 1 --> Gens/s més importants:\n")
print(head(importance_xgboost_grade_1, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 2 --> Gens/s més importants:\n")
print(head(importance_xgboost_grade_2, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 3 --> Gens/s més importants:\n")
print(head(importance_xgboost_grade_3, top_num))


# Visualització de resultants emprant gràfiques verticals amb 'ggplot'
feature_importance1 <- importance_xgboost_grade_1 %>% arrange(desc(Gain))
feature_importance2 <- importance_xgboost_grade_2 %>% arrange(desc(Gain))
feature_importance3 <- importance_xgboost_grade_3 %>% arrange(desc(Gain))

plot <- ggplot(feature_importance1[1:top_num,], aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "XGBoost - Top 15 Importància Genètica (Grau Histològic 1) segons Gain", x = "Gens", y = "Gain") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_XGBoost.png"), plot = plot, width = 15, height = 10)


plot <- ggplot(feature_importance2[1:top_num,], aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  labs(title = "XGBoost - Top 15 Importància Genètica (Grau Histològic 2) segons Gain", x = "Gens", y = "Gain") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )


# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_XGBoost.png"), plot = plot, width = 15, height = 10)


plot <- ggplot(feature_importance3[1:top_num,], aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "green") +
  coord_flip() +
  labs(title = "XGBoost - Top 15 Importància Genètica (Grau Histològic 3) segons Gain", x = "Gens", y = "Gain") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )


# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_XGBoost.png"), plot = plot, width = 10, height = 10)


# Visualització de resultants emprant diagrames de punts amb 'ggplot'
plot <- ggplot(feature_importance1[1:top_num,], aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(size = 3, color = "blue") +
  labs(title = "XGBoost - Top 15 Importància Genètica (Grau Histològic 1) segons Gain", x = "Importància (Gain)", y = "Gens")
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_XGBoost_punts.png"), plot = plot, width = 15, height = 10)

plot <- ggplot(feature_importance2[1:top_num,], aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(size = 3, color = "red") +
  labs(title = "XGBoost - Top 15 Importància Genètica (Grau Histològic 2) segons Gain", x = "Importància (Gain)", y = "Gens")
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_XGBoost_punts.png"), plot = plot, width = 15, height = 10)


plot <- ggplot(feature_importance3[1:top_num,], aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(size = 3, color = "green") +
  labs(title = "XGBoost - Top 15 Importància Genètica (Grau Histològic 3) segons Gain", x = "Importància (Gain)", y = "Gens")
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_XGBoost_punts.png"), plot = plot, width = 10, height = 10)