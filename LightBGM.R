# Paquet xgboost
if (!require("lightgbm")) install.packages("lightgbm")


library(xgboost)
library(dplyr)
library(ggplot2)
library(progress)
library(reshape2)

path_images <- "images/LightBGM/"

##-- MACHINE-LEARNING SUPERVISAT AMB REGRESSIÓ: Algorisme de LightBGM --##

cat("-------LIGHTBGM ALGORITHM----------\n")

# Variable objectiu--> grau histològic (categòrica ordinal)
# Variables independents --> conjunt de 10.600 columnes referents a gens (valors continus normalitzats entre 0 i 1)

# Crear variables binàries per a cada grau histològic
dataset$is_grade_1 <- ifelse(dataset$HISTOLOGICAL_GRADE == 1, 1, 0)
dataset$is_grade_2 <- ifelse(dataset$HISTOLOGICAL_GRADE == 2, 1, 0)
dataset$is_grade_3 <- ifelse(dataset$HISTOLOGICAL_GRADE == 3, 1, 0)


# Funció genèrica per a trobar els hiperpàrmatres òptims (aquells que maximitzen LightBGM ) emprant la validació creuada
# i la graella d'opcions de param_grid
hyperparams_search_lgbm <- function(X, y, param_grid) {
  cat("Cerca d'hiperparàmetres òptims: \n")
  
  # Inicialitzar les variables per guardar el millor AUC i els seus hiperparàmetres
  best_auc <- -Inf
  best_params <- NULL
  
  # Barra de progrés per a tenir una noció del temps
  pb <- progress_bar$new(
    format = " Avaluant hiperparàmetres :current/:total (:percent) Temps restant: \n",
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
    cat("\n Avaluant hiperparàmetres LightBGM en ", y ," : ", 
        paste(names(params), unlist(params), sep = "=", collapse = ", "), "\n")
    
    
    # Validació creuada emprant .cv 
    cv_results <- lgb.cv(
      params = params,
      data = lgb.Dataset(data = X, label = dataset[[y]]),
      nrounds = 100, # iteracions
      nfold = 3,
      eval = "auc",
      verbose = -1
    )
    
    # Obtenir el millor AUC per aquests hiperparàmetres
    current_auc <- max(unlist(cv_results$record_evals$valid$auc$eval))
    
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

# Degut a les limitacions i les restriccions de l'algorisme LightBGM amb els noms de les columnes que actuen com a característiques
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
best_hyperparams_lbgm_grade_1 <- hyperparams_search_lgbm(X,"is_grade_1",param_grid_lgbm)
best_hyperparams_lbgm_grade_2 <- hyperparams_search_lgbm(X,"is_grade_2",param_grid_lgbm)
best_hyperparams_lbgm_grade_3 <- hyperparams_search_lgbm(X,"is_grade_3",param_grid_lgbm)

# Mostra els millors hiperparàmetres per a cada grau histològic
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 1, són:",paste(names(best_hyperparams_lbgm_grade_1$best_params),":", unlist(best_hyperparams_lbgm_grade_1$best_params)),"\n")
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 2, són:",paste(names(best_hyperparams_lbgm_grade_2$best_params),":", unlist(best_hyperparams_lbgm_grade_2$best_params)),"\n")
cat("Els millors hiperparàmetres trobats per validació creuada sobre l'avaluació del grau histològic 3, són:",paste(names(best_hyperparams_lbgm_grade_3$best_params),":", unlist(best_hyperparams_lbgm_grade_3$best_params)),"\n")


# Funció genèrica per obtenir la rellevància dels diferents gens o conjunts de gens a partir de cada grau histològic diferent i
# el conjunt d'hiperparàmetres òptims trobats amb la corss validation grid Search
# La variable objectiu ve marcada pel paràmetre indicador de grau i el retorn es el conjunt de resultats d'importàncies 
# per al grau concret.
lgbm_importance_by_grade <- function(X, histological_grade, best_params) {
  
  cat("Cerca d'importàncies ",histological_grade, " \n")
  
  # L'algorisme requereix un Dataset especial per LightGBM
  # Volem saber la importància de cada gen respecte cada tipus de variable objectiu, per tant, no caldrà subdidir entre train i test
  # ja que no es volen realitzar prediccions de classificació. Tot el conjunt 'X' serà d'entrenament.
  dtrain <- lgb.Dataset(data = X, label = dataset[[histological_grade]])
  
  # Entrenem model amb la matriu i els paràmetres
  model_final <- lgb.train(params = best_params, data = dtrain, nrounds = 100, verbose = 1)
  
  # Obtenir i retornar la importància de les característiques 
  importance <- lgb.importance(model = model_final)
  
  return(importance)
  
}

# Establim el nombre de mostres que es volen mostrar en el rànking dels millors i els pitjors
top_num = 15

# A partir de la funció definida anteriorment, obtenim les importàncies de les característiques (gens) 
# per a cada classe (grau histològic: 1, 2 i 3)
# importàncies GRAU HISTOLÒGIC 1
importance_lbgm_grade_1 <- lgbm_importance_by_grade(X, "is_grade_1", best_hyperparams_lbgm_grade_1$best_params)
# importàncies GRAU HISTOLÒGIC 2
importance_lbgm_grade_2 <- lgbm_importance_by_grade(X, "is_grade_2", best_hyperparams_lbgm_grade_2$best_params)
# importàncies GRAU HISTOLÒGIC 3
importance_lbgm_grade_3 <- lgbm_importance_by_grade(X, "is_grade_3", best_hyperparams_lbgm_grade_3$best_params)


# Invertir el mapeig inicial del gens / característiques: de noms generats a noms originals
reversed_map_lgbm <- setNames(names(column_map_lgbm), column_map_lgbm)
# Reemplaçar els noms a la columna 'Feature' amb els noms originals
importance_lbgm_grade_1$Feature <- reversed_map[importance_lbgm_grade_1$Feature]
importance_lbgm_grade_2$Feature <- reversed_map[importance_lbgm_grade_2$Feature]
importance_lbgm_grade_3$Feature <- reversed_map[importance_lbgm_grade_3$Feature]

# Mostra del TOP 15 en el rànking de millors resultats d'importància de gens o conjunts de gens per a cada grau histològic
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 1
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 1 --> Gens/s més importants:\n")
print(head(importance_lbgm_grade_1, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 2 --> Gens/s més importants:\n")
print(head(importance_lbgm_grade_2, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 3 --> Gens/s més importants:\n")
print(head(importance_lbgm_grade_3, top_num))


# Visualització de resultants emprant gràfiques verticals amb 'ggplot'
feature_lbgm_importance1 <- importance_lbgm_grade_1 %>% arrange(desc(Gain))
feature_lbgm_importance2 <- importance_lbgm_grade_2 %>% arrange(desc(Gain))
feature_lbgm_importance3 <- importance_lbgm_grade_3 %>% arrange(desc(Gain))


plot <- ggplot(feature_lbgm_importance1[1:top_num,], aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "LightBGM - Top 15 Importància Genètica (Grau Histològic 1) segons Gain", x = "Gens", y = "Gain") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_LightBGM.png"), plot = plot, width = 15, height = 10)


plot <- ggplot(feature_lbgm_importance2[1:top_num,], aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  labs(title = "LightBGM - Top 15 Importància Genètica (Grau Histològic 2) segons Gain", x = "Gens", y = "Gain") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_LightBGM.png"), plot = plot, width = 15, height = 10)


plot <- ggplot(feature_lbgm_importance3[1:top_num,], aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "green") +
  coord_flip() +
  labs(title = "LightBGM - Top 15 Importància Genètica (Grau Histològic 3) segons Gain", x = "Gens", y = "Gain") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_LightBGM.png"), plot = plot, width = 15, height = 10)


# Visualització de resultants emprant diagrames de punts amb 'ggplot'
plot <- ggplot(feature_lbgm_importance1[1:top_num,], aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(size = 3, color = "blue") +
  labs(title = "LightBGM - Top 15 Importància Genètica (Grau Histològic 1) segons Gain", x = "Importància (Gain)", y = "Gens")
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )

  
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_LightBGM_punts.png"), plot = plot, width = 15, height = 10)


plot <- ggplot(feature_lbgm_importance2[1:top_num,], aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(size = 3, color = "red") +
  labs(title = "LightBGM - Top 15 Importància Genètica (Grau Histològic 2) segons Gain", x = "Importància (Gain)", y = "Gens")
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )

  
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_LightBGM_punts.png"), plot = plot, width = 15, height = 10)
  
  
plot <-ggplot(feature_lbgm_importance3[1:top_num,], aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_point(size = 3, color = "green") +
  labs(title = "LightBGM - Top 15 Importància Genètica (Grau Histològic 3) segons Gain", x = "Importància (Gain)", y = "Gens")
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank()                                # Línies de quadrícula menors
  )
  
  
# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_LightBGM_punts.png"), plot = plot, width = 15, height = 10)