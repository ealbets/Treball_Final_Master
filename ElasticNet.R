# Paquet xgboost
if (!require("glmnet")) install.packages("glmnet")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("progress")) install.packages("progress")

library(glmnet)
library(dplyr)
library(ggplot2)
library(progress)

path_images <- "images/ElasticNet/"

##-- MACHINE-LEARNING SUPERVISAT AMB REGRESSIÓ: Algorisme d'Elastic Net --##

# Variable objectiu--> grau histològic (categòrica ordinal)
# Variables independents --> conjunt de 10.600 columnes referents a gens (valors continus normalitzats entre 0 i 1)

cat("-------ELASTIC NET ALGORITHM----------\n")

# Crear variables binàries per a cada grau histològic
dataset$is_grade_1 <- ifelse(dataset$HISTOLOGICAL_GRADE == 1, 1, 0)
dataset$is_grade_2 <- ifelse(dataset$HISTOLOGICAL_GRADE == 2, 1, 0)
dataset$is_grade_3 <- ifelse(dataset$HISTOLOGICAL_GRADE == 3, 1, 0)


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

# Capturem les variables independents que són les corresponents als gens o conjunts de gens. Aquestes van de la columna 7
# a la N (10647) però n'hem d'excloure les 3 últimes variables temporals afegides en el pas anterior d'indicadors de grau específic.
# Ho convertim en matriu per a que sigui compatible amb l'algorisme.
X <- as.matrix(dataset[, 7:(ncol(dataset) - 3)])


# Definim les possibilitats d'hiperparàmetres per a cercar els òptims
# seqüència de valors alpha i valors lambda
alpha <- seq(0.1, 1, by = 0.15)
lambdes <- c(0.01, 0.05, 0.1, 0.5, 1, 10)


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


# Funció genèrica per obtenir la rellevància dels diferents gens o conjunts de gens a partir de cada grau histològic diferent i
# el conjunt d'hiperparàmetres òptims trobats amb la corss validation grid Search
# La variable objectiu ve marcada pel paràmetre indicador de grau i el retorn es el conjunt de resultats d'importàncies 
# per al grau concret.
elasticnet_importance_by_grade <- function(X, histological_grade, best_alpha, best_lambda) {
  
  cat("Cerca d'importàncies ",histological_grade, " \n")
  
  # Entrenament del model final amb els millors hiperparàmetres
  # Volem saber la importància de cada gen respecte cada tipus de variable objectiu, per tant, no caldrà subdidir entre train i test
  # ja que no es volen realitzar prediccions de classificació. Tot el conjunt 'X' serà d'entrenament
  en_model <- glmnet(
    X, dataset[[histological_grade]],
    alpha = best_alpha,
    lambda = best_lambda,
    family = "binomial"
  )
  
  # Obtenir les importàncies de les característiques a partir dels coeficients
  coeficients <- coef(en_model, s = best_lambda)
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
  
  
  
  return(feature_importance)
  
  
}

# Establim el nombre de mostres que es volen mostrar en el rànking dels millors i els pitjors
top_num = 10

# A partir de la funció definida anteriorment, obtenim les importàncies de les característiques (gens) 
# per a cada classe (grau histològic: 1, 2 i 3)
# importàncies GRAU HISTOLÒGIC 1
importance_elasticnet_grade_1 <- elasticnet_importance_by_grade(X,"is_grade_1", best_hyperparams_elasticnet_grade_1$best_alpha,best_hyperparams_elasticnet_grade_1$best_lambda)
# importàncies GRAU HISTOLÒGIC 2
importance_elasticnet_grade_2 <- elasticnet_importance_by_grade(X,"is_grade_2", best_hyperparams_elasticnet_grade_2$best_alpha,best_hyperparams_elasticnet_grade_2$best_lambda)
# importàncies GRAU HISTOLÒGIC 3
importance_elasticnet_grade_3 <- elasticnet_importance_by_grade(X,"is_grade_3", best_hyperparams_elasticnet_grade_3$best_alpha,best_hyperparams_elasticnet_grade_3$best_lambda)


# Mostra del TOP 10 en el rànking de millors resultats d'importància de gens o conjunts de gens per a cada grau histològic
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 1
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 1 --> Gens/s més importants:\n")
print(head(importance_elasticnet_grade_1, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 2
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 2 --> Gens/s més importants:\n")
print(head(importance_elasticnet_grade_2, top_num))
# TOP top_num millors importàncies genètiques GRAU HISTOLÒGIC 3
cat("Top ", top_num,  ": GRAU HISTOLÒGIC 3 --> Gens/s més importants:\n")
print(head(importance_elasticnet_grade_3, top_num))


# Visualització de resultants emprant gràfiques verticals amb 'ggplot'

# Afegim una columna per categoritzar els coeficients com a positius o negatius
importance_elasticnet_grade_1 <- importance_elasticnet_grade_1 %>%
  mutate(Signe = ifelse(Coef > 0, "Positiu (Risc)", "Negatiu (Protecció)"))

# Gràfic amb diferenciació de colors per signe del coeficient.
# blau pels coeficients positius (sentit negatiu, ja que denota risc de patiment)
# vermell pels coeficients negatius (sentit positiu, ja que denota protecció)
plot <- ggplot(importance_elasticnet_grade_1[1:top_num,], 
               aes(x = reorder(Feature, Coef), y = Coef, fill = Signe)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Positiu (Risc)" = "blue", "Negatiu (Protecció)" = "red")) +
  labs(
    title = paste0("Elastic Net - Top ", top_num, " Importància Genètica (Grau Histològic 1)"),
    subtitle = "Coeficients positius (risc) i negatius (protecció)",
    x = "Gens",
    y = "Coeficient"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank(),                               # Línies de quadrícula menors
    legend.position = "top"                                           # Llegenda a la part superior
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_ElasticNet.png"), 
       plot = plot, width = 15, height = 10)



# Afegim una columna per categoritzar els coeficients com a positius o negatius
importance_elasticnet_grade_2 <- importance_elasticnet_grade_2 %>%
  mutate(Signe = ifelse(Coef > 0, "Positiu (Risc)", "Negatiu (Protecció)"))


plot <- ggplot(importance_elasticnet_grade_2[1:top_num,], 
               aes(x = reorder(Feature, Coef), y = Coef, fill = Signe)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Positiu (Risc)" = "blue", "Negatiu (Protecció)" = "red")) +
  labs(
    title = paste0("Elastic Net - Top ", top_num, " Importància Genètica (Grau Histològic 2)"),
    subtitle = "Coeficients positius (risc) i negatius (protecció)",
    x = "Gens",
    y = "Coeficient"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank(),                               # Línies de quadrícula menors
    legend.position = "top"                                           # Llegenda a la part superior
  )


# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_ElasticNet.png"), plot = plot, width = 15, height = 10)


# Afegim una columna per categoritzar els coeficients com a positius o negatius
importance_elasticnet_grade_3 <- importance_elasticnet_grade_3 %>%
  mutate(Signe = ifelse(Coef > 0, "Positiu (Risc)", "Negatiu (Protecció)"))

plot <- ggplot(importance_elasticnet_grade_3[1:top_num,], 
               aes(x = reorder(Feature, Coef), y = Coef, fill = Signe)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Positiu (Risc)" = "blue", "Negatiu (Protecció)" = "red")) +
  labs(
    title = paste0("Elastic Net - Top ", top_num, " Importància Genètica (Grau Histològic 3)"),
    subtitle = "Coeficients positius (risc) i negatius (protecció)",
    x = "Gens",
    y = "Coeficient"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank(),                               # Línies de quadrícula menors
    legend.position = "top"                                           # Llegenda a la part superior
  )


# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_ElasticNet.png"), plot = plot, width = 10, height = 10)


# Visualització de resultants emprant diagrames de punts amb 'ggplot'
plot <- ggplot(importance_elasticnet_grade_1[1:top_num,], 
               aes(x = abs(Coef), y = reorder(Feature, abs(Coef)), color = Signe)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Positiu (Risc)" = "blue", "Negatiu (Protecció)" = "red")) +
  labs(
    title = paste0("Elastic Net - Top ", top_num, " Importància Genètica (Grau Histològic 1)"),
    subtitle = "Coeficients positius (risc) i negatius (protecció)",
    x = "Importància (|Coef|)",
    y = "Gens"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank(),                               # Línies de quadrícula menors
    legend.position = "top"                                           # Llegenda a la part superior
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_1_ElasticNet_punts.png"), plot = plot, width = 15, height = 10)


plot <- ggplot(importance_elasticnet_grade_2[1:top_num,], 
               aes(x = abs(Coef), y = reorder(Feature, abs(Coef)), color = Signe)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Positiu (Risc)" = "blue", "Negatiu (Protecció)" = "red")) +
  labs(
    title = paste0("Elastic Net - Top ", top_num, " Importància Genètica (Grau Histològic 2)"),
    subtitle = "Coeficients positius (risc) i negatius (protecció)",
    x = "Importància (|Coef|)",
    y = "Gens"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank(),                               # Línies de quadrícula menors
    legend.position = "top"                                           # Llegenda a la part superior
  )


# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_2_ElasticNet_punts.png"), plot = plot, width = 15, height = 10)


plot <- ggplot(importance_elasticnet_grade_3[1:top_num,], 
               aes(x = abs(Coef), y = reorder(Feature, abs(Coef)), color = Signe)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Positiu (Risc)" = "blue", "Negatiu (Protecció)" = "red")) +
  labs(
    title = paste0("Elastic Net - Top ", top_num, " Importància Genètica (Grau Histològic 3)"),
    subtitle = "Coeficients positius (risc) i negatius (protecció)",
    x = "Importància (|Coef|)",
    y = "Gens"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
    plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
    panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
    panel.grid.minor = element_blank(),                               # Línies de quadrícula menors
    legend.position = "top"                                           # Llegenda a la part superior
  )

# Guardem el gràfic generat
ggsave(filename = paste0(path_images,"Grau_Histologic_3_ElasticNet_punts.png"), plot = plot, width = 10, height = 10)