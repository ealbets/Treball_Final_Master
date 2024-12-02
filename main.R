# main.R
setwd("C:/Users/ernea/OneDrive/Escritorio/UOC - Data Science/SEMESTRE 5 TFM/src/Treball_Final_Master/")

##-- CARREGA, PROCESSAMENT, SELECCIÓ, TRANSFORMACIÓ, MUNTATGE I GENERACIÓ DEL CONJUNT DE DADES INICIAL--##
# Pas 1: Funcions d'ús general
source("general_functions.R")
# Pas 2: Processament fitxer dataset1
source("file_1_processment.R")
# Pas 3: Processamenr fitxer dataset2
source("file_2_processment.R")
# Pas 4: Càrrega, processament i fusió de matrius d'expressió genètica
source("load_process_and_merge.R")
# Pas 5: Selecció de caraterístiques, fusió i muntatge i generació del conjunt de dades final
source("seleccio_caracteristiques_transposicio_unio.R")
# Optional: print a message to confirm completion

##-- ESTADÍSTIQUES BÀSIQUES I ANÀLISI DEL CONJUNT DE DADES FINAL GENERAT --##
source("estadistiques_inicials.R")



##--APLICACIÓ DE TÈCNIQUES DE MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ --#

install.packages("knitr")
library(knitr)

tabla_info <- data.frame(
  "CARACTERÍSTICA" = c("TIPUS CONJUNT DE DADES", "TIPUS D'ALGORISME", "VARIABLE OBJECTIU", "TIPUS VARIABLE OBJECTIU", 
                       "VARIABLES INDEPENDENTS", "TIPUS VARIABLES INDEPENDENTS"),
  "DESCRIPCIÓ" = c("Supervisat", "Regressió", "HISTOLOGICAL_GRADE", 
                    "Quantitativa categòrica ordinal de 3 classes (1,2 i 3)", 
                    "Tot el conjunt de columnes que representen gens o conjunts de gens", 
                    "Quantitatives contínues amb valors normalitzats compresos entre 0 i 1")
)

# Mostrem la taula
kable(tabla_info, col.names = c("CARACTERÍSTICA", "DESCRIPCIÓ"), align = "l")

##-- TÈNIQUES BASADES EN ARBRES DE DECISIÓ I GRADIENBT --##
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ III: RANDOM FOREST --##
source("RandomForest.R")
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ I: XGBOOST --##
source("XGBoost.R")
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ II: LIGHTGBM --##
source("LightGBM.R")

##--TÈCNIQUES BASADES EN REGULARITZACIÓ --##
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ IV: RIDGE REGRESSION L2 --##
source("RidgeRegression_L2.R")
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ V: LASSO REGRESSION L1 --##
source("LassoRegression_L1.R")
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ VI: ELASTIC NET --##
source("ElasticNet.R")