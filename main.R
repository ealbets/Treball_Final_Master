# main.R
setwd("C:/Users/ernea/OneDrive/Escritorio/UOC - Data Science/SEMESTRE 5 TFM/src/Treball_Final_Master/")

# Especificar fitxer logs de sortida
sink("console/output_console.txt")


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


# TAULA RESUM

if (!require("gt")) install.packages("gt")
if (!require("webshot2")) install.packages("webshot2")
# Cargamos las librerías
library(gt)
library(webshot2)

# Cració taula de dades
taula_resum <- data.frame(
  "CARACTERÍSTICA" = c("TIPUS CONJUNT DE DADES", "TIPUS D'ALGORISME", "VARIABLE OBJECTIU", "TIPUS VARIABLE OBJECTIU", 
                       "VARIABLES INDEPENDENTS", "TIPUS VARIABLES INDEPENDENTS"),
  "DESCRIPCIÓ" = c("Supervisat", "Regressió", "HISTOLOGICAL_GRADE", 
                   "Quantitativa categòrica ordinal de 3 classes (1,2 i 3) convertida en 3 variables objectius diferents (pertinença a cada un dels graus) de tipus binari", 
                   "Tot el conjunt de columnes que representen gens o conjunts de gens", 
                   "Quantitatives contínues amb valors normalitzats compresos entre 0 i 1")
)

# Creció taula amb format gt
taula_formatted <- taula_resum %>%
  gt() %>%
  tab_header(
    title = "Característiques del Conjunt de Dades"
  ) %>%
  tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    table_body.border.bottom.color = "black",
    heading.border.bottom.color = "black"
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "aquamarine")  # Fons de color aquamarina
    ),
    locations = cells_title(groups = "title")  # Aplica estil al títol de la capçalera
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "lightgray"),  # Fons gris clar
      cell_text(weight = "bold")      # Text en negreta
    ),
    locations = cells_column_labels()  # Aplica l'estil als noms de les columnes
  )

# Guardar la taula como una imatge
output_path <- "data/output/taula_resum.png"
gtsave(taula_formatted, filename = output_path)

# Confirmar ruta de sortida
cat("La taula resum s'ha guardat a:", output_path, "\n")




##--APLICACIÓ DE TÈCNIQUES DE MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ --#
source("preparacio_dades_ml.R")
##-- TÈNIQUES BASADES EN ARBRES DE DECISIÓ I GRADIENT --##
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ I: XGBOOST --##
source("XGBoost.R")
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ II: LIGHTGBM --##
source("LightGBM.R")
##--TÈCNIQUES BASADES EN REGULARITZACIÓ --##
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ III: RIDGE REGRESSION L2 --##
source("RidgeRegression.R")
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ IV: LASSO REGRESSION L1 --##
source("LassoRegression.R")
##-- TÈCNICA MACHINE-LEARNING SUPERVISAT DE REGRESSIÓ V: ELASTIC NET --##
source("ElasticNet.R")

sink()
