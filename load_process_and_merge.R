if (!require("dplyr")) install.packages("dplyr")
# Load dplyr
library(dplyr)

## MATRIU EXPRESSIÓ GENÈTICA, procés obtenció GENE_SYMBOL a través de Plataforma
# Empresa de la plataforma: Affymetrix, ja conté un camp 'GENE_SYMBOL'

# Ruta a la taula de la plataforma
file_to_platform_affymetrix <- "datos/platform_microarrays/"
file_platofrm_name <- "GPL570-55999.txt"
ruta_destino <- "datos/gene_expressions_matrix/"

print("PROCESSAT PAS 1: Càrrega, processat i unió de les dades de la matriu d'expressió genètica obtinguda del conjunt de dades amb la informació proporcionada per la taula de la plataforma")

# 1. Cridem la funció 'process_platform_table' de 'general_functions' que realitza tota la carrega, processat i unió de les dades
# de la matriu d'expressió genètica obtinguda del conjunt de dades amb la informació proporcionada per la taula de la plataforma (GENE_SYMBOL)
# El nexe d'unió és el ID_REF
# Per a la matriu del conjunt de dades 1
matriu_expressio_genetica_1_final <- process_platform_table(file_to_platform_affymetrix, file_platofrm_name, matriu_expressio_genetica_1_norm)
# Per a la matriu del conjunt de dades 2
matriu_expressio_genetica_2_final <- process_platform_table(file_to_platform_affymetrix, file_platofrm_name, matriu_expressio_genetica_2_norm)


print("PROCESSAT PAS 2: Procés de neteja de les dades d'expressió genètica unides amb la informació de la plataforma que s'han obtingut en el pas anterior")

# 2. Cridem la funció 'merge_platform_with_gene_expression' que realitza el procés de neteja de les dades d'expressió genètica unides amb la
# informació de la plataforma que s'han obtingut en el pas anterior

# Per a la matriu unificada del conjunt de dades 1
print("Matriu del conjunt de dades 1:")
matriu_expressio_genetica_1_final_proc <- process_merged_gene_expression(matriu_expressio_genetica_1_final)
# Per a la matriu unificada del conjunt de dades 2
print("Matriu del conjunt de dades 2:")
matriu_expressio_genetica_2_final_proc <- process_merged_gene_expression(matriu_expressio_genetica_2_final)

# Visualizem els resultats
head(matriu_expressio_genetica_1_final_proc)
tail(matriu_expressio_genetica_1_final_proc) 

head(matriu_expressio_genetica_2_final_proc)
tail(matriu_expressio_genetica_2_final_proc)

# Obtenir dimensions de la matrius
dimensions_1 <- dim(matriu_expressio_genetica_1_final_proc)
# Calcular el número de mostres excloent la primera columna
num_mostres_1 <- dimensions_1[2] - 1  # Resta 1 per a excloure la columna GENE_SYMBOL
# Mostrar el nombre de files y columnas
cat("La matriu d'expressió genètica processada 1 té ", dimensions_1[1], "sondes genètiques i", num_mostres_1, "mostres\n")

dimensions_2 <- dim(matriu_expressio_genetica_2_final_proc)
# Calcular el número de mostres excloent la primera columna
num_mostres_2 <- dimensions_2[2] - 1  # Resta 1 per a excloure la columna GENE_SYMBOL
# Mostrar el nombre de files y columnas
cat("La matriu d'expressió genètica processada 1 té ", dimensions_2[1], "sondes genètiques i", num_mostres_2, "mostres\n")


print("PROCESSAT PAS 3: Unió de les matrius resultants en una de sola")

# Unim les dues matrius resultants de les expressions genètiques dels diferents conjunts de dades per a que quedi en una de sola
# El nexe d'unió serà GENE_SYMBOL

matriu_expressio_genetica_final_unificada <- merge(
  matriu_expressio_genetica_1_final_proc, 
  matriu_expressio_genetica_2_final_proc, 
  by = c("GENE_SYMBOL"), 
  all = TRUE
)

# Visualizem els resultats
head(matriu_expressio_genetica_final_unificada)
tail(matriu_expressio_genetica_final_unificada) 

# Obtenir dimensions de la matrius
dimensions_finals <- dim(matriu_expressio_genetica_final_unificada)
# Calcular el número de mostres excloent la primera columna
num_mostres_finals <- dimensions_finals[2] - 1  # Resta 1 per a excloure la columna GENE_SYMBOL
# Mostrar el nombre de files y columnas
cat("La matriu d'expressió genètica processada 1 té ", dimensions_finals[1], "sondes genètiques i", num_mostres_finals, "mostres\n")

print("PROCESSAT PAS 4: Bolcat de resultats en fitxer .csv")
# Escrivim els resultats a un fitxer .CSV 
# Bolquem els resultats en un fitxer .csv i elidim els 5 últims registres que no aporten cap tipus d'informació gènica
n <- nrow(matriu_expressio_genetica_final_unificada)
write_to_csv(matriu_expressio_genetica_final_unificada[1:(n - 5), ], ruta_destino, "matriu_expressio_genetica_final.csv")


## UNIÓ DATASETS parcials sense informació genètica
print("PROCESSAT PAS 5: Unió dels datasets 1 i 2 parcials (que no incluen informació genètica dels pacients")
dataset_final <- rbind(dataset1_parcial, dataset2_parcial)

head(dataset_final)
tail(dataset_final)
