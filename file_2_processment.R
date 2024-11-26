#### FILE 3 PROCESSMENT ####

### LOAD DATA ###
# Path to file
file_path_2 <- "datos/data_part_2.txt"

# Read lines
lines_file_2 <- readLines(file_path_2)

print("DATASET 2, PAS 1: Filtratge i captura de camps")
####### FILTERING #########
# Filtrem els camps que ens interessen en cada fitxer.
# Id del pacient ve a !Sample_title 
sample_title_2 <- grep("^!Sample_title", lines_file_2, value = TRUE)
# Capturem informació sobre el grau histològic del tumor en '!Sample_characteristics_ch1' (només el que contingui 'grade histo:')
sample_characteristics_2 <- grep("^!Sample_characteristics_ch1.*grade histo:", lines_file_2, value = TRUE)
# 'geo_accession' coné un identificador únic (en aquest cas GSM, per tant corresponent a una mostra individual que indentifica l'expressió genètica de l'individu)
sample_geo_accession_2 <- grep("^!Sample_geo_accession", lines_file_2, value = TRUE)
# Es captura el tipus d'expressió genètica que es sempre RNA
sample_type_2 <- grep("^!Sample_type", lines_file_2, value = TRUE)
# Es captura el país de residència de la pacient en 'Sample_contact_country'
sample_country_2 <- grep("^!Sample_contact_country", lines_file_2, value = TRUE)
# En 'sample_source_name_ch1' es captura la font biològica d'origen de la mostra (tumor de mama, biòpsia de tumor de mama...etc) 
sample_name_2 <- grep("^!Sample_source_name_ch1", lines_file_2, value = TRUE)


print("DATASET 2, PAS 2: Transformació de camps")
###### TRANSFORMATION ######
# Separar els elements pel delimitador "\t"
sample_title_split_2 <- strsplit(sample_title_2, "\t")[[1]]
sample_characteristics_split_2 <- strsplit(sample_characteristics_2, "\t")[[1]]
sample_geo_accession_split_2 <- strsplit(sample_geo_accession_2, "\t")[[1]]
sample_type_split_2 <- strsplit(sample_type_2, "\t")[[1]]
sample_country_split_2 <- strsplit(sample_country_2, "\t")[[1]]
sample_name_split_2 <- strsplit(sample_name_2, "\t")[[1]]

# Eliminar el primer element corresponent al nom de la columna o del camp
sample_title_clean_2 <- sample_title_split_2[-1]
sample_characteristics_clean_2 <- sample_characteristics_split_2[-1]
sample_geo_accession_clean_2 <- sample_geo_accession_split_2[-1]
sample_type_clean_2 <- sample_type_split_2[-1]
sample_country_clean_2 <- sample_country_split_2[-1]
sample_name_clean_2 <- sample_name_split_2[-1]

# Eliminació de les cometes entre els elements
sample_title_clean_2 <- gsub("\"", "", sample_title_clean_2)
sample_characteristics_clean_2 <- gsub("\"", "", sample_characteristics_clean_2)
sample_geo_accession_clean_2 <- gsub("\"", "", sample_geo_accession_clean_2)
sample_type_clean_2 <- gsub("\"", "", sample_type_clean_2)
sample_country_clean_2 <- gsub("\"", "", sample_country_clean_2)
sample_name_clean_2 <- gsub("\"", "", sample_name_clean_2)

# En el cas dels graus histològics, capturem només l'element numèric que ens interessa que es troba després de 'grade histo:'
tumor_grades_2 <- as.numeric(sub(".*histo: ", "", sample_characteristics_clean_2))

print("DATASET 2, PAS 3: Generació del dataset parcial sense la informació de l'expressió genètica")
# Generem el dataset2 parcial sense la informació de l'expressió genètica del pacient amb les variables escollides anteriorment. El dataframe tindrà
# - TITLE (identificador del pacient amb tumor de mama, obligatori)
# - GEO_ACCESSION (identificador de la mostra del pacient, obligatori)
# - COUNTRY (país de residència del pacient, opcional)
# - TYPE (tipus estrucutral de la mostra: ARN, obligatori)
# - NAME (font biològica d'origen de la mostra, obligatori)
# - HISTOLOGICAL_GRADE (grau histològic del tumor que serà la variable objectiu obligatòria i podrà tenir tres valors: 1, 2 o 3)
dataset2_parcial <- data.frame(
  TITLE = sample_title_clean_2,
  GEO_ACCESSION = sample_geo_accession_clean_2,
  COUNTRY = sample_country_clean_2,
  TYPE = sample_type_clean_2,
  NAME = sample_name_clean_2,
  HISTOLOGICAL_GRADE = tumor_grades_2
)

head(dataset2_parcial)
tail(dataset2_parcial)

print("DATASET 2, PAS 4: Captura i càrrega de la matriu d'expressions genètiques extreta per microarrays")
### MATRIU DE DADES DE L'EXPRESSIÓ GENÈTICA ###
# En aquest fitxer els ID_REF apunten a la taula de la plataforma 'GPL570-55999'

# Ens situem en la zona de la matriu de dades, que s'inicia en '!series_matrix_table_begin' i conclou en 'series_matrix_table_end'
inici_matriu_2 <- grep("!series_matrix_table_begin", lines_file_2)
fi_matriu_2 <- grep("!series_matrix_table_end", lines_file_2)

# Capturem el contingut entre inici i fi
# Extraiem la part desitjada
matriu_expressio_genetica_2 <- lines_file_2[(inici_matriu_2 + 1):(fi_matriu_2 - 1)]

# Ho convertim a un dataframe on l'espai '\t' serà el delimitador de columna
matriu_expressio_genetica_2 <- read.table(text = matriu_expressio_genetica_2, header = TRUE, sep = "\t", quote = "\"")

# Mostrem estadístiques bàsiques
print("** DADES :\n**")
head(matriu_expressio_genetica_2,n=15)
tail(matriu_expressio_genetica_2,n=15)
print("** ESTRUCTURA:\n **")
str(matriu_expressio_genetica_2)
print("** ESTADÍSTIQUES BÀSIQUES:\n **")
summary(matriu_expressio_genetica_2)

# Obtenir dimensions de la matriu
dimensions <- dim(matriu_expressio_genetica_2)
# Mostrar el nombre de files y columnas
cat("La matriu d'expressió genètica té ", dimensions[1], "files i", dimensions[2], "columnes.\n")

print("DATASET 2, PAS 5: Normalització de dades")
## Per tal de tenir tots els valors en la mateixa escala a l'hora d'unir els dos conjunts de dades, realitzem una normalització de valors de la matriu entre 0 i 1
## Cridem la funció creada en 'general_functions' anomenada 'normalize_matrix'
# Passem els valors numèrics (no s'inclou la primera columna dels identificadors)
matriu_expressio_genetica_2_norm <- normalize_matrix(matriu_expressio_genetica_2[, 2:ncol(matriu_expressio_genetica_2)])

# Reafegim la columna original dels identificadors
matriu_expressio_genetica_2_norm <- cbind(matriu_expressio_genetica_2[, 1, drop = FALSE], matriu_expressio_genetica_2_norm)