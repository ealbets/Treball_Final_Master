#### FILE 1 PROCESSMENT: Processat i captura dels atributs identificadors dels pacients, el grau histològic i informació addicional ####

### LOAD DATA ###
# Ruta al fitxer de dades 1
file_path_1 <- "datos/data_part_1.txt"

# Read lines
lines_file_1 <- readLines(file_path_1)

print("DATASET 1, PAS 1: Filtratge i captura de camps")
####### PAS 1: Filtratge #########
# Filtrem els camps que ens interessen en el fitxer.
# Id del pacient ve a !Sample_title
sample_title_1 <- grep("^!Sample_title", lines_file_1, value = TRUE)
# Capturem informació sobre el grau histològic del tumor en '!Sample_characteristics_ch1' (només el que contingui 'tumor_grade:')
sample_characteristics_1 <- grep("^!Sample_characteristics_ch1.*tumor grade:", lines_file_1, value = TRUE)
# 'geo_accession' coné un identificador únic (en aquest cas GSM, per tant corresponent a una mostra individual que indentifica l'expressió genètica de l'individu)
sample_geo_accession_1 <- grep("^!Sample_geo_accession", lines_file_1, value = TRUE)
# Es captura el tipus d'expressió genètica que es sempre RNA
sample_type_1 <- grep("^!Sample_type", lines_file_1, value = TRUE)
# Es captura el país de residència de la pacient en 'Sample_contact_country'
sample_country_1 <- grep("^!Sample_contact_country", lines_file_1, value = TRUE)
# En 'sample_source_name_ch1' es captura la font biològica d'origen de la mostra (tumor de mama, biòpsia de tumor de mama...etc) 
sample_name_1 <- grep("^!Sample_source_name_ch1", lines_file_1, value = TRUE)


print("DATASET 1, PAS 2: Transformació de camps")
###### PAS 2: Transformació ######
# Separar els elements resultants pel delimitador "\t"
sample_title_split_1 <- strsplit(sample_title_1, "\t")[[1]]
sample_characteristics_split_1 <- strsplit(sample_characteristics_1, "\t")[[1]]
sample_geo_accession_split_1 <- strsplit(sample_geo_accession_1, "\t")[[1]]
sample_type_split_1 <- strsplit(sample_type_1, "\t")[[1]]
sample_country_split_1 <- strsplit(sample_country_1, "\t")[[1]]
sample_name_split_1 <- strsplit(sample_name_1, "\t")[[1]]

# Per a 'sample_characteristics_2' tenim dues columnes on apareix la característica 'tumor_grade' i són complementàries, és a dir,
# els índex de posicions on surt el 'tumor_grade' en la primera és on no surt a la segona i viceversa. Quan no informa del 'tumor_grade'
# en la primera, indica el 'cm_stage / er_status'. Quan no informa del 'tumor_grade' en la segona, indica el 'er_status / pr_status'. El que es farà és:
# De la primera charactersitics_ch1 ens quedem amb el 'tumor_grade' i on NO hi hagi un 'tumor_grade' posem un string buit.
sample_characteristics_1_1 <- sample_characteristics_1[1]
sample_characteristics_1_1_split_2 <- strsplit(sample_characteristics_1_1, "\t")[[1]]
# En la segona characteristics_ch1 es procedeix de la mateixa manera
sample_characteristics_1_2 <- sample_characteristics_1[2]
sample_characteristics_1_2_split_2 <- strsplit(sample_characteristics_1_2, "\t")[[1]]

# Eliminar el primer element corresponent al nom de la columna o del camp
sample_title_clean_1 <- sample_title_split_1[-1]
sample_characteristics_1_1_clean_2 <- sample_characteristics_1_1_split_2[-1]
sample_characteristics_1_2_clean_2 <- sample_characteristics_1_2_split_2[-1]
sample_geo_accession_clean_1 <- sample_geo_accession_split_1[-1]
sample_type_1 <- sample_type_split_1[-1]
sample_country_clean_1 <- sample_country_split_1[-1]
sample_name_clean_1 <- sample_name_split_1[-1]

# Eliminació de les cometes entre els elements
sample_title_clean_1 <- gsub("\"", "", sample_title_clean_1)
sample_characteristics_1_1_clean_2 <- gsub("\"", "", sample_characteristics_1_1_clean_2)
sample_characteristics_1_2_clean_2 <- gsub("\"", "", sample_characteristics_1_2_clean_2)
sample_geo_accession_clean_1 <- gsub("\"", "", sample_geo_accession_clean_1)
sample_type_clean_1 <- gsub("\"", "", sample_type_1)
sample_country_clean_1 <- gsub("\"", "", sample_country_clean_1)
sample_name_clean_1 <- gsub("\"", "", sample_name_clean_1)

# En el cas dels graus histològics, capturem només el que es descriu en 'tumor_grade:'. Si no es conté
# la cadena 'tumor_grade:', aleshores el valor en aquesta posició es deixa a buit ja que no és el grau histològic
sample_characteristics_1_1_clean_2 <- ifelse(!grepl("tumor grade:", sample_characteristics_1_1_clean_2), "", sample_characteristics_1_1_clean_2)
sample_characteristics_1_2_clean_2 <- ifelse(!grepl("tumor grade:", sample_characteristics_1_2_clean_2), "", sample_characteristics_1_2_clean_2)
# Finalment, unificarem els resultats fent que en les posicions buides quedin els 'tumor_grade' complementaris de les 2 characteristics
sample_characteristics_1_clean_final <- paste0(sample_characteristics_1_1_clean_2, sample_characteristics_1_2_clean_2)
# Extraiem l'últim caràcter i el convertim en número. Si es troba '' es substitueix per NA 
tumor_grades_1 <- as.numeric(sub(".*grade: ", "", sample_characteristics_1_clean_final))

print("DATASET 1, PAS 3: Generació del dataset parcial sense la informació de l'expressió genètica")
# Generem el dataset1 parcial sense la informació de l'expressió genètica del pacient amb les variables escollides anteriorment. El dataframe tindrà
# - TITLE (identificador del pacient amb tumor de mama, obligatori)
# - GEO_ACCESSION (identificador de la mostra del pacient, obligatori)
# - COUNTRY (país de residència del pacient, opcional)
# - TYPE (tipus estrucutral de la mostra: ARN, obligatori)
# - NAME (font biològica d'origen de la mostra, obligatori)
# - HISTOLOGICAL_GRADE (grau histològic del tumor que serà la variable objectiu obligatòria i podrà tenir tres valors: 1, 2 o 3)
dataset1_parcial <- data.frame(
  TITLE = sample_title_clean_1,
  GEO_ACCESSION = sample_geo_accession_clean_1,
  COUNTRY = sample_country_clean_1,
  TYPE = sample_type_clean_1,
  NAME = sample_name_clean_1,
  HISTOLOGICAL_GRADE = tumor_grades_1
)

head(dataset1_parcial)
tail(dataset1_parcial)

print("DATASET 1, PAS 4: Captura i càrrega de la matriu d'expressions genètiques extreta per microarrays")
### MATRIU DE DADES DE L'EXPRESSIÓ GENÈTICA ###
# En aquest fitxer els ID_REF apunten a la taula de la plataforma 'GPL570-55999'

# Ens situem en la zona de la matriu de dades, que s'inicia en '!series_matrix_table_begin' i conclou en 'series_matrix_table_end'
inici_matriu_1 <- grep("!series_matrix_table_begin", lines_file_1)
fi_matriu_1 <- grep("!series_matrix_table_end", lines_file_1)

# Capturem el contingut entre inici i fi
# Extraiem la part desitjada
matriu_expressio_genetica_1 <- lines_file_1[(inici_matriu_1 + 1):(fi_matriu_1 - 1)]

# Ho convertim a un dataframe on l'espai '\t' serà el delimitador de columna
matriu_expressio_genetica_1 <- read.table(text = matriu_expressio_genetica_1, header = TRUE, sep = "\t", quote = "\"")

# Mostrem estadístiques bàsiques
print("** DADES :\n**")
head(matriu_expressio_genetica_1,n=15)
tail(matriu_expressio_genetica_1,n=15)
print("** ESTRUCTURA:\n **")
str(matriu_expressio_genetica_1)
print("** ESTADÍSTIQUES BÀSIQUES:\n **")
summary(matriu_expressio_genetica_1)

# Obtenir dimensions de la matriu
dimensions <- dim(matriu_expressio_genetica_1)
# Mostrar el nombre de files y columnas
cat("La matriu d'expressió genètica té ", dimensions[1], "files i", dimensions[2], "columnes.\n")

print("DATASET 1, PAS 5: Normalització de dades")
## Per tal de tenir tots els valors en la mateixa escala a l'hora d'unir els dos conjunts de dades, realitzem una normalització de valors de la matriu entre 0 i 1
## Cridem la funció creada en 'general_functions' anomenada 'normalize_matrix'
# Passem els valors numèrics (no s'inclou la primera columna dels identificadors)
matriu_expressio_genetica_1_norm <- normalize_matrix(matriu_expressio_genetica_1[, 2:ncol(matriu_expressio_genetica_1)])

# Reafegim la columna original dels identificadors
matriu_expressio_genetica_1_norm <- cbind(matriu_expressio_genetica_1[, 1, drop = FALSE], matriu_expressio_genetica_1_norm)