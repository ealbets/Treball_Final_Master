## SELECCIÓ I REDUCCIÓ DE CARACTERÍSTIQUES ##

cat("SELECCIÓ DE CARACTERÍSTIQUES a partir de la desviació estàndard\n")

# 1) En primer lloc, es calcularà la mitjana aritmètica dels valors de totes les mostres per a cada una de les files que representen els gens i es guardaran en
# una columna temporal 'Average'

cat("Càlculs desviació estàndard de cada fila:\n")
# Capturem en una variable la matriu d'expressió genètica sense la columna del GENE_SYMBOL
gene_expression_matrix <- matriu_expressio_genetica_final_unificada[, -1]
# Realitzem el càlcul de la mitjana aritmètica de cada fila a partir dels valors de les mostres i ho emmagatzemem en una columna temporal 'Average'
matriu_expressio_genetica_final_unificada$SD <- apply(gene_expression_matrix, 1, sd, na.rm = TRUE)
head(matriu_expressio_genetica_final_unificada)
tail(matriu_expressio_genetica_final_unificada)

# 2) Seguidament, s'estabilarà el llindar que servirà de filtre per determinar quines files (gens) tenen una rellevància suficientment significativa per ser
# tinguts en compte en el procés de Machine-Learning supervisat de regressió posterior i, pel contrari, quines d'elles poden ser eliminades degut als baixos valors
# de totes les seves mostres.

# Calculem la mitjana de totes les desviacions estàndard obtingudes per tenir una referència del promig general
total_sd_average <- mean(matriu_expressio_genetica_final_unificada$SD, na.rm = TRUE)
# Establim el llindar com la mitjana aritmètica global de les desviacions estàndard
llindar <- total_sd_average
cat("El llindar del filtrat és ", llindar, "\n")

cat("Aplicació de filtre: \n")
# 3) Finalment, apliquem el filtre a partir del llindar ontingut. D'aquesta manera, totes les files (gens) que tinguin una desviació estàndard superior o igual
# al llindar establert, es consideraran rellevants pels futurs models de ML. En canvi, totes aquelles que no hi arribin seran eliminades ja que significa que la majoria
# de valors en les seves mostres són similars i amb poca variabilitat i, en conseqüència, manquen de la suficient importància a l'hora de definir la variable objectiu.
matriu_expressio_genetica_final_unificada_filtrada <- matriu_expressio_genetica_final_unificada[matriu_expressio_genetica_final_unificada$SD >= total_sd_average,]

cat("Nombre de gens originalment:", nrow(matriu_expressio_genetica_final_unificada), "\n")
cat("Nombre de gens seleccionats:", nrow(matriu_expressio_genetica_final_unificada_filtrada), "\n")
cat("La matriu d'expressió genètica s'ha reduït de" ,nrow(matriu_expressio_genetica_final_unificada), " a " ,nrow(matriu_expressio_genetica_final_unificada_filtrada), " característiques.\n")

## MUNTATGE FINAL DEL CONJUNT DE DADES
cat("MUNTATGE FINAL DEL CONJUNT DE DADES\n")
# Un cop aplicada la selecció de característiques anterior, ens queda preparar el conjunt de dades per a que sigui compatible amb els models ML futurs. Caldrà realitzar les
# següents operacions:

# 1) Transposició de la matriu d'expressions genètiques: El conjunt de dades final haurà de tenir les dades gèniques en columnes i no en files, és a dir, cada fila representarà
# la informació d'un pacient diferent, els valors de les seves mostres per a cada un dels gens o conjunt de gens filtrats i finalment, el grau histològic del tumor.
# Per tant, per tal de poder tenir les dades de forma que cada fila representi un pacient caldrà transposar la matriu final d'expressions genètiques i capgirar les
# files per les columnes (els gens ara seran columnes i les mostres dels pacients seran files).

# Transposició matriu expressions genètiques
cat("Transposició matriu d'expressions genètiques\n")
# Eliminem la columna temporal 'Average' necessària en el punt anterior
matriu_expressio_genetica_final_unificada_filtrada <- matriu_expressio_genetica_final_unificada_filtrada[, !names(matriu_expressio_genetica_final_unificada_filtrada) %in% "SD"]
# Canviar el nom de la columna GENE_SYMBOL i posar-los com a nomsn (identificadors) de columna GENE_SYMBOL
rownames(matriu_expressio_genetica_final_unificada_filtrada) <- matriu_expressio_genetica_final_unificada_filtrada$GENE_SYMBOL
# Eliminació de la columna 'GENE_SYMBOL' ja que ara hi haurà tantes columnes com diferents GENE_SYMBOLS únics hi hagi
matriu_expressio_genetica_final_unificada_filtrada <- matriu_expressio_genetica_final_unificada_filtrada[, -1]  
# Transpoció de la matriu per a que les files siguin mostres de pacients i les columnes siguin els gens o conjunts de gens
matriu_expressio_genetica_final_unificada_filtrada_t <- t(matriu_expressio_genetica_final_unificada_filtrada)
# Conversió a dataframe
matriu_expressio_genetica_final_unificada_filtrada_t <- as.data.frame(matriu_expressio_genetica_final_unificada_filtrada_t)
# Finalment, creem la columna 'GEO_ACCESSION' que es correspon amb l'identificador de les mostres dels pacients i hi traslladem els GEO_ACCESSION que es troben ara com a
# identificadors de les files
matriu_expressio_genetica_final_unificada_filtrada_t$GEO_ACCESSION <- rownames(matriu_expressio_genetica_final_unificada_filtrada_t)

cat("Unió de les expressions genètiques amb la informació dels pacients i el grau histològic del tumor\n")
# 2) Unió de la matriu d'expressions genètiques unificada i transposada, amb el dataset que conté la informació dels pacients i els graus histològics dels pacients a partir
# del GEO_ACCESSION com a variable relacionant.
dataset_final_unificat <- merge(dataset_final, matriu_expressio_genetica_final_unificada_filtrada_t, by = "GEO_ACCESSION", all.x = TRUE)

# Eliminem els registres que no tenen la informació mínima obligatòria requerida: grau histològic
dataset_final_unificat <- dataset_final_unificat[!is.na(dataset_final_unificat$HISTOLOGICAL_GRADE), ]

# Obtenir dimensions finals del dataset
dimensions_finals_dataset <- dim(dataset_final_unificat)
# Mostrar el nombre de files y columnas
cat("El conjunt de dades final i definitiu té ", dimensions_finals_dataset[1], "files (pacients) i", dimensions_finals_dataset[2], "columnes\n")

# Finalment, escribim els resultats finals i definitius del conjunt de dades en un .csv
write_to_csv(dataset_final_unificat, ruta_destino, "dataset_final.csv")
