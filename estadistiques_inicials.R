cat(" ESTADÍSTIQUES BÀSIQUES I ANÀLISI DEL CONJUNT DE DADES FINAL GENERAT: \n")
dataset <- dataset_final_unificat

# Estadístiques descriptives i resums de columnes
# Excloem totes les columnes referents a gens o conjunt de gens
cat("Estadístiques \n")
summary(dataset[1:6])
cat("Estructura interna: \n")
str(dataset)

##-- ESTADÍSTIQUES PER GRAU HISTOLÒGIC, VARIABLE OBJECTIU --##
cat(" ESTADÍSTIQUES PER GRAU HISTOLÒGIC, VARIABLE OBJECTIU: \n")
# Pas 1: Comptar freqüències de cada categoria en la columna histologic_grade
counts <- table(dataset$HISTOLOGICAL_GRADE)
cat("Recompte de registres per a cada grau histològic: \n")
# Convertim els resultats a un data frame
counts_df <- as.data.frame(counts)
# Renombrament columnes
colnames(counts_df) <- c("GRAU HISTOLOGIC", "NUM. REGISTRES")
# Visualització taula
print(counts_df)

# Pas 2: Càlcul de percentatges
total_rows = nrow(dataset)
percentatges <- (counts / total_rows) * 100

cat("Percentatges de representació de registres per a cada grau histològic: \n")
# Convertim els resultats a un data frame
percentatges_df <- as.data.frame(percentatges)
# Renombrament columnes
colnames(percentatges_df) <- c("GRAU HISTOLOGIC", "% REGISTRES")
# Visualització taula
print(percentatges_df)

# Loop que recorre les files de la taula de percentatges
for(i in 1:nrow(percentatges_df)) {
  cat("Les pacients amb tumor de mama de grau histològic ", 
      percentatges_df$`GRAU HISTOLOGIC`[i], 
      " representen el ", 
      percentatges_df$`% REGISTRES`[i], 
      "% del total\n")
}

# Pas 3: Mostra d'histogrames segons variable objectiu: GRAU HISTOLOGIC

#par(mar = c(4, 4, 2, 1))  # Ajuste de márgenes (opcional)
cat("Generació d'histograma que mostra les freqüències de pacients de tumors de mama en cada un dels tres graus histològics \n")
# Creació histograma de freqüencies de pacients segons els diferents graus histològics
hist_obj  <- hist(dataset$HISTOLOGICAL_GRADE, 
                  breaks = seq(0.5, 3.5, by = 1),  # Punts de tall
                  main = "Distribució de freqüències tumorals segons els Graus Histològics", 
                  xlab = "GRAU HISTOLÒGIC", 
                  ylab = "Freqüència", 
                  col = "red", 
                  border = "black", 
                  xaxt = "n")  # Evita los valors per defecte en eix X

# Inclusió de labels en l'eix X
axis(1, at = 1:3, labels = c("GRAU 1", "GRAU 2", "GRAU 3"))

# Inclusió informació percentatges en les barres
text(x = hist_obj$mids,  # Posicions X
     y = hist_obj$counts + 1,  
     labels = paste0(round(percentatges, 1), "%"),  # Format de percentatges
     cex = 1.2,  # Mida text
     col = "black")  # Color text


cat("RESULTATS:\n")
cat("Els casos de tumors de grau histològic 2 són els més freqüents en el nostre conjunt de dades destacant amb un 48% i 96 registres.\n")
cat("El segueix de prop els tumors de grau histològic 3 amb un 40% (80 registres) i, en últim lloc, i amb una diferència considerablement pronunciada, tenim els casos de tumors de grau histològic 1 amb un 12% i 24 registres.\n\n")

cat("S'aprecia un desbalanceig considerable del nombre de tumors de grau histològic 1 (12%) respecte als de grau 2 (48%) i 3 (40%). Aquest fet pot afectar a l'hora de definir amb precisió els gens o conjunts de gens més involucrats, ja que, a excepció de la resta, es posseeixen pocs casos per provar.\n")
cat("Pel que respecta a la diferència entre el nombre d'elements de grau histològic 2 amb el d'elements de grau histològic 3, aquesta no és considerable ni molt significativa.\n")
cat("Per tant, en aquests dos casos no hi ha un desbalanceig molt pronunciat i els resultats són bastant equitatius.\n")
cat("El nombre de casos és més nombrós i, en conseqüència, els resultats poden ser més precisos i estar més ben ajustats.\n")