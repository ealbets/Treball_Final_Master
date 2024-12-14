## FUNCIONS GENÈRIQUES
if (!require("caret")) install.packages("caret")

library(ggplot2)
library(caret)

## Funció de càrrega, processament i unió amb les dades gèniques de la plataforma
# 1. Carrega i llegeix el contingut de la taula de la ruta especificada
# 2. Captura i reanomena la variable que ens interessa, l'odentificador i el símbol del gen: GENE_SYMBOL
# 3. Uneix la matriu d'expressions gèniques passada per paràmetre amb el seus GENE_SYMBOL que figuren en la taula de la plataforma emprant l'ID_REF com a nexe.
process_platform_table <- function(origin_path, platform_name, gene_expression_matrix) {
  
  # filePath
  filePath <- paste0(origin_path, platform_name)
  # Read lines
  paltform_file <- readLines(filePath)
  paltform_file <- read.table(text = paltform_file, header = TRUE, sep = "\t", quote = "\"")
  
  # Selecció només de la informació de l'ID i el Gene.Symbol
  gene_info_selected <- paltform_file[, c("ID", "Gene.Symbol")]
  # Substituir '///' per ', ' en la columna Gene.Symbol, pels que tenen més d'un símbol
  gene_info_selected$Gene.Symbol <- gsub("///", ", ", gene_info_selected$Gene.Symbol)
  # Reanomenem GEN_SYMBOL
  colnames(gene_info_selected)[colnames(gene_info_selected) == "Gene.Symbol"] <- "GENE_SYMBOL"
  
  # Unió de les dues taules emprant la columna ID_REF de matriu_expressio_genetica_4 y ID de gene_info_selected
  matriu_genetica_final <- merge(gene_expression_matrix, gene_info_selected[, c("ID", "GENE_SYMBOL")],
                          by.x = "ID_REF", by.y = "ID", all.x = TRUE)
  
  
  # Reorganitzar les columnes per a que ID i GENE_SYMBOL apareixin primer
  matriu_genetica_final <- matriu_genetica_final[, c("ID_REF", "GENE_SYMBOL", setdiff(names(matriu_genetica_final), c("ID_REF", "GENE_SYMBOL")))]
  
  return(matriu_genetica_final)

}

## Funció de processament i neteja de la matriu final d'expressió genètica unida a les dades de la plataforma
# 1. Comprova si existeixen registres repetits que coincideixen amb el mateix GENE_SYMBOL. En cas afirmatiu:
#  - S'observa si existeixen files l'identificador de les quals no es troba present entre els existents de la plataforma i per tant, no tenen relació. En cas afirmatiu:
#     * Els eliminem
#  - S'agrupen totes les files per GENE_SYMBOL i a les repetides s'assigna el valor de la mitja aritmètica com a valor resultant

process_merged_gene_expression <- function(gene_expression_matrix) {
  
  # Comprovació de repeticions en la columna de GENE_SYMBOL per tal d'aplicar la mitja en cas afirmatiu
  repeated_genes <- gene_expression_matrix %>%
    group_by(GENE_SYMBOL) %>%
    summarise(count = n()) %>%
    filter(count > 1)

  if (nrow(repeated_genes) > 0) {
    cat("S'han trobat símbols de gens repetits en la columna GENE_SYMBOL, procedim a realitzar les mitjanes aritmètiques dels resultats repetits\n")
  
    # Comprovem si existeixen registres sense relació amb l'identificador de la plataforma
    missing_genes <- gene_expression_matrix %>%
      filter(is.na(GENE_SYMBOL) | GENE_SYMBOL == "")
  
    if (nrow(missing_genes) > 0) {
      cat("Tenim", nrow(missing_genes)," registres NO associats a cap identificador de la plataforma que seran eliminats.")
    
    # Eliminació registres no associats a cap identificador de la plataforma
    gene_expression_matrix <- gene_expression_matrix %>%
      filter(!is.na(GENE_SYMBOL) & GENE_SYMBOL != "")
    
      cat("Els registres sense identificador de plataforma s'han eliminat correctament.\n")
    
  }  else {
      cat("Tots els registres estan associats amb un identificador.\n")
  }
  
    # Realitzem l'agrupació per columna GENE_SYMBOL i la mitjana aritmètica dels valors que coincideixen amb el mateix GENE_SYMBOL, d'aquesta manera
    # - Desapareixen els registres amb gene_symbol repetit quedant-ne només un com amb els valors mitjos de tots els coincidents com a resultant
    gene_expression_matrix <- gene_expression_matrix %>%
      group_by(GENE_SYMBOL) %>%       # Agrupem per GENE_SYMBOL
      summarise(across(starts_with("GSM"), ~ round(mean(.x, na.rm = TRUE), 3))) %>%  # Calcula la mitja de les columnes que comencen per "GSM" (geo accession)
      ungroup()
  
  } else {
    cat("No s'han trobat símbols de gens repetits en la columna GENE_SYMBOL.\n")
  }
  
  return(gene_expression_matrix)
  
}

## Funció de bolcatge dels resultats en un fitxer .csv a una ruta destí
write_to_csv <- function(gene_expression_matrix_processed, destination_path, filename) {
  
  # Bolquem els resultats en un fitxer .csv
  ruta_destino <- paste0(destination_path,filename)
  write.csv(gene_expression_matrix_processed, ruta_destino, row.names = FALSE)
  
}


## Funció de normalització
# Realitza la normalització de tots els valors continguts en una matriu de tipus llista @matrix passada per paràmetre. 
# La normalització posa tots els valors en un interval entre 0 i 1 tenint en compte el màxim i el mínim
# Funció per normalitzar una matriu entre 0 i 1
normalize_matrix <- function(matrix_list) {
  # Comprovar que la matriu d'entrada és una llista
  if (!is.list(matrix_list)) {
    stop("L'entrada ha de ser una llista.")
  }
  
  # Aplanar la llista en un sol vector per calcular el mínim i el màxim
  vector_values <- unlist(matrix_list)
  min_val <- min(vector_values, na.rm = TRUE)
  max_val <- max(vector_values, na.rm = TRUE)
  
  # Funció de normalització entre 0 i 1
  normalize_value <- function(x) {
    (x - min_val) / (max_val - min_val)
  }
  
  
  # Aplicar la normalització a cada element de la llista mantenint la seva estructura
  matrix_norm <- lapply(matrix_list, function(sublist) {
    # Normalitzar cada subllista o vector individualment
    sapply(sublist, normalize_value)
  })
  
  # Retornar la matriu normalitzada
  return(matrix_norm)
}

# Funció genèrica que retorna la divisió del conjunt de dades d'entrenament i de test: x_train, x_test, y_train, y_test
# Paràmetres: dataset, conjunt de dades
#             histological_grade, nom de la columna objectiu
# La divisió es realitza de forma arbitrària coma 2/3 per a train (70%) i 1/3 per a test (30%)
train_test_datasets_generation <- function(X, dataset, histological_grade) {
  #Divisió del conjunt de dades en subconjunts de train (entrenament, 70% per conveni 2/3) i test (proves, 30% 1/3)
  set.seed(123)  # Reproducibilitat

  # variable objectiu
  target <- dataset[[histological_grade]]
  
  # Generació dataset de train al 70% aleatori i amb una divisió estratificada: crea índexs per a train
  # per tal de poder tenir una proporció equilibrada d'ambdues classes (0,1) en els dos subconjunts de dades train i test.
  # La funció de 'caret' --> createDataPartition s'encarrega d'aquest punt
  train_indexs <- createDataPartition(y = target, p = 0.7, list = FALSE)
  
  # dades de train
  X_train <- X[train_indexs, ]
  y_train <- target[train_indexs]
  
  # dades de test
  X_test <- X[-train_indexs, ]
  y_test <- target[-train_indexs]
  
  cat("Dimensions del conjunt de dades d'entrenament: ", dim(X_train), " rows x cols \n")
  cat("Dimensions del conjunt de dades de test: ", dim(X_test), " rows x cols \n")
  
  cat("Distribució de classes en el conjunt de train:\n")
  print(table(y_train))
  
  cat("Distribució de classes en el conjunt de test:\n")
  print(table(y_test))

  return (list(X_train=X_train,
        X_test=X_test,
        y_train=y_train,
        y_test=y_test))

}

# Funció genèrica que crea un diagrama de barres horitzonals a partir de la llibreria ggplot2 per visualitzar
# el rànking de resultats de les importàncies dels gens en cada una de les tècniques de machine-learning diferents
# aplicades i per a un determinat grau histològic específic
# Paràmetres:
# - data: dades que contenen les importàncies de cada característica (gens) 
# - X: columna de dades en eix X
# - y: columna de dades en eix Y
# - color: color de les barres horitzontals
# - title: títol de la gràfica
# - top_num: nombre d'elements (gens) a visualitzar en les y's
# - xlablel: etiqueta en l'eix 'X'
# - ylabel: etiqueta en l'eix 'Y'
create_horitzontal_barchart_plot <- function(data, X, Y, color, title, top_num, xlabel, ylabel) {
  
  plot <- ggplot(data[1:top_num,], aes(x = reorder(.data[[X]], .data[[Y]]), y = .data[[Y]])) +
    geom_bar(stat = "identity", fill = color) +
    coord_flip() +
    labs(title = title, x = xlabel, y = ylabel) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
      plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
      panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
      panel.grid.minor = element_blank()                                # Línies de quadrícula menors
    )
  
  return(plot)
}


# Funció genèrica que crea un diagrama de punts a partir de la llibreria ggplot2 per visualitzar
# el rànking de resultats de les importàncies dels gens en cada una de les tècniques de machine-learning diferents
# aplicades i per a un determinat grau histològic específic
# Paràmetres:
# - data: dades que contenen les importàncies de cada característica (gens) 
# - X: columna de dades en eix X
# - y: columna de dades en eix Y
# - color: color de les barres horitzontals
# - title: títol de la gràfica
# - top_num: nombre d'elements (gens) a visualitzar en les y's
# - xlablel: etiqueta en l'eix 'X'
# - ylabel: etiqueta en l'eix 'Y'
create_point_chart_plot <- function(data, X, Y, color, title, top_num, xlabel, ylabel) {
  
  plot <- ggplot(data[1:top_num,], aes(x = .data[[Y]], y = reorder(.data[[X]], .data[[Y]]))) +
    geom_point(size = 3, color = color) +
    labs(title = title, x = xlabel, y = ylabel)
  theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
      plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
      panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
      panel.grid.minor = element_blank()                                # Línies de quadrícula menors
    )
  
  return(plot)
}


# Funció genèrica que crea un diagrama de barres horitzonals amb valors positius i/o negatius a partir de la llibreria ggplot2 
# per visualitzar el rànking de resultats de les importàncies dels gens en cada una de les tècniques de machine-learning 
# basades en la regularització i a partir d'un grau histològic específic
# Paràmetres:
# - data: dades que contenen les importàncies de cada característica (gens) 
# - X: columna de dades en eix X
# - y: columna de dades en eix Y
# - z: distinció del signe
# - color1: color de les barres horitzontals amb valors positius
# - color2: color de les barres horitzontals amb valors negatius
# - title: títol de la gràfica
# - subtitle: subtitol
# - top_num: nombre d'elements (gens) a visualitzar en les y's
# - xlablel: etiqueta en l'eix 'X'
# - ylabel: etiqueta en l'eix 'Y'
create_horitzontal_barchart_with_sign_plot <- function(data, X, Y, z, color1, color2, title, subtitle, top_num, xlabel, ylabel) {
  
  plot <- ggplot(data[1:top_num,], aes(x = reorder(.data[[X]], .data[[Y]]), y = .data[[Y]], fill= .data[[z]])) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Positiu (Risc)" = color1, "Negatiu (Protecció)" = color2)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlabel,
      y = ylabel
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
      plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
      panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
      panel.grid.minor = element_blank(),                               # Línies de quadrícula menors
      legend.position = "top"                                           # Llegenda a la part superior
    )
  
  return(plot)
}
  
# Funció genèrica que a partir d'una matriu de confusió com a paràmetre, retorna els valors de:
# - Sensibilitat
# - Especificitat
# - Exactitud
# - Precisió
calcula_metrics <- function(confusion_matrix) {
    
    # capturem
    # veritables positius: encerts positius
    TP <- confusion_matrix[2, 2]
    # falsos negatius: errors en etiquetar valors negatius que són positius
    FN <- confusion_matrix[2, 1]
    # falsos positius: errors en etiquetar valors positius que són negatius
    FP <- confusion_matrix[1, 2]
    # veritables negatius: encerts negatius
    TN <- confusion_matrix[1, 1]
    
    # fòrmules sensibilitat, especificitat, exactitiud i precisió
    
    # controlant indeterminació
    if (TP == 0 && FN == 0) {
      sensibilitat <- 0
    } else {
      sensibilitat <- TP / (TP + FN)
    }

    especificitat <- TN / (TN + FP)
    exactitud <- (TP + TN) / sum(confusion_matrix)
    precisio <- TP / (TP + FP)
    
    return(list(
      Sensibilitat = sensibilitat,
      Especificitat = especificitat,
      Exactitud = exactitud,
      Precisio = precisio
    ))
}

# Graficar corbes ROC per grau histològic
# -title (títol del gràfic)
# -data (conjunt de dades en format data frame que conté les especificitats i sensibilitats extretes per roc)
plot_roc <- function(data, title) {
  ggplot(data, aes(x = Specificity, y = Sensitivity)) +
    geom_line(color = "blue", size = 1) +
    geom_abline(linetype = "dashed", color = "gray") +
    labs(title = title, x = "1 - Especificitat", y = "Sensibilitat") +
  theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "white"), # Fons del gràfic
      plot.background = element_rect(fill = "white", color = "white"),  # Fons del quadre del gràfic
      panel.grid.major = element_line(color = "grey80"),                # Línies de la quadrícula
      panel.grid.minor = element_blank(),                               # Línies de quadrícula menors
      legend.position = "top"                                           # Llegenda a la part superior
    )
}


# Creació de la matriu de confusió amb inicialització per omplir valors buits
# Funció auxiliar per corregir el tema de no generar valors quan no hi ha prediccions per a una classe concreta
# Paràmetres: 
# - y_pred_labels (valors etiquetats resultants de les prediccions)
# - y_test (etiquetes reals originals del conjunt de dades de test)
# - classes (etiquetes possibles dins les categories, en aquest cas sempre serà binària de 0 / 1)
generate_confusion_matrix <- function(y_pred_labels, y_test, classes = c(0, 1)) {
  
  # Inicialitzar una matriu de confusió buida amb totes les combinacions possibles
  conf_matrix <- matrix(0, nrow = length(classes), ncol = length(classes),
                         dimnames = list(Predicted = classes, Actual = classes))
  
  # Generar la matriu de confusió amb les dades disponibles
  temp_conf_matrix <- table(Predicted = factor(y_pred_labels, levels = classes),
                      Actual = factor(y_test, levels = classes))
  
  # Omplir la matriu inicial amb els valors del `table()`
  conf_matrix[rownames(temp_conf_matrix), colnames(temp_conf_matrix)] <- temp_conf_matrix
  
  return(conf_matrix)
}
  
