## FUNCIONS GENÈRIQUES

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