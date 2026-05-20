#' Load and Pre-process Single-Cell CEF Data
#'
#' This function reads a .cef file (Common Exchange Format), extracts basic metadata 
#' about the number of genes and cells, and returns a cleaned, transposed data frame 
#' ready for downstream analysis.
#'
#' @param file_path A string providing the path to the .cef.txt file.
#' @param verbose Logical. If TRUE (default), prints the gene and cell counts to the console.
#'
#' @return A data frame where rows represent cells and columns represent genes. 
#'         Numeric values are automatically converted back to numeric type.
#' @export
#'
#' @examples
#' # df <- process_cef_data("path/to/GSE76381_EmbryoMoleculeCounts.cef.txt")
process_cef_data <- function(file_path, verbose = TRUE) {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Original RMarkdown logic
  df_entry <- read.table(
    file_path,
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Optional metadata report before removing rows
  if (verbose) {
    n_genes <- suppressWarnings(as.numeric(df_entry[1, 5]))
    n_cells <- suppressWarnings(as.numeric(df_entry[1, 6]))
    
    message("--- Data Summary ---")
    message("Number of genes: ", n_genes)
    message("Number of cells: ", n_cells)
    message("--------------------")
  }
  
  df_entry <- df_entry[-c(1, 5), ]
  
  name <- df_entry[, 1]
  
  df <- as.data.frame(
    t(df_entry[, -1]),
    stringsAsFactors = FALSE
  )
  
  colnames(df) <- name
  
  return(df)
}