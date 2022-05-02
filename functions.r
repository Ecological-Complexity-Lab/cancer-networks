# ------------------------------
# a file comtianinf functions that repeat and are common 
# to multiple script files
# ------------------------------

# load the cancer network matrices -------
load_cancer_mats <- function(excel_path="HPC/binari_validated_corrs.xlsx") {
  # do the convert for every cancer
  sheet_names <- excel_sheets(excel_path)
  
  networks <- list()
  for (name in sheet_names) {
    x <- as.data.frame(read_excel(excel_path, sheet = name, col_names = TRUE))
    x2 <- x[,-1]
    rownames(x2) <- x[,1]
    networks[[name]] <- x2
  }
  
  return(networks)
}

# build a cancer union matrix -------
get_cancer_union <- function(networks) { 
  # assumption: all tables have the same structure
  union_table <- networks[[1]]
  union_table[,] <- 0 
  
  for (i in 1:length(networks)) {
    net <- networks[[i]]
    union_table <- union_table + net
  }

  return(union_table)
}



