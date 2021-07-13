#------------------------------
# build a table that shows the percent of the folding potential that is
# activated for each chaperonr for every cancer.
# 
#------------------------------


#------ includes ------
library(readxl)
library(ggplot2)
library(pheatmap)

#-------------- load the networks from an excel file --------
excel_path <- "HPC/binari_validated_corrs.xlsx"

# do the convert for every cancer
sheet_names <- excel_sheets(excel_path)

networks <- list()
for (name in sheet_names) {
  x <- as.data.frame(read_excel(excel_path, sheet = name, col_names = TRUE))
  x2 <- x[,-1]
  rownames(x2) <- x[,1]
  networks[[name]] <- x2
}

#----------- create union of the cancer networks -------------
union_table <- networks[[1]]
union_table[,] <- 0 

for (i in 1:length(networks)) {
  net <- networks[[i]]
  union_table <- union_table + net
}
union_table_mat <- as.matrix(union_table)

#----------- flat the union and calc folding potential -------
binari_union <- 1*(union_table>0)

chap_potential <- rowSums(binari_union, na.rm = FALSE, dims = 1)
chap_potential

#----------- calc chaperon degree per cancer -----------------
# get the list of proteins ENSID as rownames, SYBOL as first row
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
Symbol <- chaps_meta[, 4]
degree_table_chap <- as.data.frame(Symbol)

# go over the cancers
for (i in 1:length(networks)) {
  netname <- names(networks)[i]
  net <- networks[[i]]
  
  # go over the proteins
  chap_degrees <- c()
  for (chap in 1:nrow(chaps_meta)) {
    SYMB <- chaps_meta[chap, "Symbol"]
    # calc relative degree
    
    chap_degrees[SYMB] <- sum(net[SYMB, ])
  }
  
  # put each cancer in a new column
  curr_cancer_degrees <- data.frame(chap_degrees)
  new_col_names <- c(colnames(degree_table_chap), netname)
  
  degree_table_chap <- merge(degree_table_chap, curr_cancer_degrees,
                             by.x="Symbol", by.y="row.names", all.x=TRUE)
  degree_table_chap
  
  colnames(degree_table_chap) <- new_col_names
}


#----------- calculate percent of potential used ----------
row.names(degree_table_chap) <- degree_table_chap$Symbol
folding_percent <- degree_table_chap[,-1]

for (chapp in chaps_meta$Symbol) {
  chap_pot <- chap_potential[chapp]
  folding_percent[chapp, ] <- folding_percent[chapp,]/chap_pot
}

write.csv(folding_percent*100, file = "output/chap_folding_percent.csv")

pheatmap(folding_percent, cutree_rows = 2, cutree_cols = 2,
         filename = "output/folding_percentage.pdf",
         clustering_method = "ward.D", clustering_distance_rows = "manhattan",
         main = "Folding percentage by cancer")
