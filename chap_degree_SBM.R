#------------------------------
# Use Stochastic Block Modeling to fit chaperons and cancers to different 
# groups. it will be based on chaperon degree per cancer.
#------------------------------

#------ includes ------
library(blockmodels)
library(tidyverse)
library(magrittr)
library(readxl)

#------------- load the networks from an excel file --------
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

#------------- find chaperon degree per cancer -----------------
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
row.names(degree_table_chap) <- degree_table_chap[,1]
degree_table_chap <- degree_table_chap[,-1]

#------------ activate SBM on degree table -------------
library(blockmodels)

# chaps  = rows    = plants
# cancer = columns = pollinators
sbm_model <- BM_poisson('LBM', as.matrix(degree_table_chap)) # The model on a network with discrete edge values

#LBM -> for bipartite networks
#SBM -> for adjacent networks

# sbm_model$plotting <- '' # Set plots of the convergence process to off
# sbm_model$explore_max <- 10 # Optional: set a maximum number of groups to search for
sbm_model$estimate()

best_run <- which.max(sbm_model$ICL)

# parse results
# Get the memberships of nodes into clusters
x <- sbm_model$memberships[[best_run]]
# Matrix with cluster affiliation. High values indicate that a node is in a cluster
x$plot() # Left plot bacteria by group, right plot metabolites by group
groups_chaps <- x$Z1 # Groups of chaps (rows chaps, columns are groups)
groups_cancers <- x$Z2 # Groups of cancers (rows cancers, columns are groups)
group_id_chaps <- apply(groups_chaps, 1, FUN = which.max) # The node belongs to a cluster to which it has highest probability
group_id_cancers <- apply(groups_cancers, 1, FUN = which.max) # The node belongs to a cluster to which it has highest probability

# Get the final result
groups_chaps <- tibble(node=rownames(degree_table_chap), group=group_id_chaps) %>% mutate(type='Chaperon')
groups_cancers <- tibble(node=colnames(degree_table_chap), group=group_id_cancers) %>% mutate(type='Cancer')
results <- bind_rows(groups_chaps, groups_cancers)

write.csv(results, file = "output/SBM_grouping.csv", row.names = FALSE)
