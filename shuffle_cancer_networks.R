#------------------------------
# shuffle cancer networks to validate the results 
# got by the observed networks
#------------------------------

#-------- includes --------
library(readxl)
library(bipartite)
library(ggplot2)
library(tidyr)

#-------- consts --------
N_SIM = 1000 # number of shuffles to preform
SEED = 49

#------- functions --------
get_potentials <- function(shuffss, sim_id) {
  #create union of the cancer networks
  union_table <- shuffss[[1]][,,sim_id]
  union_table[,] <- 0 
  
  for (j in names(shuffss)) {
    net <- shuffss[[j]][,,sim_id]
    union_table <- union_table + net
  }
  # calculate the potential
  binari_union <- 1*(union_table>0)
  chap_potential <- rowSums(binari_union, na.rm = FALSE, dims = 1)
  
  return(chap_potential)
}

get_realized_niche <- function(shuffss, i, potential) {
  chap_names <- rownames(shuffss[[1]][,,1])
  
  # calc chaperon degree per cancer
  degree_table_chap <- as.data.frame(chap_names)
  
  # go over the cancers
  for (netname in sheet_names) {
    net <- shuffss[[netname]][,,i]
    
    # go over the proteins
    chap_degrees <- c()
    for (chap in chap_names) {
      # calc degree
      chap_degrees[chap] <- sum(net[chap, ])
    }
    
    # put each cancer in a new column
    curr_cancer_degrees <- data.frame(chap_degrees)
    new_col_names <- c(colnames(degree_table_chap), netname)
    
    degree_table_chap <- merge(degree_table_chap, curr_cancer_degrees,
                               by.x="chap_names", by.y="row.names", all.x=TRUE)
    
    colnames(degree_table_chap) <- new_col_names
  }
  
  # rearrange row names
  row.names(degree_table_chap) <- degree_table_chap$chap_names
  folding_percent <- degree_table_chap[,-1]
  
  # calculate percent of potential used
  for (chapp in chap_names) {
    folding_percent[chapp, ] <- folding_percent[chapp,]/potential[chapp]
  }
  
  return(folding_percent)
}

calculate_ev_nestedness <- function(B){
  # It is faster to calculate the ev for smaller matrices. Because the leading
  # ev of BB^T and B^TB is the same, we first check how to produce A.
  if (nrow(B)<ncol(B)){
    A <- B%*%t(B)
  } else {
    A <- t(B)%*%B
  }
  ev_max <- max(eigen(A, symmetric = T, only.values = T)$values) 
  # Not calculating eigenvectors speeds calculatoins remarkably.
  return(ev_max)
}

#-------- load the networks from an excel file --------
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

#-------- shuffle networks n times and save the results --------
shuffs <- list()
for (name in sheet_names) {
  x <- networks[[name]]

  null <- vegan::nullmodel(x, method = 'curveball')
  shuffled_matrices <- simulate(null, nsim = N_SIM, burnin = 5000, seed = SEED)
  
  shuffs[[name]] <- shuffled_matrices
}

realized_shuffs <- list()

for (i in 1:N_SIM) {
  pot <- get_potentials(shuffs, i)
  real <- get_realized_niche(shuffs, i, pot)
  realized_shuffs[[i]] <- real # saving the results
}

#-------- observed vs shuffled nestedness --------
# compare the nestedness of observed realized niche to that of the shuffled networks

# read the observed
obs_real <- read.csv("output/chap_folding_percent.csv", row.names = 1)
ev_obs <- calculate_ev_nestedness(as.matrix(obs_real))

# turn dfs to matrices
rows.cols <- dim(realized_shuffs[[1]])
sheets <- length(realized_shuffs)
shuf_real_mats <- array(unlist(realized_shuffs), dim = c(rows.cols, sheets))
colnames(shuf_real_mats) <- colnames(realized_shuffs[[1]])
row.names(shuf_real_mats) <- row.names(realized_shuffs[[1]])

# get eigen values for shuffled
ev_shuffled <- apply(shuf_real_mats, MARGIN = 3, FUN = calculate_ev_nestedness)

# Calculate p-value
p_value <- sum(ev_shuffled>ev_obs)/N_SIM

p <- ggplot(as.data.frame(ev_shuffled), aes(x = ev_shuffled)) +
  geom_histogram(colour = "darkblue", aes(fill = ..count..), bins = 30) +
  scale_x_continuous(name = "Simulated eigenvalue") +
  scale_y_continuous(name = "Count") +
  ggtitle("Histogram of simulated eigenvalue for chaperones realized niche", 
          subtitle = paste("observed eigenvalue:", round(ev_obs,3))) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
p
ggsave("output/nestedness_shuffled_realized_niche.pdf")
