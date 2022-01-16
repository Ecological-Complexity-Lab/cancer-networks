#------------------------------
# shuffle cancer networks to validate the results 
# got by the observed networks
# Here "generalism" is what i thought to call "folding percentage" - 
# the percent out of all the proteins.
#------------------------------

#-------- includes --------
library(readxl)
library(bipartite)
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(reshape2)

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

get_generalism <- function(shuffss, i) {
  chap_names <- rownames(shuffss[[1]][,,1])
  prot_names <- colnames(shuffss[[1]][,,1])
  
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
  
  # calculate percentage of all proteins used
  percentage <- degree_table_chap[,-1]/length(prot_names)
  
  return(percentage) 
}

get_similarity_per_cancer <- function(shuffss, i) {
  all_simlrs <- matrix(0, nrow = 105, ncol = 0)
  
  # go over the cancers
  for (netname in names(shuffss)) {
    net <- shuffss[[netname]][,,i]
    
    # generate similarity between the chaperons in this cancer
    res <- vegdist(net, method="jaccard")
    simlr <- 1-res
    
    all_simlrs <- cbind(all_simlrs, simlr)
  }
  
  colnames(all_simlrs) <- names(shuffss)

  return(all_simlrs)
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

#-------- shuffle networks n times and save the results - curveball --------
shuffs <- list()
for (name in sheet_names) {
  x <- networks[[name]]

  null <- vegan::nullmodel(x, method = 'curveball')
  shuffled_matrices <- simulate(null, nsim = N_SIM, burnin = 5000, seed = SEED)
  
  shuffs[[name]] <- shuffled_matrices
}

realized_shuffs <- list()
cancr_jsccard <- list()

for (i in 1:N_SIM) {
  pot <- get_potentials(shuffs, i)
  real <- get_realized_niche(shuffs, i, pot)
  cncr_simlr <- get_similarity_per_cancer(shuffs, i)
  
  realized_shuffs[[i]] <- real # saving the results
  cancr_jsccard[[i]] <- cncr_simlr 
}

#-------- shuffle networks n times and save the results - r00 --------
shuffs_r00 <- list()
for (name in sheet_names) {
  x <- networks[[name]]
  
  null <- vegan::nullmodel(x, method = 'r00')
  shuffled_matrices <- simulate(null, nsim = N_SIM, burnin = 5000, seed = SEED)
  
  shuffs_r00[[name]] <- shuffled_matrices
}

generalism_shuffs <- list()
for (i in 1:N_SIM) {
  of_all <- get_generalism(shuffs_r00, i)
  generalism_shuffs[[i]] <- of_all # saving the results
}

#-------- observed vs shuffled nestedness of realized niche --------
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


#-------- observed vs shuffled nestedness of generalism --------
# compare the nestedness of observed generalism to that of the shuffled networks

# read the observed
obs_of_all <- read.csv("output/chap_folding_percent_of_all.csv", row.names = 1)
ev_all_obs <- calculate_ev_nestedness(as.matrix(obs_of_all))

# turn dfs to matrices
rows.cols <- dim(generalism_shuffs[[1]])
sheets <- length(generalism_shuffs)
shuf_gen_mats <- array(unlist(generalism_shuffs), dim = c(rows.cols, sheets))
colnames(shuf_gen_mats) <- colnames(generalism_shuffs[[1]])
row.names(shuf_gen_mats) <- row.names(generalism_shuffs[[1]])

# get eigen values for shuffled
ev_shuffled <- apply(shuf_gen_mats, MARGIN = 3, FUN = calculate_ev_nestedness)

# Calculate p-value
p_value <- sum(ev_shuffled>ev_all_obs)/N_SIM

p <- ggplot(as.data.frame(ev_shuffled), aes(x = ev_shuffled)) +
  geom_histogram(colour = "darkblue", aes(fill = ..count..), bins = 30) +
  scale_x_continuous(name = "Simulated eigenvalue") +
  scale_y_continuous(name = "Count") +
  ggtitle("Histogram of simulated eigenvalue for chaperones' generalism", 
          subtitle = paste("observed eigenvalue:", round(ev_all_obs,3))) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
p
ggsave("output/nestedness_shuffled_generalism.pdf")

# when using curveball this is all the same because the shuffling doesn't change the degrees. 
# this is why i used r00 for this analysis

#-------- observed vs shuffled similarity distribution --------
cancer_simlr <- read.csv("output/jaccard_values_per_cancer.csv")
mlt_obs <- as.data.frame(melt(cancer_simlr))

# turn dfs to matrices
rows.cols <- dim(cancr_jsccard[[1]])
sheets <- length(cancr_jsccard)
shuf_sim_mats <- array(unlist(cancr_jsccard), dim = c(rows.cols, sheets))

mlt_sim <- as.data.frame(melt(shuf_sim_mats))


p1<-ggplot(mlt_sim, aes(x=value)) + 
    geom_histogram(aes(y = stat(count) / sum(count)), 
                   bins=30, color="black") + 
    labs(title = "shuffled networks", x="Jaccard similarity index")
p1

p2<-ggplot(mlt_obs, aes(x=value)) + 
    geom_histogram(aes(y = stat(count) / sum(count)), 
                 bins=30, color="black") + 
    labs(title = "Observed network", x="Jaccard similarity index")
p2


combine_dfs <- mlt_sim %>% select(kind = Var1, value)
temp <- mlt_obs %>% select(kind = variable, value)
combine_dfs["kind"] <- "shuff"
temp["kind"] <- "obs"
combine_dfs <- rbind(temp, combine_dfs)

p3 <- ggplot(combine_dfs%>% group_by(kind), aes(x=value, fill=kind)) + 
  geom_histogram(aes(y = stat(density)),
                 alpha=0.5, position = 'identity',
                 bins=30) + 
  labs(title = "Jaccard distribution in observed vs shuffled networks", 
       x="Jaccard similarity index")
p3
ggsave("output/nestedness_shuffled_jaccard.pdf")
