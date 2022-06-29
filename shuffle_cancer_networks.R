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

source("functions.r")

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

get_similarity_per_chapron <- function(shuffss, i) {
  a_net <- shuffss[[1]][,,1]
  prots <- colnames(a_net)
  chaps <- rownames(a_net)
  
  # build cancerXprot matrices
  # build empty df
  empty_net <- data.frame(matrix(nrow=length(sheet_names), 
                                 ncol = length(prots)))
  colnames(empty_net) <- prots
  row.names(empty_net) <- sheet_names
  
  chap_nets <- list()
  for (chap in chaps) {
    chap_net <- empty_net
    
    for (cancr in sheet_names) {
      cncr_df <- shuffss[[cancr]][,,i]
      
      chap_net[cancr,] <- cncr_df[chap,]
    }
    
    # add to network list
    chap_nets[[chap]] <- chap_net
  }
  
  all_simlrs <- matrix(0, nrow = 66, ncol = 0)
  for (chapp in names(chap_nets)) {
    net <- chap_nets[[chapp]]
    res <- vegdist(net, method="jaccard")
    simlr <- 1-res
    
    all_simlrs <- cbind(all_simlrs, simlr)
  }
  
  colnames(all_simlrs) <- names(chap_nets)
  return(all_simlrs)
}

calculate_ev_nestedness_old <- function(B){
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

calculate_ev_nestedness <- function(B){
  m <- nrow(B)
  n <- ncol(B)
  A <- matrix(0,m+n,m+n)
  A[(m+1):(n+m),1:m] <- t(B)
  A[1:m,(m+1):(n+m)] <- B
  stopifnot(isSymmetric(A))
  ev_max <- max(eigen(A, symmetric = T, only.values = T)$values) # Not calculating eigenvectors speeds calculations remarkably.
  return(ev_max)
}

#-------- load the networks from an excel file --------
networks <- load_cancer_mats()

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
chap_jsccard <- list()

for (i in 1:N_SIM) {
  pot <- get_potentials(shuffs, i)
  real <- get_realized_niche(shuffs, i, pot)
  cncr_simlr <- get_similarity_per_cancer(shuffs, i)
  chap_simlr <- get_similarity_per_chapron(shuffs, i)
  
  realized_shuffs[[i]] <- real # saving the results
  cancr_jsccard[[i]] <- cncr_simlr
  chap_jsccard[[i]] <- chap_simlr
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

#-------- saving the shuffling results --------
save(generalism_shuffs, realized_shuffs, cancr_jsccard, chap_jsccard, 
     file = "output/data/shuff_results.RData")

# load the r objects as is.
load("output/data/shuff_results.RData")

#-------- observed vs shuffled nestedness of realized niche --------
# compare the nestedness of observed realized niche to that of the shuffled networks

# read the observed
obs_real <- read.csv("output/chap_realized_niche.csv", row.names = 1)
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

write_csv(as.data.frame(c(ev_obs, ev_shuffled)), 
          "output/data/rn_vs_shuff_eigenvalues.csv")

pp <- ggplot(as.data.frame(ev_shuffled), aes(x = ev_shuffled)) +
  geom_histogram(colour = "darkblue", aes(fill = ..count..), bins = 30) +
  scale_x_continuous(name = "Simulated eigenvalue") +
  scale_y_continuous(name = "Count") +
  ggtitle("Histogram of simulated eigenvalue\nfor chaperones realized niche") +
  geom_vline(xintercept = ev_obs, linetype="dashed", color = "red", size=1) + 
  theme(plot.title = element_text(hjust = 0.5))
pp
ggsave("output/paper_figures/nestedness_shuffled_realized_niche.pdf",
       plot = pp, width = 3.2, height = 2.8, units = "in")


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

write_csv(as.data.frame(c(ev_all_obs, ev_shuffled)), 
          "output/data/generalizm_vs_shuff_eigenvalues.csv")

p <- ggplot(as.data.frame(ev_shuffled), aes(x = ev_shuffled)) +
  geom_histogram(colour = "darkblue", aes(fill = ..count..), bins = 30) +
  scale_x_continuous(name = "Simulated eigenvalue") +
  scale_y_continuous(name = "Count") +
  ggtitle("Histogram of simulated eigenvalue\nfor chaperones' generalism") +
  geom_vline(xintercept = ev_all_obs, linetype="dashed", color = "red", size=1) + 
  theme(plot.title = element_text(hjust = 0.5))
p
ggsave("output/paper_figures/nestedness_shuffled_generalism.pdf", 
       plot = p, width = 3.2, height = 2.8, units = "in")

# when using curveball this is all the same because the shuffling doesn't change the degrees. 
# this is why i used r00 for this analysis

#-------- observed vs shuffled similarity distribution - per cancer --------
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
combine_dfs["kind"] <- "Shuffled"
temp["kind"] <- "Observed"
combine_dfs <- rbind(temp, combine_dfs)

write_csv(combine_dfs, "output/data/cancer_jaccard_with_shuff.csv")

p3 <- ggplot(combine_dfs%>% group_by(kind), aes(x=value, fill=kind)) + 
        geom_histogram(aes(y = stat(density)),
                       alpha=0.5, position = 'identity',
                       bins=30) + 
        labs(x="Jaccard similarity index",
             y="Density", fill="Population") +
        paper_figs_theme
p3
ggsave("output/paper_figures/shuffled_jaccard_per_cancer.pdf", p3)


#-------- observed vs shuffled similarity distribution - per chap --------
chp_simlr <- read.csv("output/jaccard_values_per_chap.csv")
mlt_obs <- as.data.frame(melt(chp_simlr))

# turn dfs to matrices
rows.cols <- dim(chap_jsccard[[1]])
sheets <- length(chap_jsccard)
shuf_sim_mats <- array(unlist(chap_jsccard), dim = c(rows.cols, sheets))

mlt_sim <- as.data.frame(melt(shuf_sim_mats))


p4<-ggplot(mlt_obs, aes(x=value)) + 
  geom_histogram(#aes(y = stat(count) / sum(count)),
                 bins=30, color="black") + 
  labs(title = "Observed network", x="Jaccard similarity index")
p4

p5<-ggplot(mlt_sim, aes(x=value)) + 
  geom_histogram(aes(y = stat(count) / sum(count)), 
                 bins=30, color="black") + 
  labs(title = "Shuffled networks", x="Jaccard similarity index")
p5

combine_dfs <- mlt_sim %>% select(kind = Var1, value)
temp <- mlt_obs %>% select(kind = variable, value)
combine_dfs["kind"] <- "shuff"
temp["kind"] <- "obs"
combine_dfs <- rbind(temp, combine_dfs)

write_csv(combine_dfs, "output/data/chap_jaccard_with_shuff.csv")

p3 <- ggplot(combine_dfs%>% group_by(kind), aes(x=value, fill=kind)) + 
        geom_histogram(aes(y = stat(density)),
                       alpha=0.5, position = 'identity',
                       bins=30) + 
        labs(x="Jaccard similarity index",
             y="Density", fill="Population") +
        paper_figs_theme
p3
ggsave("output/paper_figures/shuffled_jaccard_per_chap.pdf")


# calculate distribution median for the figure caption
combine_dfs <- read_csv("output/data/chap_jaccard_with_shuff.csv")

combine_dfs %>% group_by(kind) %>%
  summarise(min=min(value), median=median(value), mean=mean(value), max=max(value))


# ---- z-score testing for the difference - per chap ----
# here we are testing for each chap whether it is
# more similar to itself then found by random.
chp_simlr <- read.csv("output/jaccard_values_per_chap.csv")

# need a shuffed mat of: (Jaccard_shuff)
# chap - mean_shuff_jaccard - sd_suff_jaccard - shuff_id
Jaccard_shuff <- melt(chap_jsccard) %>% select(chap=Var2, jaccard=value, run=L1)
Jaccard_shuff <- Jaccard_shuff %>% group_by(chap, run) %>%
                      summarise(run_mean=mean(jaccard), run_sd=sd(jaccard))

# and need a mat of: (Jaccard_obs)
# chap - mean_jaccard - sd_jaccard
Jaccard_obs <- melt(chp_simlr) %>% 
               select(chap=variable, jaccard=value) %>%
               group_by(chap) %>%
               summarise(obs_mean=mean(jaccard), obs_sd=sd(jaccard))

# calculate z-core and significance
PF_J_z_score <- 
  Jaccard_shuff %>%
  group_by(chap) %>%
  summarise(shuff_mean=mean(run_mean), shuff_sd=sd(run_mean)) %>% 
  inner_join(Jaccard_obs) %>%
  mutate(z=(obs_mean-shuff_mean)/shuff_sd) %>% 
  mutate(signif=case_when(z>1.96 ~ 'above', # Obs is more than the shuffled
                          z< -1.96 ~ 'below', # Obs is lower than the shuffled
                          z<=1.96 | z>=-1.96 ~ 'not signif'))

# What proportion of ASVs have a statistical significant PF_J?
PF_J_z_score %>% 
  group_by(signif) %>% 
  summarise(n=n(),prop=n/nrow(PF_J_z_score))

write.csv(PF_J_z_score, "output/data/chap_similarity_z_score.csv")

#-- same but ׳with doing a single step ----- 
# doing a summarise for the jaccard value across all runs and chaperones in one step
# and not doing two steps like the last section 
# (summarise per run, then per chap)
Jaccard_single_shuff <- melt(chap_jsccard) %>% select(chap=Var2, jaccard=value, run=L1)

PF_J_z_score <- 
  Jaccard_single_shuff %>%
  group_by(chap) %>%
  summarise(shuff_mean=mean(jaccard), shuff_sd=sd(jaccard)) %>% 
  inner_join(Jaccard_obs) %>%
  mutate(z=(obs_mean-shuff_mean)/shuff_sd) %>% 
  mutate(signif=case_when(z>1.96 ~ 'above', # Obs is more than the shuffled
                          z< -1.96 ~ 'below', # Obs is lower than the shuffled
                          z<=1.96 | z>=-1.96 ~ 'not signif'))

# ---- z-score testing for the distributon difference - per cancer ----
# check if all the values in all the cancers are in general lower then 
# all the values in all the cancers on the shuffled
cncr_simlr <- read.csv("output/jaccard_values_per_cancer.csv")
mlt_obs <- as.data.frame(melt(cncr_simlr))
Jaccard_cancr_shuff <- melt(cancr_jsccard) %>% select(cancer=Var2, jaccard=value, run=L1)



obs_mean <- mean(mlt_obs$value)
shuf_mean <- mean(Jaccard_cancr_shuff$jaccard)
shuf_sd <- sd(Jaccard_cancr_shuff$jaccard)

z <- (obs_mean-shuf_mean)/shuf_sd

result <- case_when(z>1.96 ~ 'above', # Obs is more than the shuffled
                    z< -1.96 ~ 'below', # Obs is lower than the shuffled
                    z<=1.96 | z>=-1.96 ~ 'not signif')

