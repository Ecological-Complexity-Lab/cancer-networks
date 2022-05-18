#------------- stability_check.r ---------------------
# test the stability of each cancer network underportubation
# remove a chap and see how many die, in iterations.
#
#------------------------------

cancer_sample_size_order <- c("KIRP", # 288
                               "LIHC", # 371
                               "STAD", # 375
                               "COAD", # 478
                               "PRAD", # 498
                               "HNSC", # 500
                               "LUSC", # 502
                               "THCA", # 502
                               "LUAD", # 533
                               "KIRC", # 538
                               "UCEC", # 551
                               "BRCA") # 1102
cancer_percent_of_valid_interactions <- c("BRCA", # 22%
                                          "LUSC", # 31%
                                          "LUAD", # 33%
                                          "HNSC", # 48%
                                          "KIRC", # 49%
                                          "UCEC", # 49%
                                          "COAD", # 50%
                                          "LIHC", # 58%
                                          "PRAD", # 58%
                                          "STAD", # 61%
                                          "THCA", # 61%
                                          "KIRP") # 100%

#-------- includes --------
library(readxl)
library(bipartite)
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(reshape2)

source("functions.r")

#------ local functions ------
# taken from the stability game
# This function removes a single species from either the rows or the columns
single_extinct <- function(m, margin=NULL, x){
  if (margin==1){m[x,] <- 0}
  if (margin==2){m[,x] <- 0}
  return(m)
}

# This function counts the number of remaining species in either rows or columns
record_2nd_extinctions <- function(m, margin=NULL){
  # Notice that the margins are the opposite!!!
  # Because if removal is from rows we need to quantify columns and vice versa
  if (margin==2){x <- sum(rowSums(m)==0)}
  if (margin==1){x <- sum(colSums(m)==0)}
  return(x)
}

# This is the main function
extinct <- function(m, margin, extinction_sequence){
  # m is the matrix, margin is from where to remove (1 rows 2 columns)
  
  num_extinct <- 0 # Initialize the number of species going extinct
  # Loop through the extinction sequence
  for (e in extinction_sequence){
    # print(e)
    m <- single_extinct(m,margin,e)
    num_extinct <- c(num_extinct, record_2nd_extinctions(m, margin))
  }
  # Produce the results
  if (margin==1){
    results <- data.frame(num_removed=0:nrow(m), num_extinct=num_extinct)
    results$prop_removed <- results$num_removed/nrow(m)
    results$prop_remain <- 1-results$num_extinct/ncol(m)
  }
  if (margin==2){
    results <- data.frame(num_removed=0:ncol(m), num_extinct=num_extinct)
    results$prop_removed <- results$num_removed/ncol(m)
    results$prop_remain <- 1-results$num_extinct/nrow(m)
  }
  return(results)
}

extinction_per_matric <- function(mat, removal_order) {
  # preapare to a matrix
  A <- data.matrix(mat)
  A[is.na(A)] <- 0 # remove NA values
  
  df <- extinct(A, 1, extinction_sequence = removal_order) # Run extinction function on rows
  
  # Calculate area under the curve
  y <- df[,2]
  y <- (sum(y) - cumsum(y))/sum(y)
  x <- df$prop_removed
  # x <- (object[, "no"]/max(object[, "no"]))
  ext.curve <- splinefun(x, y)
  ext.area <- integrate(ext.curve, 0, 1)
  R <- as.numeric(ext.area[[1]])
  # Plot
  p <- ggplot(df, aes(prop_removed, prop_remain))+
    geom_point(size=2, color='red')+
    geom_line(size=1, color='red')+
    labs(title=paste('Your score: ',round(100*R,1)), x='Prop. removed', y='Prop. remained') +
    theme_bw() # Make a plot
  return(list(R=R, df=df, plt=p))
}


#-------- run --------
networks <- load_cancer_mats()

#-------- play for all cancer ----
all_R_values <- NULL
all_dfs <- NULL

# remove from most connected to least connected
for (cncr in names(networks)) {
  net <- data.matrix(networks[[cncr]])
  dead_order <- order(rowSums(net), decreasing = TRUE) # from most connected to least connected
  
  res <- extinction_per_matric(net, dead_order)
  res$df$cancer <- cncr
  res$df$run_type <- "high_to_low"
  tbl <- tibble(cancer=cncr, 
                under_curve=res$R, 
                run_type="high_to_low", 
                sim=0, 
                chap_removed=paste(rownames(net)[dead_order], collapse='->' ))
  all_R_values <- rbind(all_R_values, tbl)
  all_dfs <- rbind(all_dfs, res$df)
}

# remove from least connected to most connected
for (cncr in names(networks)) {
  net <- data.matrix(networks[[cncr]])
  dead_order <- order(rowSums(net)) # from least connected to most connected
  
  res <- extinction_per_matric(net, dead_order)
  res$df$cancer <- cncr
  res$df$run_type <- "low_to_high"
  tbl <- tibble(cancer=cncr, 
                under_curve=res$R, 
                run_type="low_to_high", 
                sim=0, 
                chap_removed=paste(rownames(net)[dead_order], collapse='->' ))
  all_R_values <- rbind(all_R_values, tbl)
  all_dfs <- rbind(all_dfs, res$df)
}

# remove chaperons by random and check
for (i in 1:500) {
  for (cncr in names(networks)) {
    net <- data.matrix(networks[[cncr]])
    dead_order <- sample(x = 1:nrow(net), size = nrow(net), replace = FALSE)# random
    
    res <- extinction_per_matric(net, dead_order)
    res$df$cancer <- cncr
    res$df$run_type <- paste("random_", i, sep="")
    tbl <- tibble(cancer=cncr, 
                  under_curve=res$R, 
                  run_type="random", 
                  sim=i, 
                  chap_removed=paste(rownames(net)[dead_order], collapse='->'))
    all_R_values <- rbind(all_R_values, tbl)
    all_dfs <- rbind(all_dfs, res$df)
  }
}

write_csv(all_R_values, "output/data/stability_results.csv")
write_csv(all_dfs, "output/data/stability_all_steps.csv")

p3 <- ggplot(all_R_values, aes(x=under_curve, fill=run_type)) + 
  geom_histogram(aes(y = stat(density)),
                 alpha=0.5, position = 'dodge',
                 bins=20) + 
  labs(x="area under extinction curve")
p3

all_R_values %>% filter(cancer=="BRCA") 

## between cancers
# x axis ordered by cancer sample size
pdf("output/stability.pdf")
ggplot(all_R_values, aes(x=fct_relevel(cancer, cancer_sample_size_order), 
                         y=under_curve)) + 
  geom_boxplot() + ggtitle("stability. chaps according to sample size")


# x axis ordered by pecent of interactions remain valid in the normalization
ggplot(all_R_values, aes(x=fct_relevel(cancer, cancer_percent_of_valid_interactions), 
                         y=under_curve)) + 
  geom_boxplot() + ggtitle("stability. chaps according to valid %")

# between types of extinctions
ggplot(all_R_values, aes(x=run_type, 
                         y=under_curve)) + 
  geom_boxplot() + ggtitle("stability. by removal type")
dev.off()
#-> no significant change between types, 
#   maybe because of the large difference between chaperones and proteins and the way they are connected
