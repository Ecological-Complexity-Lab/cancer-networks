library(bipartite)

# Network for example - rows are plants, columns are pollinators
memmott1999

# Look at the edge data
qplot(memmott1999, color='red', binwidth=1)

dfun(web = memmott1999) # d' index of specialization

memmott1999_bin <- 1*(memmott1999>0)

library(blockmodels)
library(tidyverse)
library(magrittr)

sbm_model <- BM_bernoulli('LBM', memmott1999_bin) # The model on a binary network
sbm_model <- BM_poisson('LBM', memmott1999) # The model on a network with discrete edge values

# sbm_model$plotting <- '' # Set plots of the convergence process to off
# sbm_model$explore_max <- 10 # Optional: set a maximum number of groups to search for
sbm_model$estimate()

# parse results
# Get the memberships of nodes into clusters
x <- sbm_model$memberships[[best_run]]
# Matrix with cluster affiliation. High values indicate that a node is in a cluster
x$plot() # Left plot bacteria by group, right plot metabolites by group
groups_Plants <- x$Z1 # Groups of plants (rows plants, columns are groups)
groups_Pollin <- x$Z2 # Groups of pollinators (rows pllinators, columns are groups)
group_id_Plants <- apply(groups_Plants, 1, FUN = which.max) # The node belongs to a cluster to which it has highest probability
group_id_Pollin <- apply(groups_Pollin, 1, FUN = which.max) # The node belongs to a cluster to which it has highest probability

# Get the final result
groups_Plants <- tibble(node=rownames(memmott1999), group=group_id_Plants) %>% mutate(type='Plant')
groups_Pollin <- tibble(node=colnames(memmott1999), group=group_id_Pollin) %>% mutate(type='Pollinator')
results <- bind_rows(groups_Plants, groups_Pollin)

