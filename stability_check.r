#------------- stability_check.r ---------------------
# test the stability of each cancer network underportubation
# remove a chap and see how many die, in iterations.
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

# returns the order of chaperons to remove from the matrix,
# it first orders by module number and then descending by degree.
order_by_module <- function(nett, chap_attr, module_order = c(3,2,1)) {
  chap_degs <- as.data.frame(rowSums(nett))
  chap_info <- merge(chap_attr, t(t(chap_degs)), by=0) %>% 
    select(Row.names, module, degree = "rowSums(nett)") %>% 
    #arrange(module, desc(degree))
    arrange(match(module, module_order), desc(degree))
  index_order <- match(chap_info$Row.names, row.names(nett))
  return(index_order)
}

# removing for all the cancers, by a specific module order
# uses the function above
remove_by_module <- function(networks, chap_attrib, m_order=c(1,2,3),
                             all_R_values, all_dfs) {
  order_num <- m_order[1]*100 + m_order[2]*10 + m_order[3]
  # remove from most connected to least connected - by module!
  for (cncr in names(networks)) {
    net <- data.matrix(networks[[cncr]])
    net <- net[,colSums(net) > 0] # remove proteins without edges
    
    dead_order <- order_by_module(net, chap_attrib,module_order = m_order)
    
    res <- extinction_per_matric(net, dead_order)
    res$df$cancer <- cncr
    res$df$run_type <- "by_module"
    res$df$sim <- order_num
    tbl <- tibble(cancer=cncr, 
                  under_curve=res$R, 
                  run_type="by_module", 
                  sim=order_num, 
                  chap_removed=paste(rownames(net)[dead_order], collapse='->' ))
    all_R_values <- rbind(all_R_values, tbl)
    all_dfs <- rbind(all_dfs, res$df)
  }
  return(list(Rs=all_R_values, process=all_dfs))
}


#-------- load --------
networks <- load_cancer_mats()

#-------- play for all cancer ----
all_R_values <- NULL
all_dfs <- NULL

# remove from most connected to least connected
for (cncr in names(networks)) {
  net <- data.matrix(networks[[cncr]])
  net <- net[,colSums(net) > 0] # remove proteins without edges
  
  dead_order <- order(rowSums(net), decreasing = TRUE) # from most connected to least connected
  
  res <- extinction_per_matric(net, dead_order)
  res$df$cancer <- cncr
  res$df$run_type <- "high_to_low"
  res$df$sim <- 0
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
  net <- net[,colSums(net) > 0]
  
  dead_order <- order(rowSums(net)) # from least connected to most connected
  
  res <- extinction_per_matric(net, dead_order)
  res$df$cancer <- cncr
  res$df$run_type <- "low_to_high"
  res$df$sim <- 0
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
    net <- net[,colSums(net) > 0]
    
    dead_order <- sample(x = 1:nrow(net), size = nrow(net), replace = FALSE)# random
    
    res <- extinction_per_matric(net, dead_order)
    res$df$cancer <- cncr
    res$df$run_type <- "random"
    res$df$sim <- i
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


# ------ Visualize -------
all_R_values <- read.csv("output/data/stability_results.csv")
all_dfs <- read.csv("output/data/stability_all_steps.csv")

p3 <- ggplot(all_R_values, aes(x=under_curve, fill=run_type)) + 
  geom_histogram(aes(y = stat(density)),
                 alpha=0.5, position = 'dodge',
                 bins=20) + 
  labs(x="area under extinction curve")
p3


## between cancers
# x axis ordered by cancer sample size
pdf("output/stability.pdf")
ggplot(all_R_values, aes(x=fct_relevel(cancer, cancer_sample_size_order), 
                         y=under_curve, color=run_type)) + 
  geom_boxplot() + ggtitle("stability, chaps according to sample size") +
  labs(x=element_blank(), y="area under extinction curve")

# x axis ordered by pecent of interactions remain valid in the normalization
ggplot(all_R_values, aes(x=fct_relevel(cancer, cancer_percent_of_valid_interactions), 
                         y=under_curve, color=run_type)) + 
  geom_boxplot() + ggtitle("stability, chaps according to valid %")+
  labs(x=element_blank(), y="area under extinction curve")

# between types of extinctions
ggplot(all_R_values, aes(x=run_type, 
                         y=under_curve, 
                         color=run_type)) + 
  geom_boxplot() + labs(x=element_blank(), y="area under extinction curve") +
  facet_wrap(vars(cancer), ncol = 4) + paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

# prepare the random so there's only one value per step
df <- all_dfs %>% filter(run_type!="random") %>% 
  select(num_removed, cancer, prop_remain, prop_removed, run_type)
collapsed_random <- all_dfs %>% filter(run_type=="random") %>%
  group_by(num_removed, cancer) %>%
  summarise(mean_left=mean(prop_remain), mean_removed=mean(prop_removed), run_type="random") %>%
  select(num_removed, cancer, prop_remain=mean_left, prop_removed=mean_removed, run_type)

all_types <- rbind(df, collapsed_random)
# plot the collapse while it is happening per cancer
ggplot(all_types, aes(prop_removed, prop_remain, color=run_type))+
  geom_point(size=2)+
  geom_line(size=1)+
  labs(x="% of chaperons removed", 
       y="% of proteins remained",
       color="Removal type")+
  facet_wrap(vars(cancer), ncol = 4) +
  paper_figs_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1))

dev.off()
#-> is there a significant change between types? 
#   if not, maybe it's because of the large difference between the number 
#   of chaperons and proteins and the way they are connected


# ------ play for all cancers - by module -------
# read chap attributes (module numbers)
chap_attrib <- read.csv("output/data/chap_attributes.csv", row.names = 1) %>% 
               select(module, amount)
by_mdl_R_values <- NULL
by_mdl_dfs <- NULL

# generate all different module orders
combs <- expand.grid(1:3, 1:3, 1:3) %>% 
         filter((Var1!=Var2) & (Var2!=Var3) & (Var3!=Var1))

for (i in 1:nrow(combs)) {
  ord <- as.numeric(combs[i,])
  res <- remove_by_module(networks, chap_attrib, ord,
                          by_mdl_R_values, by_mdl_dfs)
  by_mdl_R_values <- res$Rs
  by_mdl_dfs      <- res$process
}

by_mdl_R_values$run_name <- paste(by_mdl_R_values$run_type,
                                  by_mdl_R_values$sim, sep = "_")
by_mdl_dfs$run_name      <- paste(by_mdl_dfs$run_type,
                                  by_mdl_dfs$sim, sep = "_")

# save results
write_csv(by_mdl_R_values, "output/data/stability_results_by_module.csv")
write_csv(by_mdl_dfs, "output/data/stability_all_steps_by_module.csv")

# ------make data ready for visualization ------
# read for visualization
by_mdl_R_values <- read.csv("output/data/stability_results_by_module.csv")
by_mdl_dfs <- read.csv("output/data/stability_all_steps_by_module.csv")

# adding random and high_to_low
all_R_values <- read.csv("output/data/stability_results.csv") %>%
  mutate(run_name=run_type) %>% 
  filter((run_type=="high_to_low")|(run_type=="random")|(run_type=="low_to_high"))
all_dfs <- read.csv("output/data/stability_all_steps.csv") %>%
  mutate(run_name=run_type) %>% 
  filter((run_type=="high_to_low")|(run_type=="random")|(run_type=="low_to_high"))

# both
both_R_vals <- rbind(by_mdl_R_values, all_R_values)
both_dfs <- rbind(by_mdl_dfs, all_dfs)

# check correlate -------
# for r X nestedness colSum correlation
median_rndm <- both_R_vals %>% select(cancer, under_curve, run_name) %>%
  filter(run_name=="random") %>% group_by(cancer) %>%
  summarise(under_curve=median(under_curve), run_name="median_random") 
# prepare the random so there's only one value per step
df_no_rand <- both_R_vals %>% filter(run_name!="random") %>% 
  select(cancer, under_curve, run_name)

med_rand <- rbind(df_no_rand, median_rndm)

# get nestedness data
rn <- read.csv("output/chap_realized_niche.csv", row.names = 1)
rn_sums <- as.data.frame(colSums(rn)) %>% tibble::rownames_to_column("cancer")
cancr_sums <- med_rand %>% 
      left_join(rn_sums, by="cancer") %>% rename(rn_sum = `colSums(rn)`)

# calculate correlation co-efficient and P value
corrs <- cancr_sums %>%
  group_by(run_name) %>%
  summarize(cor=cor(under_curve, rn_sum, method = "spearman"), 
            p=cor.test(under_curve,rn_sum)$p.value,
            estimate=cor.test(under_curve,rn_sum)$estimate)

write.csv(corrs, "output/data/stability_r_vs_rn_col_sums_corrs.csv")

# for r X nestedness col MEAN correlation
# prepare nestedness data
rn <- read.csv("output/chap_realized_niche.csv", row.names = 1)
rn_means <- as.data.frame(colMeans(rn)) %>% tibble::rownames_to_column("cancer")
cancr_means <- med_rand %>% 
  left_join(rn_means, by="cancer") %>% rename(rn_mean = `colMeans(rn)`)

# calculate correlation co-efficient and P value
corrs2 <- cancr_means %>%
  group_by(run_name) %>%
  summarize(cor=cor(under_curve, rn_mean, method = "spearman"), 
            p=cor.test(under_curve,rn_mean)$p.value,
            estimate=cor.test(under_curve,rn_mean)$estimate)

write.csv(corrs2, "output/data/stability_r_vs_rn_col_mean_corrs.csv")



# -------- Visualize - removal by module -------
pdf("output/stability_by_module.pdf", 12, 6)
# between types of extinctions
ggplot(both_R_vals, aes(x=run_name, 
                            y=under_curve, 
                            color=run_name)) + 
  geom_boxplot() + labs(x=element_blank(), y="area under extinction curve") +
  facet_wrap(vars(cancer), ncol = 4) + paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))


# by removal method
both_R_vals %>% select(cancer, under_curve, run_name) %>%
  filter(run_name != 'by_module_321') %>% 
  filter(run_name != 'low_to_high') %>% 
  filter(run_name != 'by_module_312') %>% 
  filter(run_name != 'by_module_132') %>% 
  filter(run_name != 'by_module_231') %>% 
  mutate(cancer=factor(cancer, levels=cancer_nestedness_order)) %>%
  ggplot(aes(x=cancer, y=under_curve, color=run_name)) + 
  geom_boxplot() + 
  labs(x=element_blank(), y="area under extinction curve") +
  facet_wrap(vars(run_name), ncol = 4) + paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

# correlate R with nestedness colSum
cancr_sums %>% 
  filter(run_name != 'by_module_321') %>% 
  filter(run_name != 'low_to_high') %>% 
  filter(run_name != 'by_module_312') %>% 
  filter(run_name != 'by_module_132') %>% 
  filter(run_name != 'by_module_231') %>% 
  ggplot(aes(x=rn_sum, y=under_curve, color=run_name)) + 
  geom_point() +
#  stat_cor(label.y = 0.85, method = "spearman", 
#           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_smooth(method = "lm", se=FALSE) +
  labs(x="Cancer realized niche sum", y="area under extinction curve") +
  facet_wrap(vars(run_name), ncol = 4) + paper_figs_theme

cancr_means %>% 
  filter(run_name != 'by_module_321') %>% 
  filter(run_name != 'low_to_high') %>% 
  filter(run_name != 'by_module_312') %>% 
  filter(run_name != 'by_module_132') %>% 
  filter(run_name != 'by_module_231') %>% 
  ggplot(aes(x=rn_mean, y=under_curve, color=run_name)) + 
  geom_point() +
  #  stat_cor(label.y = 0.85, method = "spearman", 
  #           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_smooth(method = "lm", se=FALSE) +
  labs(x="Cancer realized niche mean", y="area under extinction curve") +
  facet_wrap(vars(run_name), ncol = 4) + paper_figs_theme
 
# prepare the random so there's only one value per step
df <- both_dfs %>% filter(run_name!="random") %>% 
  select(num_removed, cancer, prop_remain, prop_removed, run_name)
collapsed_random <- both_dfs %>% filter(run_name=="random") %>%
  group_by(num_removed, cancer) %>%
  summarise(mean_left=mean(prop_remain), mean_removed=mean(prop_removed), run_name="random") %>%
  select(num_removed, cancer, prop_remain=mean_left, prop_removed=mean_removed, run_name)

all_types <- rbind(df, collapsed_random)

# prepare the random so there's only one value per cancer
df <- both_R_vals %>% filter(run_name!="random") %>% 
  select(cancer, under_curve, run_name)
collapsed_random <- both_R_vals %>% filter(run_name=="random") %>%
                    group_by(cancer, run_name) %>%
                    summarise(mean_R=mean(under_curve)) %>%
                    select(cancer, under_curve=mean_R, run_name)

slim_rs <- rbind(df, collapsed_random)

# plot the collapse while it is happening per cancer
all_types %>% 
  filter(run_name != 'by_module_321') %>% 
  filter(run_name != 'low_to_high') %>% 
  filter(run_name != 'by_module_312') %>% 
  filter(run_name != 'by_module_132') %>% 
  filter(run_name != 'by_module_231') %>% 
  mutate(cancer=factor(cancer, levels=cancer_nestedness_order)) %>% 
  left_join(slim_rs, by=c("cancer"="cancer", "run_name"="run_name")) %>%
  mutate(y_txt= case_when(run_name=='by_module_123' ~ 0.04,
                          run_name=='by_module_213' ~ 0.16,
                          run_name=='high_to_low' ~ 0.28, 
                          run_name=='random' ~ 0.40))%>%
ggplot(aes(prop_removed, prop_remain, color=run_name))+
  geom_point(size=2)+
  geom_line(size=1)+
  labs(x="% of chaperons removed", 
       y="% of proteins remained",
       color="Removal type") +
  # Add the R to each cancer and each removal plot in the white-space using: 
  geom_text(aes(x=0.07, y=y_txt, label=round(under_curve, digits = 3)), stat = "unique") +

  facet_wrap(vars(cancer), ncol = 4) +
  paper_figs_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1))

dev.off()


# ------ play for all cancers - by functionality -------
# read chap attributes (module numbers)
chap_attrib <- read.csv("output/data/chap_attributes.csv", row.names = 1) %>% 
  select(module, amount, function., complax, is_co_chap)

both <- list(Rs=NULL, steps=NULL)

# generate order by functionality
res_from_ext <- function(nett, index_order, run_title, bothh) {
  res <- extinction_per_matric(nett, index_order)
  res$df$cancer <- cncr
  res$df$run_type <- run_title
  res$df$sim <- 0
  tbl <- tibble(cancer=cncr, 
                under_curve=res$R, 
                run_type=run_title, 
                sim=0, 
                chap_removed=paste(rownames(nett)[index_order], collapse='->'))
  bothh$Rs <- rbind(bothh$Rs, tbl)
  bothh$steps <- rbind(bothh$steps, res$df)
  
  return(bothh)
}

# run for all the functionalist
for (cncr in names(networks)) {
  net <- data.matrix(networks[[cncr]])
  net <- net[,colSums(net) > 0]

  # make each removal order seperatly
  chap_degs <- as.data.frame(rowSums(net))
  chap_info <- merge(chap_attrib, t(t(chap_degs)), by=0) %>% 
    select(Row.names, module, degree = "rowSums(net)", complax, function., is_co_chap)
  
  ord1 <- chap_info %>% arrange(is_co_chap, desc(degree)) # first chaps -> co-chaps
  ord2 <- chap_info %>% arrange(complax, desc(degree)) # by complax (then degree)
  ord3 <- chap_info %>% arrange(function., desc(degree)) # folding -> protease
  ord4 <- chap_info %>% arrange(desc(function.), desc(degree)) # folding <- protease
  
  index_order1 <- match(ord1$Row.names, row.names(net))
  index_order2 <- match(ord2$Row.names, row.names(net))
  index_order3 <- match(ord3$Row.names, row.names(net))
  index_order4 <- match(ord4$Row.names, row.names(net))
  
  both <- res_from_ext(net, index_order1, "chap_to_cochap", both)
  both <- res_from_ext(net, index_order2, "by_complax", both)
  both <- res_from_ext(net, index_order3, "fold_to_cutting", both)
  both <- res_from_ext(net, index_order4, "cutting_to_fold", both)
}

both$Rs$run_name <- both$Rs$run_type
both$steps$run_name <- both$steps$run_type


# check the correlation of this with rn sum -----

med_rand <- both$Rs 

# for r X nestedness col MEAN correlation
# prepare nestedness data
rn <- read.csv("output/chap_realized_niche.csv", row.names = 1)
rn_means <- as.data.frame(colMeans(rn)) %>% tibble::rownames_to_column("cancer")
cancr_means <- med_rand %>% 
  left_join(rn_means, by="cancer") %>% rename(rn_mean = `colMeans(rn)`)

# calculate correlation co-efficient and P value
corrs3 <- cancr_means %>%
  group_by(run_name) %>%
  summarize(cor=cor(under_curve, rn_mean, method = "spearman"), 
            p=cor.test(under_curve,rn_mean)$p.value,
            estimate=cor.test(under_curve,rn_mean)$estimate)

write.csv(corrs3, "output/data/stability_r_vs_rn_col_mean_corrs_by_func.csv")



# visualize the results - by protein functionality ------
pdf("output/stability_by_func.pdf", 12, 6)
# between types of extinctions
ggplot(both$Rs, aes(x=run_name, 
                        y=under_curve, 
                        color=run_name)) + 
  geom_boxplot() + labs(x=element_blank(), y="area under extinction curve") +
  facet_wrap(vars(cancer), ncol = 4) + paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))


# by removal method
both$Rs %>% select(cancer, under_curve, run_name) %>%
  mutate(cancer=factor(cancer, levels=cancer_nestedness_order)) %>%
  ggplot(aes(x=cancer, y=under_curve, color=run_name)) + 
  geom_boxplot() + 
  labs(x=element_blank(), y="area under extinction curve") +
  facet_wrap(vars(run_name), ncol = 4) + paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

# plot the collapse while it is happening per cancer
both$steps %>% 
  mutate(cancer=factor(cancer, levels=cancer_nestedness_order)) %>% 
  left_join(both$Rs, by=c("cancer"="cancer", "run_name"="run_name")) %>%
  mutate(y_txt= case_when(run_name=='chap_to_cochap' ~ 0.04,
                          run_name=='by_complax' ~ 0.16,
                          run_name=='fold_to_cutting' ~ 0.28, 
                          run_name=='cutting_to_fold' ~ 0.40))%>%
  ggplot(aes(prop_removed, prop_remain, color=run_name))+
  geom_point(size=2)+
  geom_line(size=1)+
  labs(x="% of chaperons removed", 
       y="% of proteins remained",
       color="Removal type") +
  # Add the R to each cancer and each removal plot in the white-space using: 
  geom_text(aes(x=0.07, y=y_txt, label=round(under_curve, digits = 3)), stat = "unique") +
  
  facet_wrap(vars(cancer), ncol = 4) +
  paper_figs_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1))

cancr_means %>% 
  ggplot(aes(x=rn_mean, y=under_curve, color=run_name)) + 
  geom_point() +
  #  stat_cor(label.y = 0.85, method = "spearman", 
  #           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_smooth(method = "lm", se=FALSE) +
  labs(x="Cancer realized niche mean", y="area under extinction curve") +
  facet_wrap(vars(run_name), ncol = 4) + paper_figs_theme

dev.off()


# ------ compare R to evenness and ------- 

# get info
networks_clean <- lapply(networks, function(x) x[,colSums(x) > 0])
df <- data.frame(cancer=names(networks_clean))
df$mean_deg <- unlist(lapply(networks_clean, function(x) mean(colSums(x))))
df$diver <- unlist(lapply(networks_clean, 
                          function(x) vegan::diversity(colSums(x), 
                                                       index = 'shannon')))
df$deliminator <- unlist(lapply(networks_clean, 
                          function(x) log(length(colSums(x)))))
df <- df %>% mutate(evenness = diver/deliminator) %>% 
             mutate(evenness_mean = mean_deg*(diver/deliminator)) %>% 
             select(cancer, mean_deg, evenness, evenness_mean)

# is it okay to compare when the number of "species" is not the same across cancers?
Rs <- both_R_vals %>% filter(run_type!="random") %>% 
  select(cancer, under_curve, run_name) %>%
  pivot_wider(names_from = run_name, values_from = under_curve)

all_columns_to_compare <- df %>% left_join(Rs, by="cancer")



# testing correlations:
library(corrplot)
library(Hmisc)

# corr mean_deg
cor(all_columns_to_compare[,5:ncol(all_columns_to_compare)], 
    all_columns_to_compare$mean_deg)
cor.test(all_columns_to_compare[,5:ncol(all_columns_to_compare)], 
         all_columns_to_compare$mean_deg)

# corr evenness
cor(all_columns_to_compare[,5:ncol(all_columns_to_compare)], 
    all_columns_to_compare$evenness)

# corr evenness_mean
cor(all_columns_to_compare[,5:ncol(all_columns_to_compare)], 
    all_columns_to_compare$evenness_mean)
 
pdf("output/test_evenness.pdf")
cor_5 <- rcorr(as.matrix(all_columns_to_compare[,2:ncol(all_columns_to_compare)]))
M <- cor_5$r
p_mat <- cor_5$P
corrplot(M, type = "upper", order = "hclust", 
         p.mat = p_mat, sig.level = 0.01)
dev.off()
 