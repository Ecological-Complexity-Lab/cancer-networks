# -------------- SBM_analysis.r -----------
# Prepares the format used to run SBM analysis, and later process the results.
# this analysis - SBM and link prediction was done by Xei and is present in the 
# relevant folders of the current repo - including raw output.

# includes ---------
library(readxl)
library(ggplot2)
library(reshape2)
library(tidyverse)

source("functions.r")

# consts ------
layer_order <- c("BRCA", "COAD", "HNSC", "KIRC", "KIRP", 
                 "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC") # abc order

# Runs -----------
# make the data to be in the format compatible with the MultiTansor tool -----
# this data was given to Xei, to run the SBM and link prediction.

# lists
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
prots_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]


# read the multilayer data
networks <- load_cancer_mats()

# convert to the right format:
# E node_name1 node_name2 <value_layer_1> <value_layer_2> ... <value_layer_12>
formatted_edglist <- NULL

for (chap in chaps_meta$Symbol) {
  for (prot in 1:nrow(prots_meta)) {
    prottt <- prots_meta[prot, ]
    row <- tibble(E="E", node1=chap, node2=prottt$Symbol)
    for (cncr in layer_order) {
      val <- networks[[cncr]][chap,prottt$ENSID]
      row[cncr] = val
    }
    formatted_edglist <- rbind(formatted_edglist, row)
  }
}

# save the formatted the multilayer adjacency matrix (as edgelist)
write_delim(formatted_edglist, col_names = FALSE,
            file = "output/data/adjacency_edgelist.dat", delim = " ")
# save to be used in other stuff then multitensor as well:\
write.csv(formatted_edglist[,2:ncol(formatted_edglist)], 
          file = "output/data/adjacency_edgelist.csv", row.names = FALSE)

# plot Xei's results ----------------
bi_file <- "output/data/bipartite_membership.csv"

membership_vercors <- read.table(bi_file, sep=",", header=TRUE,
                                 stringsAsFactors=FALSE, quote="", fill=FALSE)


long_formt <- melt(membership_vercors)

ggplot(long_formt, aes(x=value, X, fill=1))+geom_tile()+facet_wrap(~variable)


ggplot(long_formt, aes(variable, X, fill= as.factor(value))) + 
  geom_tile() +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())


# we only care about the first 2 community runs - so plot them nicely
long_formt2 <- melt(membership_vercors[, 1:4])

ggplot(long_formt2, aes(x=reorder(X,value),y=value, fill=1))+
  geom_tile()+
  facet_grid(variable~.)+ 
  ylab("Modul ID") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = 'none')

# Plot Xei's membership probabilities ----
# ploting only the probabilities for 2 communities
mem_prob <- "data_for_lp&cd/community_breast_colon_head_kidneyc_kidneyp_liver_lunga_lungs_prostate_stomach_thyroid_uterine__u_K2.dat"
membs <- read.table(mem_prob, sep=" ", header=FALSE, skip = 1,
                                 stringsAsFactors=FALSE, quote="", fill=FALSE)

membs <- membs[1:15,]
names(membs) <- c("chap", "1", "2")
# The order pf chaperons needs to be confirmed, it is taken from Xei's code without proper explanation
membs$Name <- c('HSPA9', 'HSPD1', 'HSPE1', 'SPG7', 'CLPP', 'DNAJA3', 'LONP1',
                     'TRAP1', 'YME1L1', 'CLPX', 'AFG3L2', 'HTRA2', 'DNAJC19', 'GRPEL2',
                     'HSCB') 

pop <- membs %>% select(chap, Name, "1", "2")

# save to be used in plotting
write.csv(pop, file = "output/data/chap_membership_matrix.csv", row.names = FALSE)

## plot Xei's link prediction results: -----
# string conversion - 
converter <- tibble(code=layer_order,
             title=c("breast",
                     "colon",
                     "head",
                     "kidneyc",
                     "kidneyp",
                     "liver",
                     "lunga",
                     "lungs",
                     "prostate",
                     "stomach",
                     "thyroid",
                     "uterine"))

#read link prediction data:
lp_file <- "code_for_link_prediction&community_detection/link_prediction.csv"
lp_data <- read.table(lp_file, sep=",", header=TRUE,
                      stringsAsFactors=FALSE, quote="", fill=FALSE)
# prepare data format
diagnl <- lp_data[1:12,] %>% add_column(second=lp_data$X[1:12]) %>% 
  select(first=X, second, AUC)
non_dig <- lp_data[13:nrow(lp_data),] %>% separate(X, c('first','second'))
back_togather <- rbind(diagnl, non_dig)

to_save <- back_togather %>% 
  left_join(converter, by=c("first"="title")) %>% 
  left_join(converter, by=c("second"="title")) %>%
  select(From=code.x, To=code.y, AUC, everything())

write.csv(to_save, file = "output/data/link_prediction_plot_ready.csv", row.names = F)

bb <- ggplot(back_togather, aes(x=from, y=to, fill=AUC)) + 
    geom_tile() + 
    scale_fill_gradient(low = "lightyellow", high = "navyblue") +
    paper_figs_theme + 
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          panel.border = element_blank())
bb

# calculate the jaccard between cancer monolayers: -----
# read data:
aj_file <- "output/data/adjacency_edgelist.csv"
aj_data <- read.csv(aj_file, header=TRUE,
                      stringsAsFactors=FALSE, fill=FALSE)

jaccard_results <- NULL
for (cncr1 in layer_order) {
  for (cncr2 in layer_order) {
    only_them <- aj_data[,c(cncr1, cncr2)]
    only_them$sum <- rowSums(only_them)
    both <- nrow(only_them %>% filter(sum == 2))
    unio <- nrow(only_them %>% filter(sum > 0))
    
    jaccard_results <- rbind(jaccard_results,
                             tibble(cancer1=cncr1, 
                                    cancer2=cncr2, 
                                    jaccard=both/unio))
  }  
}

# save to be used in plotting
write.csv(jaccard_results, 
          file = "output/data/cancer_edges_jaccard.csv", row.names = FALSE)
