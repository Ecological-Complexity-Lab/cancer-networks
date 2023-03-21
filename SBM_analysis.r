# -------------- SBM_analysis.r -----------
#
#
#

# includes ---------
library(readxl)
library(vegan)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(igraph)

source("functions.r")

# consts ------
layer_order <- c("BRCA", "COAD", "HNSC", "KIRC", "KIRP", 
                 "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC") # abc order

# Runs -----------
# make the data to be in the format compatible with the MultiTansor tool -----
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

# run the MultiTensor tool ------------------
# save the format in the tool's data folder
write_delim(formatted_edglist, col_names = FALSE,
            file = "../../MultiTensor/data/adjacency_cancer.dat", delim = " ")

# run once -----
setwd("../../MultiTensor/python/")
call <- "python2 main.py -l=12 -k=3 -a=\"adjacency_cancer.dat\""
system(call)
setwd("../../git_root/cancer_neworks/")

# run for a range of community numbers ----
#check likelihood across runs
setwd("../../MultiTensor/python/")
for (i in 1:12) {
  
  call <- paste("python2 main.py -l=12 -k=", i, " -u=1 -a=\"adjacency_cancer.dat\" ", sep = "")
  print(call)
  system(call)
}
setwd("../../git_root/cancer_neworks/")

