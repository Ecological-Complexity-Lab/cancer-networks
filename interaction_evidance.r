# -------------- interaction_evidance.r -----------
# This script is meant to help validate the network found in the paper by 
# by comparing the interactions we found to the ones that exist in the STRING database.
# - definition:
# an interaction is found if it exist in at least one of the cancers


# includes ---------
library(readxl)
library(ggplot2)
library(reshape2)
library(tidyr)
library(tidyverse)
library(dplyr)

source("functions.r")

# read data ------
# prot metadata:
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
#prots_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]

# reading all the data from the STRING database site (network and aliases) 
pairs <- read.table("9606.protein.links.full.v11.5.txt", sep=" ", header=TRUE, 
                    stringsAsFactors=FALSE, quote="", fill=FALSE)
alians <- read.table("9606.protein.aliases.v11.5.txt", sep="\t", header=FALSE, 
                     stringsAsFactors=FALSE, quote="", fill=TRUE)
alians  <- alians %>% filter(V3 == "Ensembl_gene")
prots_meta <- prots_meta %>% left_join(alians[,1:2], by=c("ENSID" = "V2"))

# filter only the interactions between mito-prots
mito_pairs <- pairs[(pairs$protein2 %in% prots_meta$V1) & 
                    (pairs$protein1 %in% prots_meta$V1), ]

# make it use symbols
ch_id <- prots_meta %>% select(V1, Symbol)
mito_pairs_ch <- mito_pairs %>% 
                 left_join(ch_id, by=c("protein1" = "V1")) %>%
                 left_join(ch_id, by=c("protein2" = "V1")) %>%
                 rename(prot1=Symbol.x, prot2=Symbol.y) %>% 
                 select( prot1, prot2, everything(), -protein1, -protein2,)

write.csv(mito_pairs_ch, file = "output/data/STRING_prot_pairs.csv", row.names = FALSE)
mito_pairs_ch <- read.csv(file = "output/data/STRING_prot_pairs.csv")

# filter for this to be an the bipartite network
biprt <- mito_pairs_ch[xor((mito_pairs_ch$prot1 %in% chaps_meta$Symbol), 
                           (mito_pairs_ch$prot2 %in% chaps_meta$Symbol)), ]


write.csv(biprt, file = "output/data/STRING_bipartit.csv", row.names = FALSE)
biprt <- read.csv(file = "output/data/STRING_bipartit.csv")


# how many of our interactions is found in the STRING database

# TODO go over pairs and see - if they have an interaction in our system, does that interaction have evidance in the DB?


