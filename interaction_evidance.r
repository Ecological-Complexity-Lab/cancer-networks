# -------------- interaction_evidance.r -----------
# This script is meant to help validate the network found in the paper by 
# by comparing the interactions we found to the ones that exist in the STRING database.
# - definition:
# an interaction is found if it exist in at least one of the cancers


# includes ---------
library(readxl)
library(ggplot2)
library(reshape2)
library(tidyverse)

source("functions.r")

# functions -----
update_counts <- function(line, evd) {
  for (var in names(line[,3:ncol(line)])) {
    if (line[var] > 0) {
      evd[var] = evd[var]+1
    }
  }
  return(evd)
}


# read data ------
# prot metadata:
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
clients_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]

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


# look at filtered data -----
# check how many of our interactions is found in the STRING database:

# read the two networks
biprt <- read.csv(file = "output/data/STRING_bipartit.csv") # string
mln <- read.delim("output/data/adjacency_edgelist.dat", sep = " ", header = FALSE) # our network

# go over pairs and see - if they have an interaction in our system, 
#                         does that interaction have evidence in the STRING db?
evidances <- integer(ncol(biprt)-2)
names(evidances) <- names(biprt[,3:ncol(biprt)])

coex_intr <- 0
na_pairs <- 0
for (chap in chaps_meta$Symbol) {
  print(paste(chap, "- start"))
  for (prot in clients_meta$Symbol) {
    intr <- mln[(mln$V2 == chap) & (mln$V3 == prot),]
    if (sum(intr[,4:ncol(intr)]) > 0) { # if this interaction exists in our network
      coex_intr = coex_intr + 1
      # find the relevant pair in the STRING db
      string_line <- biprt[(biprt$prot1 == chap) & (biprt$prot2 == prot),]
      
      if (nrow(string_line) == 0){ # if the interaction is not in STRING db
        na_pairs = na_pairs + 1
      } else {
        evidances <- update_counts(string_line, evidances)
      }
    } 
  }
  
}
print("All done.")

rrr <- melt(data.frame(as.list(evidances)))
rrr$percentage <- 100*rrr$value/coex_intr
write.csv(rrr, file = "output/data/backedup_interactions.csv", row.names = FALSE)

# plot the comparison results ------
rrr <- read.delim("output/data/backedup_interactions.csv", sep = ",", header = TRUE) # our network

ggplot(rrr, aes(x=variable, y=percentage)) + geom_point() +
  ylab("interactions with supporded evidence (%)") + paper_figs_theme_no_legend +
theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
      axis.title.x=element_blank()) 






# TODO tests to delete later
nrow(mln)

intr <- mln[(mln$V2 == "HSPD1") & (mln$V3 == "RMDN3"),]

biprt[(biprt$prot1 == "CYC1") & (biprt$prot2 == "CLPP"),]
biprt[(biprt$prot1 == "CLPP") & (biprt$prot2 == "CYC1"),]

sum(intr[,4:ncol(intr)])


