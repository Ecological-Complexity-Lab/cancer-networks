#! /gpfs0/shai/projects/R4/R-4.2.0/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.2.0/lib64/R/library")
print(.libPaths())
print(sessionInfo())

#------------ chap_stringdb_affirmation.r -------------
# this script calculates the the percentage of affirmation in randome proteins
# to later compare this to the observed percentage and see if it is enriched.

# required data:
# 1. full STRING database for humans
# 2. STRING alians data
# 3. local mln data
# 4. protein metadata

# includes ---------
library(readxl)
library(ggplot2)
library(reshape2)
library(tidyverse)

# functions -----
# get per chap - percentage of local that was found in experimental/ database
get_backed_percent_line <- function(evid_all, chaperone, to_compare) {
  relevant_entries <- evid_all %>% filter(Chap == chaperone) %>% 
    filter(evidence == "local" | evidence == to_compare)
  bbb <- relevant_entries %>% group_by(Chap, Prot) %>% summarise(n=n())
  
  chap_intr <- as.numeric(table(relevant_entries$evidence)["local"])
  supported <- as.numeric(table(bbb$n)["2"])
  
  chap_percnt <- 100*supported/chap_intr
  
  return(tibble(Chap=chaperone, evidence=to_compare,value=chap_percnt))
}

# get the number of random proteins that have evidence for a correlation - 
# when sampling random <obs_prots> proteins from all over the homo-sapiens data set
percentage_producer <- function(chap_pairs, obs_prots, all_prots) {
  # sample from all the genome X genes
  rand_prot <- sample(1:nrow(all_prots), obs_prots, FALSE)
  
  # check percentage of them that have experimental evidence - do this 10K times
  exp_count <- 0
  db_count <- 0
  for (id in rand_prot) {
    name <- all_prots[id, 1]
    uu <- chap_pairs[(chap_pairs$protein2 == name),]
    if (nrow(uu) > 0){
      exp_count <- exp_count + ifelse(uu$experiments>0, yes=1, no=0)
      db_count <- db_count + ifelse(uu$database>0, yes=1, no=0)
    }
  }
  return(c(exp_count, db_count))
}

# Handle args -----
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  chaperon <- as.character(args[1])
}
#chaperon <- "HSPD1"


# consts -------
n_rands <- 10000
str_db_file <- "9606.protein.links.full.v11.5.txt"
str_alias_file <- "9606.protein.aliases.v11.5.txt"
prots_file <- "Mito_genes.tab"
mln_file <- "adjacency_edgelist.dat"
output_file <- "STRINGdb_rand_percentage.csv"

# run algorithm -----
print("--- Read data used for the analysis --- ")
# read data:
mln <- read.delim(mln_file, sep = " ", header = FALSE) # our network
pairs <- read.table(str_db_file, sep=" ", header=TRUE, 
                    stringsAsFactors=FALSE, quote="", fill=FALSE)

# read metadata:
prots_meta <- read.table(prots_file, sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
alians <- read.table(str_alias_file, sep="\t", header=FALSE, 
                     stringsAsFactors=FALSE, quote="", fill=TRUE)
alians  <- alians %>% filter(V3 == "Ensembl_gene")
prots_meta <- prots_meta %>% left_join(alians[,1:2], by=c("ENSID" = "V2"))

# get chaperon name as a STRING ID
str_chap <- prots_meta %>% filter(Symbol == chaperon)
str_chap <- str_chap$V1

print("--- count potential interaction in the local co-expression mln --- ")
# count number of potential interactions for chaperon
prots <- unique(mln$V3)

coex_intr <- 0
for (prot in prots) {
  intr <- mln[(mln$V2 == chaperon) & (mln$V3 == prot),]
  if (sum(intr[,4:ncol(intr)]) > 0) { # if this interaction exists in our network
    coex_intr = coex_intr + 1
  } 
}

# get pairs of chap with all the proteins the human genome
pairs <- pairs[(pairs$protein1 == str_chap), ]

print("--- Randomly sample from the whole human genome and check affirmation percentage --- ")
# simulate random sampling of pairs for chaperon - to get a distribution
exp_list <- numeric(n_rands)
db_list <- numeric(n_rands)
for (i in 1:n_rands) {
  a <-  100*percentage_producer(pairs, coex_intr, alians)/coex_intr
  exp_list[i] <- a[1]
  db_list[i] <- a[2]
}

e <- data.frame(Chap = chaperon, evidence = "experiments", 
                value=exp_list, type="shuff")
d <- data.frame(Chap = chaperon, evidence = "databases", 
                value=db_list, type="shuff")
all <- rbind(e, d)

# save into a file
print("--- write results to a file ---")
write_csv(all, file = output_file, append = T)

