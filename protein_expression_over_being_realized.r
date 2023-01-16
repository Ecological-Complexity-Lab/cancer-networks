#------------------------------
# This scripts finds the correlation between a protein's expression to his
# percent of being realized by a chaperon
#
#------------------------------

#-------- includes --------
library(readxl)
library(tidyr)

source("functions.r")

#-------- load data --------
# load expression data
prot_exp <- read.csv("output/prot_median_expressions.csv", row.names = 1)

# load cancer networks
networks <- load_cancer_mats()

#-------- test predictability --------
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chap <- "HSPD1"
prot <- "ABCD1"
# TODO get ESID for this prot 


net <- networks[["BRCA"]]

df <- prot_exp %>% select(expression = all_of(prot))

# TODO we may not need this at all. wait for the meeting.
