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
# update for protein pair the what string data we have
update_counts <- function(line, evd) { 
  for (var in names(line[,3:ncol(line)])) {
    if (line[var] > 0) {
      evd[var] = evd[var]+1
    }
  }
  return(evd)
}

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




# read metadata ------
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
clients_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]

# read and process STRING data ------
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


# look at filtered data - get total evidence coverage from STRING db ----- 
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


# high resolution evidence data -------
#read data:
biprt <- read.csv(file = "output/data/STRING_bipartit.csv") # string
mln <- read.delim("output/data/adjacency_edgelist.dat", sep = " ", header = FALSE) # our network

# go over pairs and build the evidence tibble
evid_local <- NULL
evid_string <- NULL
evid_genma <- NULL

na_pairs <- 0
for (chap in chaps_meta$Symbol) {
  print(paste(chap, "- start"))
  for (prot in clients_meta$Symbol) {
    intr <- mln[(mln$V2 == chap) & (mln$V3 == prot),]
    if (sum(intr[,4:ncol(intr)]) > 0) { # if this interaction exists in our network
      row <- tibble(Chap=chap, Prot=prot, DB="local", evidence="local", weight=1)
      evid_local <- rbind(evid_local, row)
      
      # find the relevant pair in the STRING db
      string_line <- biprt[(biprt$prot1 == chap) & (biprt$prot2 == prot),]
      
      if (nrow(string_line) == 0){ # if the interaction is not in STRING db
        na_pairs = na_pairs + 1
      } else {
        # create a row in the tibble for each category
        for (var in names(string_line)[3:ncol(string_line)]) {
          val <- as.numeric(string_line[var])
          if (val > 0) {
            row <- tibble(Chap=chap, Prot=prot, DB="STRING", evidence=var, weight=val)
            evid_string <- rbind(evid_string, row)
          }
        }
      }
    } 
  }
}
print("All done.")


evid_all <- rbind(evid_local, evid_string, evid_genma)
write.csv(evid_all, file = "output/data/pairs_evidance.csv", row.names = FALSE)

# plot amount of interactions in general by category ------
# read data:
evid_all <- as_tibble(read.csv(file = "output/data/pairs_evidance.csv"))

#plot:
to_plot <- evid_all %>% group_by(Chap, DB, evidence) %>% summarise(n=n()) 

# remove transferred, only look at known for humans
to_plot <- dplyr::filter(to_plot, !grepl('_transferred', evidence))

ggplot(to_plot, aes(x=Chap, y=n)) + 
  geom_point() + facet_wrap(~evidence) +
  ylab("number of interactions") + paper_figs_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.title.x=element_blank()) 

ggplot(to_plot, aes(x=evidence, y=n)) + 
  geom_point() + facet_wrap(~Chap) +
  ylab("number of interactions") + paper_figs_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.title.x=element_blank()) 


# get *per chap* - percentage of local that was found in experimental/ database
perc_tibble <- NULL
for (chanl in c("experiments", "database")) {
  for (chapp in chaps_meta$Symbol) {
    row <- get_backed_percent_line(evid_all, chapp, chanl)
    perc_tibble <- rbind(perc_tibble, row)
  }
}

write.csv(perc_tibble, file = "output/data/STRING_affirm_percentage.csv", row.names = FALSE)


# plot percent of local interactions that was affirmed by STRING
ggplot(perc_tibble, aes(x=Chap, y=value, color=evidence)) + 
  geom_point(size=3) +
  ylab("% of interactions affirmed") + paper_figs_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.title.x=element_blank()) 

# Validation enrichment analysis ------------
# process - compare the union of each chap with experimental, to random 
# *** this process was done over the HPC ***

# Plot STRING db validation results -------
# check the observed percentage compared to the distribution - is it enriched?
# load observed:
perc_tibble <- as_tibble(read.csv(file = "output/data/STRING_affirm_percentage.csv"))
perc_tibble$type <- "obs"

# read random data
perc <- read.csv(file = "HPC/DATA/STRINGdb_rand_percentage.csv", header = F) # string
names(perc) <- c("Chap", "evidence", "value", "type")
exp_perc <- perc %>% filter(evidence == "experiments")
db_perc <- perc %>% filter(evidence == "databases")

# plot random distribution per chap + observed as a point.
all_exp <- rbind(perc_tibble %>% filter(evidence == "experiments"), exp_perc)
all_db <- rbind(perc_tibble %>% filter(evidence == "database"), db_perc)

ggplot() +
  geom_boxplot(data=all_exp, aes(x=factor(Chap), y=value)) + 
  geom_point(data=all_exp[all_exp$type == "obs",], 
             aes(x=factor(Chap), y=value), color="red", size = 3) + 
  paper_figs_theme + ylab("% of interactions affirmed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.title.x=element_blank())

ggplot() +
  geom_boxplot(data=all_db, aes(x=factor(Chap), y=value)) + 
  geom_point(data=all_db[all_db$type == "obs",], 
             aes(x=factor(Chap), y=value, color=factor(type))) + 
  paper_figs_theme_no_legend + ylab("% of interactions affirmed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.title.x=element_blank())


# abstractly check validation by list from specific papers: ----
# TODO
# build a relevant functions

# fetch the data from HSPD1 paper

# fetch the data from TRAP1 paper

# fetch the data from CLPP paper







