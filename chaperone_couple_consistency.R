#------------------------------
# chaperones' couple therapy
# ron over all the chaperones that have a common protein, generate a table that 
# shows how consistant each chaperone-chaperone-client triangle is.
#------------------------------

#-------- includes --------
library(readxl)
library(tidyverse)

#-------- functions --------
calc_obs_couple_percent <- function(nets, chap1, chap2, prot) {
  count_cncr <- 0
  # make this to be as above
  for (name in names(networks)) {
    net <- networks[[name]]
    both <- net[chp1,prot] + net[chp2,prot] > 1
    count_cncr <- count_cncr+both
  }
  return(count_cncr/length(networks))
}

calc_sims_couple_percent <- function(nets, chap1, chap2, prot) {
  # TODO add randomizing to validate the results
  return(1)
}


#-------- prepare protein meta data --------
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
prots_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]


#-------- load the networks from an excel file --------
excel_path <- "HPC/binari_validated_corrs.xlsx"

# do the convert for every cancer
sheet_names <- excel_sheets(excel_path)

networks <- list()
for (name in sheet_names) {
  x <- as.data.frame(read_excel(excel_path, sheet = name, col_names = TRUE))
  x2 <- x[,-1]
  rownames(x2) <- x[,1]
  networks[[name]] <- x2
}


#-------- build the table --------
couple_df <- data.frame(chaperone1=character(),
                        chaperone2=character(),
                        protein=character(),
                        prot_name=character(),
                        obs_percent=double(),
                        random_percent=double(),
                        stringsAsFactors=FALSE)

for (i in 1:(nrow(chaps_meta)-1)){
  chp1 <- chaps_meta$Symbol[i]
  
  for (j in (i+1):nrow(chaps_meta)) {
    chp2 <- chaps_meta$Symbol[j]
    
    for (k in 1:nrow(prots_meta)) {
      prt <- prots_meta$ENSID[k]
      prt_s <- prots_meta$Symbol[k]
      obs <- calc_obs_couple_percent(networks, chp1, chp2, prt)
      if (obs>0) {
        sim <- calc_sims_couple_percent(networks, chp1, chp2, prt)
        
        couple_df %>% add_row(chaperone1=chp1, 
                              chaperone2=chp2,
                              protein=prt,
                              prot_name=prt_s,
                              obs_percent=obs,
                              random_percent=sim) -> couple_df
      }
    }
  }
}

write.csv(couple_df, file = "output/couple_folding_persistancy.csv")

