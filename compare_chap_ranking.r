#------------------------------
# this script compares the order of two lists
# it is used to compare the order of chaperones between:
# 1. realised niche - by row sum, high to low
# 2. log10 expression level - by row sum, high to low
# 3. folding percentage out of all proteins - by row sum, high to low
# 4. folding potential - high to low
#------------------------------

#-------- includes --------
library(tidyverse)


#-------- assign lists --------
realized_niche <- c("CLPP", "HSPE1", "YME1L1","CLPX", "HSPD1",
                    "HSPA9", "HSCB", "AFG3L2", "DNAJC19", "LONP1",
                    "DNAJA3", "TRAP1", "HTRA2", "GRPEL2","SPG7")
expression <- c("HSPD1", "HSPA9", "HSPE1", "LONP1", "CLPP",
                "YME1L1", "TRAP1", "DNAJA3","AFG3L2", "CLPX",
                "HSCB", "DNAJC19", "HTRA2", "SPG7", "GRPEL2")
foldind_percent <- c("HSPE1","CLPP", "HSPD1", "DNAJC19", "DNAJA3", 
                     "AFG3L2", "HSPA9", "HSCB", "TRAP1", "HTRA2",
                     "YME1L1", "LONP1", "CLPX", "GRPEL2", "SPG7")
folding_potential <- c("DNAJC19", "DNAJA3", "TRAP1", "HTRA2", "HSPD1", 
                       "HSPE1", "AFG3L2", "HSPA9", "GRPEL2", "LONP1", 
                       "SPG7", "HSCB", "CLPP", "CLPX", "YME1L1")

#-------- rank the lists --------
niche_df <- data.frame(Chaperone = realized_niche, 
                       Rank = 1:length(realized_niche)) 

percent_df <- data.frame(Chaperone = foldind_percent) %>% 
           left_join(niche_df,by = "Chaperone")

exprsn_df <- data.frame(Chaperone = expression) %>% 
           left_join(niche_df,by = "Chaperone")

potential_df <- data.frame(Chaperone = folding_potential) %>% 
  left_join(niche_df,by = "Chaperone")


#-------- rank the lists --------
# use Kendall rank correlation test to see how similar is their order
niche_percent_res<-cor.test(niche_df$Rank, percent_df$Rank, method="kendall")
niche_percent_res

niche_pot_res<-cor.test(niche_df$Rank, potential_df$Rank, method="kendall")
niche_pot_res

niche_exp_res<-cor.test(niche_df$Rank, exprsn_df$Rank, method="kendall")
niche_exp_res
