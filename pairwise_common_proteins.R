#------------------------------
# build a table that shows the number of proteins that are common for each
# couple of chaperons. in each cancer, and for the union of the cancers
#------------------------------


#------ includes ------
library(readxl)
library(ggplot2)
library(tibble)
library(tidyverse)
library(reshape2)
source("functions.r")

#------ functions ------
get_intersection_from_mat <- function(nett) { # network MUST be binary!
  # in the function with: nett
  chaps <- rownames(nett)
  chap_mat <- matrix(0, nrow = nrow(nett), ncol = nrow(nett), 
                     dimnames = list(chaps, chaps))
  
  for (chap1 in 1:(length(chaps)-1)) {
    for (chap2 in (chap1+1):length(chaps)) {
      two <- nett[c(chap1,chap2),]
      common <- ncol(two[,colSums(two)>1])
      # if there is only one column a list is returned, so the ncol() returns a NULL
      common <- ifelse(is.null(common), 1, common)
      
      chap_mat[chap2, chap1] <- common
    }
  }
  chap_mat[upper.tri(chap_mat, diag = TRUE)] <- NA
  
  return(chap_mat)
}

reorder_chapchap <- function(chpchp){
  # make it as a full mat
  diag(chpchp) <- 1
  chpchp[upper.tri(chpchp)]  <- t(chpchp)[upper.tri(chpchp)]
  
  # Use intersecton percentage between variables as distance
  dd <- as.dist((1-chpchp)/2)
  hc <- hclust(dd)
  chpchp <- chpchp[hc$order, hc$order]
  chpchp[upper.tri(chpchp, diag = TRUE)] <- NA
  return(chpchp)
}

#----------- load the networks from an excel file --------
networks <- load_cancer_mats()

#----------- create union of the cancer networks -------------
union_table <- get_cancer_union(networks)

# turn it to binari
binari_union <- 1*(union_table>0)

#----------- create intersection networks -------------
# get intersections per cancer
intersections <- list()
intersections_df <- NULL
for (cncr in names(networks)) {
  nett <- networks[[cncr]]
  inter_mat <- get_intersection_from_mat(nett)
  inter_ordered <- reorder_chapchap(inter_mat/n_proteins)
  intersections[[cncr]] <- inter_ordered
  
  inter_long <- melt(inter_ordered, na.rm = TRUE)
  inter_long$cancer <- cncr
  intersections_df <- rbind(intersections_df, inter_long)
}
write.csv(intersections_df, file = "output/data/chapchap_intersections.csv", row.names = FALSE)

# get intersections for the union of the cancers
union_chpchap_mat <- get_intersection_from_mat(binari_union)

#----------- show matrix as tile plot -------------
n_proteins <- ncol(binari_union)

# union plot
union_chpchap_mat_ordered <- reorder_chapchap(union_chpchap_mat/n_proteins)
union_chpchap <- melt(union_chpchap_mat_ordered, na.rm = TRUE)
write.csv(union_chpchap, file = "output/data/chapchap_intersections_union.csv", row.names = FALSE)

p <- ggplot(data = union_chpchap, aes(Var2, Var1, fill = value))+
  geom_tile(color = "Black")+
  scale_fill_gradient(low = "cornsilk", high = "blue4", 
                      limit = c(0,max(union_chpchap$value)),
                      name="Intersection") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) + coord_fixed()
ggsave("output/paper_figures/chap_intersection_cancer_union.pdf",
       plot = p, width = 5)
       

# each cancers
p2 <- ggplot(data = intersections_df, aes(Var2, Var1, fill = value))+
  geom_tile(color = "Black")+
  scale_fill_gradient(low = "cornsilk", high = "blue4", 
                      limit = c(0,max(intersections_df$value)),
                      name="Intersection") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  coord_fixed() + facet_wrap(~ cancer)
ggsave("output/chap_intersection_per_cancer.pdf", plot = p2)
