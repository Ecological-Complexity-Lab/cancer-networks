#------------------------------
# use similarity indexes to find the differences 
# by comparing between the chaperones in each cancer
#------------------------------

#-------- includes --------
library(readxl)
library(vegan)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)

#-------- functions --------
long_format_from_dist_triangle <- function(similr, cancer_name) {
  # this function gets a distances triangle and turns it into a long format
  x <- as.matrix(similr)
  r_ind <- cbind(as.data.frame(rownames(x)), 1:15)
  names(r_ind) <- c("chap_row", "row")
  c_ind <- cbind(as.data.frame(rownames(x)), 1:15)
  names(c_ind) <- c("chap_col", "col")
  
  ind <- which(upper.tri(x, diag = FALSE), arr.ind = TRUE)
  sims_cancer <- as.data.frame(cbind(ind, value=x[ind])) %>% 
    left_join(y=r_ind, by=c("row")) %>%
    left_join(y=c_ind, by=c("col")) %>% 
    select(row=chap_row, col=chap_col, similarity=value) %>% 
    add_column(cancer=cancer_name)
  return(sims_cancer)
}

#-------- load the networks from an excel file, and calculate similarity --------
excel_path <- "HPC/binari_validated_corrs.xlsx"

# do the convert for every cancer
sheet_names <- excel_sheets(excel_path)

networks <- list()
sim_triangles <- list()
all_simlrs <- matrix(0, nrow = 105, ncol = 0)
all_simlrs_long <- NULL
for (name in sheet_names) {
  x <- as.data.frame(read_excel(excel_path, sheet = name, col_names = TRUE))
  x2 <- x[,-1]
  rownames(x2) <- x[,1]
  networks[[name]] <- x2
  
  # generate similarity between the chaperons in this cancer
  res <- vegdist(x2, method="jaccard")
  simlr <- 1-res
  
  all_simlrs <- cbind(all_simlrs, simlr)
  sim_triangles[[name]] <- simlr
  
  long <- long_format_from_dist_triangle(simlr, name)
  all_simlrs_long <- rbind(all_simlrs_long, long)
}

colnames(all_simlrs) <- sheet_names
write.csv(all_simlrs, file = "output/jaccard_values_per_cancer.csv", row.names = FALSE)
write.csv(all_simlrs_long, file = "output/jaccard_values_per_cancer_long_format.csv", row.names = FALSE)

#-------- visualize the results --------
all_simlrs <- read.csv("output/jaccard_values_per_cancer.csv")

# melt the data to a long format
mlt_sim <- as.data.frame(melt(all_simlrs))
g1 <- ggplot(mlt_sim, aes(x=variable,y=value))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  xlab("") + ylab("Similarity index")
g1

# distribute similarities
p<-ggplot(mlt_sim, aes(x=value)) + 
  geom_histogram(bins=40, color="black") + 
  labs(x="Jaccard similarity index")
p

# distribute similarities by cancer
p2<-ggplot(mlt_sim, aes(x=value, fill= variable)) + 
  geom_histogram(bins=20, color="black", position="identity", alpha=0.5) + 
  labs(x="Jaccard similarity index", fill = "Cancer")
p2

ggsave("output/paper_figures/cancer_jaccard_boxplot.pdf", g1)

#-------- investigate who is in the upper tail --------
all_simlrs_long <- read.csv("output/jaccard_values_per_cancer_long_format.csv")

in_tail <- all_simlrs_long %>% filter(similarity > 0.45)

# distribution in tail  by cancer
p2<-ggplot(in_tail, aes(x=similarity, fill= cancer)) + 
  geom_histogram(bins=20, color="black", position="stack", alpha=0.5) + 
  labs(x="Jaccard similarity index", fill = "Cancer")
p2

n_distinct(c(in_tail$col, in_tail$row))
table(c(in_tail$col, in_tail$row))


rrr <- igraph::graph_from_data_frame(in_tail, directed = FALSE)

png(filename = "output/cancer_similarity_tail_chapchap_network.png")
plot(rrr)
dev.off()


