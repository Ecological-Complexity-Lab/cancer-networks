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

#-------- load the networks from an excel file, and calculate similarity --------
excel_path <- "HPC/binari_validated_corrs.xlsx"

# do the convert for every cancer
sheet_names <- excel_sheets(excel_path)

networks <- list()
sim_triangles <- list()
all_simlrs <- matrix(0, nrow = 105, ncol = 0)

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
}

colnames(all_simlrs) <- sheet_names
write.csv(all_simlrs, file = "output/jaccard_values_per_cancer.csv", row.names = FALSE)


#-------- visualize the results --------
# melt the data to a long format
mlt_sim <- as.data.frame(melt(all_simlrs))
g1 <- ggplot(mlt_sim, aes(x=Var2,y=value))+
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
p2<-ggplot(mlt_sim, aes(x=value, fill= Var2)) + 
  geom_histogram(bins=20, color="black", position="identity", alpha=0.5) + 
  labs(x="Jaccard similarity index", fill = "Cancer")
p2

