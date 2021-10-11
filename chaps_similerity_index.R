#------------------------------
# use similarity indexes to find the differences between chaperons
# 
#------------------------------

#-------- includes --------
library(readxl)
library(vegan)
library(ggplot2)
library(reshape2)

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


#-------- build cancerXproteins network per chap --------
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
prots_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]

# build empty df
empty_net <- data.frame(matrix(nrow=length(sheet_names), 
                               ncol = length(prots_meta$ENSID)))
colnames(empty_net) <- prots_meta$ENSID
row.names(empty_net) <- sheet_names


chap_nets <- list()
for (chap in chaps_meta$Symbol) {
  chap_net <- empty_net
  
  for (cancr in sheet_names) {
    cncr_df <- networks[[cancr]]
    
    chap_net[cancr,] <- cncr_df[chap,]
  }
  
  # add to network list
  chap_nets[[chap]] <- chap_net
}


#-------- use vegan::vegdist() --------
# Jaccard index is computed as 2B/(1+B), where B is Brayâ--Curtis *dissimilarity*.
# so use the 1-P value to use the similarity.

all_simlrs <- matrix(0, nrow = 66, ncol = 0)
for (chapp in names(chap_nets)) {
  net <- chap_nets[[chapp]]
  res <- vegdist(net, method="jaccard")
  simlr <- 1-res
  
  all_simlrs <- cbind(all_simlrs, simlr)
}

colnames(all_simlrs) <- names(chap_nets)


#-------- visualize the results --------
# melt the data to a long format
mlt_sim <- as.data.frame(melt(all_simlrs))
g1 <- ggplot(mlt_sim, aes(x=Var2,y=value))+
      geom_boxplot(outlier.colour="black", outlier.shape=16,
                   outlier.size=2, notch=FALSE)
g1


#-------- load realized niche to visualize together --------
fold_per <- t(read.csv("output/chap_folding_percent.csv", row.names = 1))

mlt_per <- as.data.frame(melt(fold_per))
g2 <- ggplot(mlt_per, aes(x=Var2,y=value))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)
g2

# prpepare data to combine
mlt_per$Var1 <- "Realized Niche"
mlt_sim$Var1 <- "Similarity Index"
all <- rbind(mlt_per, mlt_sim)
  
# build the plot to work togather
g3 <- ggplot(all, aes(x=Var2, y=value, fill=Var1)) +
      geom_boxplot(outlier.colour="black", outlier.shape=16,
                   outlier.size=2, notch=FALSE) +
      scale_y_continuous(name = "Realized Niche (%)",# Add a second axis:
                         sec.axis = sec_axis(trans=~., 
                                             name="Similerity (b diversity")) +
      labs(fill="Indexes:", x=NULL) + theme(legend.position="top", 
                                    axis.text.x=element_text(angle=45, hjust=1))
g3


ggsave("output/similarity_and_percent_boxplot.pdf", g3)


# try a violin chart
g4 <- ggplot(all, aes(x=Var2, y=value, fill=Var1)) +
  geom_violin() +
  scale_y_continuous(name = "Realized Niche (%)",# Add a second axis:
                     sec.axis = sec_axis(trans=~., 
                                         name="Similerity (b diversity")) +
  labs(fill="Indexes:", x=NULL) + theme(legend.position="top", 
                                        axis.text.x=element_text(angle=45, hjust=1))
g4

ggsave("output/similarity_and_percent_violin.pdf", g4)


#------- test correlation between median similarity and folding ----------
f <- apply(fold_per,2,median)
s <- apply(all_simlrs,2,median)
f <- f[order(names(f))]
s <- s[order(names(s))]

st <- cor.test(s, f, method="spearman", exact=FALSE)
obs_pval <- st$p.value
obs_rval <- st$estimate

st