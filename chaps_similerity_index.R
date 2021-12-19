#------------------------------
# use similarity indexes to find the differences between chaperons
# by comparing between the cancers for each chaperone
#------------------------------

#-------- includes --------
library(readxl)
library(vegan)
library(ggplot2)
library(reshape2)
library(tidyr)

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
# Jaccard index is computed as 2B/(1+B), where B is Bray-Curtis *dissimilarity*.
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

#-------- build a DF for dependency visualization --------
potentials <- read.csv("output/folding_potential.csv", row.names = 1)
expression <- read.csv("output/chap_median_expressions.csv", row.names = 1) # TODO

# make a new df with: 
# name | potential | median_similarity | 25q | 75q | med_realized | 25q | 75q
f_s <- t(do.call(cbind, lapply(as.data.frame(fold_per), summary)))
s_s <- t(do.call(cbind, lapply(as.data.frame(all_simlrs), summary)))
e_s <- t(do.call(cbind, lapply(as.data.frame(log10(t(expression))), summary)))


f_s <- as.data.frame(f_s) %>% dplyr::select("1st Qu.", "Median" ,"3rd Qu.") %>% 
  dplyr::rename(fold_q1 = `1st Qu.`, fold_med = Median, fold_q3 = `3rd Qu.`)
s_s <- as.data.frame(s_s) %>% dplyr::select("1st Qu.", "Median" ,"3rd Qu.") %>% 
  dplyr::rename(sim_q1 = `1st Qu.`, sim_med = Median, sim_q3 = `3rd Qu.`)
e_s <- as.data.frame(e_s) %>% dplyr::select("1st Qu.", "Median" ,"3rd Qu.") %>% 
  dplyr::rename(exp_q1 = `1st Qu.`, exp_med = Median, exp_q3 = `3rd Qu.`)

all_stats <- merge(potentials, f_s, by="row.names", all=TRUE)
all_stats <- merge(all_stats, s_s, by.x="Row.names" ,by.y="row.names", all=TRUE)
all_stats <- merge(all_stats, e_s, by.x="Row.names" ,by.y="row.names", all=TRUE)
rownames(all_stats) <- all_stats$Row.names
all_stats <- all_stats[2:length(all_stats)]
all_stats

# -------- visualize the two (niche ~ similarity) as a dependency --------
sml <- all_stats$sim_med
m_q1 <- all_stats$sim_q1
m_q3 <- all_stats$sim_q3
fol <- all_stats$fold_med
f_q1 <- all_stats$fold_q1
f_q3 <- all_stats$fold_q3
png(filename = "output/similarity_potential_scatter.png")
plot(y = fol, x = sml,
     xlim=range(c(sml-m_q1, sml+m_q3)),
     ylim=range(c(fol-f_q1, fol+f_q3)),
     pch=19, ylab="realized folding", xlab="similarity",
     main="jaccard vs realized folding dapancdency"
)
# we draw arrows with very special "arrowheads" as whiskers
arrows(sml-m_q1, fol, sml+m_q3, fol, length=0.05, angle=90, code=3) # siml whiskers
arrows(sml, fol-f_q1, sml, fol+f_q3,  length=0.05, angle=90, code=3) # fold whiskers
dev.off()

# -------- visualize the similarity vs folding potential --------
pot <- all_stats$potential
med <- all_stats$sim_med
m_q1 <- all_stats$sim_q1
m_q3 <- all_stats$sim_q3
png(filename = "output/similarity_realized_scatter.png")
plot(y = pot, x = med,
     xlim=range(c(med-m_q1, med+m_q3)),
     pch=19, ylab="Potential", xlab="similarity",
     main="jaccard vs folding potential dapancdency"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(med-m_q1, pot, med+m_q3, pot, length=0.05, angle=90, code=3)
text(med, pot, all_stats$Row.names, pos=1)
dev.off()

#----------- visualize (niche ~ potential) as a dependency ----------
# this doesn't really belong here logically but all that is needed for it is here.
pot <- all_stats$potential
med <- all_stats$fold_med
m_q1 <- all_stats$fold_q1
m_q3 <- all_stats$fold_q3
png(filename = "output/similarity_realized_scatter.png")
plot(y = pot, x = med,
     xlim=range(c(med-m_q1, med+m_q3)),
     pch=19, ylab="Potential", xlab="realized folding",
     main="realized folding vs folding potential dapancdency"
)
arrows(med-m_q1, pot, med+m_q3, pot, length=0.05, angle=90, code=3)
text(med, pot, all_stats$Row.names, pos=1)
dev.off()


#----------- visualize (similarity ~ expression) as a dependency ----------
sml <- all_stats$sim_med
m_q1 <- all_stats$sim_q1
m_q3 <- all_stats$sim_q3
exp <- all_stats$exp_med
e_q1 <- all_stats$exp_q1
e_q3 <- all_stats$exp_q3
png(filename = "output/similarity_expression_scatter.png")
plot(y = exp, x = sml,
     xlim=range(c(sml-m_q1, sml+m_q3)),
     ylim=range(c(exp-e_q1, exp+e_q3)),
     pch=19, ylab="Log10 expression", xlab="similarity",
     main="jaccard vs log10 expression levels dapancdency"
)
# we draw arrows with very special "arrowheads" as whiskers
arrows(sml-m_q1, exp, sml+m_q3, exp, length=0.05, angle=90, code=3) # siml whiskers
arrows(sml, exp-e_q1, sml, exp+e_q3,  length=0.05, angle=90, code=3) # expression whiskers
dev.off()

#------- test correlation between median similarity and expression ----------
st <- cor.test(sml, exp, method="spearman", exact=FALSE)
obs_pval <- st$p.value
obs_rval <- st$estimate

st