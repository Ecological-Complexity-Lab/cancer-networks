#------------------------------
# build a table that shows the percent of the folding potential that is
# activated for each chaperon for every cancer.
# 
#------------------------------


#------ includes ------
library(readxl)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

#----------- load the networks from an excel file --------
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

#----------- calc chaperon degree per cancer -----------------

# get the list of proteins ENSID as rownames, SYBOL as first row
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
prots_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]
Symbol <- chaps_meta[, 4]
degree_table_chap <- as.data.frame(Symbol)

# go over the cancers
for (i in 1:length(networks)) {
  netname <- names(networks)[i]
  net <- networks[[i]]
  
  # go over the proteins
  chap_degrees <- c()
  for (chap in 1:nrow(chaps_meta)) {
    SYMB <- chaps_meta[chap, "Symbol"]
    # calc relative degree
    
    chap_degrees[SYMB] <- sum(net[SYMB, ])
  }
  
  # put each cancer in a new column
  curr_cancer_degrees <- data.frame(chap_degrees)
  new_col_names <- c(colnames(degree_table_chap), netname)
  
  degree_table_chap <- merge(degree_table_chap, curr_cancer_degrees,
                             by.x="Symbol", by.y="row.names", all.x=TRUE)
  degree_table_chap
  
  colnames(degree_table_chap) <- new_col_names
}
# save the table to a file
write.csv(degree_table_chap, "output/degree_chaperons.csv", row.names = FALSE)

#----------- create union of the cancer networks -------------
union_table <- networks[[1]]
union_table[,] <- 0 

for (i in 1:length(networks)) {
  net <- networks[[i]]
  union_table <- union_table + net
}
union_table_mat <- as.matrix(union_table)

#----------- flat the union and calc folding potential -------
binari_union <- 1*(union_table>0)

chap_potential <- rowSums(binari_union, na.rm = FALSE, dims = 1)

# divide by total protein number to find the percent from all mito chaps
puf <- as.data.frame(chap_potential)/nrow(prots_meta)
puf$chaperons <- rownames(puf)

g <- ggplot(puf, aes(x=reorder(chaperons, -chap_potential),y=chap_potential)) +
     geom_bar(stat="identity", fill="steelblue", width=0.5)+
     labs(title="Chaperon folding potential", x ="chaperons",
          y = "folding potential (%)") + coord_flip()

ggsave("output/chap_fold_potential.pdf", g)

puf <- puf[,c(2,1)]
names(puf) <- c("name", "potential")
puf["amount"] <- puf$potential * nrow(prots_meta)
write.csv(puf, file = "output/folding_potential.csv", row.names = FALSE)

#----------- calculate percent of potential used ----------
row.names(degree_table_chap) <- degree_table_chap$Symbol
folding_percent <- degree_table_chap[,-1]

for (chapp in chaps_meta$Symbol) {
  chap_pot <- chap_potential[chapp]
  folding_percent[chapp, ] <- folding_percent[chapp,]/chap_pot
}

write.csv(folding_percent, file = "output/chap_folding_percent.csv")

#----------- generate heatmap from data ----------
pheatmap(folding_percent, cutree_rows = 2, cutree_cols = 2,
         clustering_method = "ward.D", clustering_distance_rows = "manhattan",
         filename = "output/folding_percentage_clustered_heatmap.pdf",
         color = colorRampPalette(brewer.pal(n = 7, name = "YlGnBu"))(100),
         main = "Folding percentage by cancer", angle_col = 45)

# reorder heat map by row sum and col sum:
chp_sum <- rowSums(folding_percent)
chap_ordered <- t(folding_percent[order(chp_sum,decreasing=T),])
cncr_sum <- rowSums(chap_ordered)
chap_ordered <- t(chap_ordered[order(cncr_sum,decreasing=T),])

pheatmap(chap_ordered, cluster_rows = F, cluster_cols = F,
         filename = "output/folding_percentage_nestedness_heatmap.pdf",
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100),
         main = "Folding percentage by cancer", angle_col = 45 )

# test color: display.brewer.pal(n=7,name="YlGnBu")
