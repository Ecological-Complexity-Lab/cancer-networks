#------------------------------
# read the expression levels for each chaperon, 
# and protein, to have an expression summery for them in each cancer.
#------------------------------

#-------- includes --------
library(readxl)
library(tools)
library(pheatmap)
library(plyr)

#-------- make a name 'dictionary' ---------
# to convert long name to short
keys <- c("Breast Invasive Carcinoma", 
          "Colon_Adenocarcinoma", 
          "Head_and_Neck_Squamous_Cell_Carcinoma", 
          "Kidney_Renal_Clear_Cell_Carcinoma",
          "Kidney_Renal_Papillary_Cell_Carcinoma",
          "Liver_Hepatocellular_Carcinoma",
          "Lung_Adenocarcinoma",
          "Lung_Squamous_Cell_Carcinoma",
          "Prostate_Adenocarcinoma",
          "Stomach_Adenocarcinoma",
          "Thyroid_Carcinoma",
          "Uterine_Corpus_Endometrial_Carcinoma")
cancer_names <- c("BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LIHC",
                  "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")
names(cancer_names) <- keys

#-------- functions --------
exp_data_for_gene_list_ENSID <- function(gene_list_data) {
  # read each file in the path and build a expXcancer for each gene
  folder_path <- "HPC/DATA/expression"
  files <- list.files(folder_path)
  
  ex_tables <- list()
  for (file in files) {
    chap_exp <- read.table(paste(folder_path,"/",file, sep = ""), 
                           header = FALSE)
    # filter only the chap data
    chap_exp <- chap_exp[chap_exp$V1 %in% gene_list_data$ENSID, ]
    rownames(chap_exp) <- gene_list_data$Symbol[match(chap_exp$V1, gene_list_data$ENSID)]
    
    cncr_name <- cancer_names[file_path_sans_ext(file)]
    ex_tables[[cncr_name]] <- chap_exp[,-1]
  }
  return(ex_tables)
}

#-------- load chaperons expression data from excel files --------
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)

exp_tables <- exp_data_for_gene_list_ENSID(chaps_meta)

#-------- create a table of mean/median expressions for chap per cancer ------
exp_mean <- data.frame(matrix(nrow=length(chaps_meta$Symbol), 
                              ncol = length(exp_tables)))
row.names(exp_mean) <- chaps_meta$Symbol
colnames(exp_mean) <- names(exp_tables)

exp_median <- data.frame(matrix(nrow=length(chaps_meta$Symbol), 
                                ncol = length(exp_tables)))
row.names(exp_median) <- chaps_meta$Symbol
colnames(exp_median) <- names(exp_tables)

for (cncr in names(exp_tables)) {
  exp_tbl <- as.matrix(exp_tables[[cncr]])
  
  for (chp in chaps_meta$Symbol) {
    avrg <- mean(exp_tbl[chp,])
    medn <- median(exp_tbl[chp,])
    
    exp_mean[chp, cncr] <- avrg
    exp_median[chp, cncr] <- medn
  }
}

write.csv(exp_mean, file = "output/chap_mean_expressions.csv")
write.csv(exp_median, file = "output/chap_median_expressions.csv")


#-------- visualize chaps as heat map --------
pheatmap(exp_median, clustering_method = "ward.D", clustering_distance_rows = "manhattan",
         filename = "output/median_chap_expression_heatmap.pdf",
         color=colorRampPalette(c("white", "orange"))(100),
         main = "Median Chaperon Expression", angle_col = 315,display_numbers = TRUE)

pheatmap(exp_mean, clustering_method = "ward.D", clustering_distance_rows = "manhattan",
         filename = "output/mean_chap_expression_heatmap.pdf",
         color=colorRampPalette(c("white", "chartreuse3"))(100),
         main = "Mean Chaperon Expression", angle_col = 315,display_numbers = TRUE)


# reorder heat map by row sum and col sum:
chp_sum <- rowSums(exp_median)
chap_ordered <- t(exp_median[order(chp_sum,decreasing=T),])
cncr_sum <- rowSums(chap_ordered)
chap_ordered <- t(chap_ordered[order(cncr_sum,decreasing=T),])


pheatmap(log10(chap_ordered), cluster_rows = F, cluster_cols = F,
         color=colorRampPalette(c("white", "orange"))(100),
         filename = "output/Median_log10_chap_expression_heatmap.pdf",
         main = "Log10 of Median Chaperon Expression", angle_col = 315,
         display_numbers = TRUE)


#-------- load protein expression data from excel files -------
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
prots_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]

exp_p_tables <- exp_data_for_gene_list_ENSID(prots_meta)

#-------- create a table of mean/median expressions for protein per cancer ------
exp_mean_prot <- data.frame(matrix(nrow=length(prots_meta$Symbol), 
                              ncol = length(exp_p_tables)))
row.names(exp_mean_prot) <- prots_meta$Symbol
colnames(exp_mean_prot) <- names(exp_p_tables)

exp_med_prot <- data.frame(matrix(nrow=length(prots_meta$Symbol), 
                                ncol = length(exp_p_tables)))
row.names(exp_med_prot) <- prots_meta$Symbol
colnames(exp_med_prot) <- names(exp_p_tables)

for (cncr in names(exp_p_tables)) {
  exp_tbl <- as.matrix(exp_p_tables[[cncr]])
  
  for (prt in prots_meta$Symbol) {
    avrg <- mean(exp_tbl[prt,])
    medn <- median(exp_tbl[prt,])
    
    exp_mean_prot[prt, cncr] <- avrg
    exp_med_prot[prt, cncr] <- medn
  }
}

write.csv(t(exp_mean_prot), file = "output/prot_mean_expressions.csv")
write.csv(t(exp_med_prot), file = "output/prot_median_expressions.csv")

#-------- visualize prots as histogram per cancer --------
exp_med_prot <- read.csv("output/prot_median_expressions.csv", row.names = 1)
prot_exp_nolog <- melt(t(exp_med_prot))
prot_exp <- melt(t(log10(exp_med_prot)))
head(prot_exp)
ggplot(prot_exp_nolog, aes(x=value, fill=Var2, color=Var2)) +
  geom_histogram(position="identity", alpha=0.5)+
  xlab("protein expression") + ylab("Count")


ggplot(prot_exp, aes(x=value, fill=Var2, color=Var2)) +
  geom_histogram(position="identity", alpha=0.5)+
  xlab("Log10 of protein expression") + ylab("Count")

mu <- ddply(prot_exp_nolog, "Var2", summarise, grp.mean=mean(value))

big_only <- filter(prot_exp, value > 2 | value == 2)
ggplot(big_only, aes(x=value, fill=Var2, color=Var2)) +
  geom_histogram(position="identity", alpha=0.5)+
  xlab("Log10 of protein expression") + ylab("Count") + 
  ggtitle("Protein expression of 100 and above")

#-------- visualize prots as heatmap --------

pheatmap(exp_med_prot, clustering_method = "ward.D", clustering_distance_rows = "manhattan",
         filename = "output/median_prot_expression_heatmap.pdf",
         color=colorRampPalette(c("white", "orange"))(100),
         main = "Median Protein Expression", angle_col = 315, show_rownames = FALSE)

# reorder heat map by row sum and col sum:
prt_sum <- rowSums(exp_med_prot)
prot_ordered <- t(exp_med_prot[order(prt_sum,decreasing=T),])
cncr_sum <- rowSums(prot_ordered)
prot_ordered <- t(prot_ordered[order(cncr_sum,decreasing=T),])

pheatmap(prot_ordered, cluster_rows = F, cluster_cols = F,
         color=colorRampPalette(c("white", "orange"))(100),
         filename = "output/Median_prot_expression_heatmap_unclustered.pdf",
         main = "Log10 of Median Protein Expression", angle_col = 315, show_rownames = FALSE)

# use log10
logs <- log10(prot_ordered)
logs[is.infinite(logs)] <- -3
pheatmap(logs, cluster_rows = F, cluster_cols = F,
         color=colorRampPalette(c("white", "orange"))(100),
         filename = "output/Median_log10_prot_expression_heatmap.pdf",
         main = "Log10 of Median Protein Expression", angle_col = 315, show_rownames = FALSE)
