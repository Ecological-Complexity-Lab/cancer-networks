#------------------------------
# read the expression levels for each chaperon, 
# and compare it to its the potential folding
#------------------------------

#-------- includes --------
library(readxl)
library(tools)
library(pheatmap)



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


#-------- load chaperons expression data from excel files --------
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)

# read each file in the path and build a expXcancer for each chap
folder_path <- "HPC/DATA/expression"
files <- list.files(folder_path)

exp_tables <- list()
for (file in files) {
  chap_exp <- read.table(paste(folder_path,"/",file, sep = ""), 
                         header = FALSE)
  # filter only the chap data
  chap_exp <- chap_exp[chap_exp$V1 %in% chaps_meta$ENSID, ]
  rownames(chap_exp) <- chaps_meta$Symbol[match(chap_exp$V1, chaps_meta$ENSID)]
  
  cncr_name <- cancer_names[file_path_sans_ext(file)]
  exp_tables[[cncr_name]] <- chap_exp[,-1]
}


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


#-------- visualize as heat map --------
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
