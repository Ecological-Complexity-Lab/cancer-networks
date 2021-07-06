# This  script uses expression data to create a binari network
# by going through every correlation and validating it using bootstraps - 
# meaning to take a subset of samples and check that the correlation is still 
# significant, 1000 times.


# bash intro
#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())
Sys.setlocale(category = "LC_ALL", locale = "") # done so that write.csv works, but doesnt seem to really change any locale
 
library(magrittr)
library(tidyverse)

# read args
if (length(commandArgs(trailingOnly=TRUE))==0) {
   stop('No arguments were found!')
} else {
   args <- commandArgs(trailingOnly=TRUE)
   cancer_id <- as.numeric(args[1])
   if ((is.na(cancer_id)) || (cancer_id<1) || (cancer_id>12)) {
     stop('Argumet must be a number in the range 1-12')
   }
 }

# ------------ consts -----------------
set.seed(seed = 14412)
JOB_ID <- Sys.getenv("JOB_ID") 

SAMPLE_NUM <- 288
SIM_NUM <- 1000
EXP_FOLDER_LOC <- 'DATA/expression/'
OUTPUT_FOLDER <- 'DATA/binari/'

TISSUE_LIST <- c('Breast Invasive Carcinoma',
                 'Colon_Adenocarcinoma',
                 'Head_and_Neck_Squamous_Cell_Carcinoma',
                 'Kidney_Renal_Clear_Cell_Carcinoma',
                 'Kidney_Renal_Papillary_Cell_Carcinoma',
                 'Liver_Hepatocellular_Carcinoma',
                 'Lung_Adenocarcinoma',
                 'Lung_Squamous_Cell_Carcinoma',
                 'Prostate_Adenocarcinoma',
                 'Stomach_Adenocarcinoma',
                 'Thyroid_Carcinoma',
                 'Uterine_Corpus_Endometrial_Carcinoma')
tissue <- TISSUE_LIST[cancer_id]

# ------------ functions ---------------------------

LOG <- function(s, appnd = T) {
    write_lines(s, paste("logs/", JOB_ID, tissue,'log.txt', sep='_'), append = appnd)
}

Calc_bonf <- function(alfa, nrows, ncols) {
  p <- alfa/(nrows*ncols)
  return(p)
}

bootstrap_corrs_for_tissue <- function(x, y, num_bootings) {
  # get a correlation and return a list of p values received from 
  # the correlations of the bootstraps 
  p_list <- numeric(num_bootings)
  
  for (i in 1:num_bootings) {
    samples_to_include <- sample(1:length(x), size = SAMPLE_NUM, replace = FALSE)
    
    sub_x <- x[samples_to_include]
    sub_y <- y[samples_to_include]
    
    # calc subset corr test
    corr_result <- cor.test(sub_x, sub_y, method = "spearman", exact=FALSE)
    
    # if the correlation test failed
    if (is.na(corr_result$estimate)) { LOG("Bootstrap correlation failed.") }
    
    # save the r value
    p_list[i] <- corr_result$p.value
  }
  # return the list of Ps we found
  return(p_list)
}

Is_obs_p_valid <- function(obs_p, x, y, booting_num, bonf_value, should_show=FALSE) {
  # use bootstraps to validate the specific correlation that we check
  # the func returns 1 if the correlation is valid and 0 otherwise
  # valid = significant 95% of the simulations
  corr_p_list <- bootstrap_corrs_for_tissue(x, y, booting_num)
  corr_p_list <- sort(corr_p_list)
  
  # plot if required
  if (should_show) {
    hist(corr_p_list)
    abline(v=bonf_value, col="blue")
    abline(v=obs_p, col="red")
  }
  
  # calculate the proportion of the significant simulations 
  signf_ps <- length(corr_p_list[corr_p_list < bonf_value])
  percent <- signf_ps/booting_num

  # test the observed compared to the simulations
  output <- 1
  if (percent < 0.95) { output <- 0 } 

  return(output)
}

get_valid_corrs_per_tissue <- function(tisssue, chaps, prots) {
  # this function returns a binari table of the valid correlations that pass
  # validation using bootstraps tests
  
  # get expression data (all samples)
  exppath <- paste(EXP_FOLDER_LOC, tisssue, ".tsv", sep = "")
  exp <- read.table(exppath, header = FALSE, row.names = 1)
  
  LOG("Data was read from the expression file.")
  
  # counter to document how many corrs violated the null hypothesis 
  sim_number <- SIM_NUM
  count_diff_corrs <- 0
  all_passing_Rs <- 0
  bonf_value <- Calc_bonf(0.05, nrow(chaps), nrow(prots))
  
  # create empty table
  m <- matrix(0, nrow = nrow(chaps), ncol = nrow(prots))
  
  binari_tbl <- data.frame(m)
  rownames(binari_tbl) <- chaps$Symbol
  colnames(binari_tbl) <- prots$ENSID

  LOG(paste("Tissue name: ", tisssue, ".", sep = ""))
  LOG(paste("Number of samples: ", ncol(exp), ".", sep = ""))
  LOG(paste("Expected number of samples in subset: ", SAMPLE_NUM, ".", sep = ""))
  LOG(paste("Number of bootstrapings: ", sim_number, ".", sep = ""))
  LOG(paste("Bonferroni cutoff value: ", bonf_value, ".", sep = ""))
  
  # for every chapXprot combination, check if a passing R value is like bootstraps corrs
  for (chap in 1:nrow(chaps)) {
    chapID <- chaps[chap, "ENSID"]
    chapSymb <- chaps[chap, "Symbol"]
    x <- unlist(exp[rownames(exp) == chapID, ])
    
    for (prot in 1:nrow(prots)) {
      # prepare prot data for corr
      protID <- prots[prot, "ENSID"]
      protSymb <- prots[prot, "Symbol"]
      y <- unlist(exp[rownames(exp) == protID, ])
      
      if (length(x) != length(y)) { 
        LOG(paste("Error in test correlation of:", chapSymb, "and", protSymb)) 
      }
      
      # calc observed corr test
      corr_result <- cor.test(x, y, method = "spearman", exact=FALSE)
      
      # if the correlation test failed
      if (is.na(corr_result$estimate)) {next}
      
      # get observed values
      obs_r <- corr_result$estimate
      obs_p <- corr_result$p.value
      
      # if the correlation is not significant then go to the next loop iteration
      if (obs_p > bonf_value) { next }
      # if the correlation is not positive then go to the next loop iteration
      if (obs_r < 0) { next }
      
      # chack if obs r is different (1) or not (0) to a bootstrap average
      is_valid <- Is_obs_p_valid(obs_p, x, y, sim_number, bonf_value)
      binari_tbl[chapSymb, protID] <- is_valid
      
      # counters
      count_diff_corrs <- count_diff_corrs + is_valid
      all_passing_Rs <- all_passing_Rs + 1
    }
  }
  
  diff_percent <- (count_diff_corrs/all_passing_Rs)*100
  
  LOG(paste("all passing correlations: ", all_passing_Rs, ".", sep = ""))
  LOG(paste("validated correlations: ", count_diff_corrs, ".", sep = ""))
  LOG(paste("validated percentage: ", diff_percent, "% ", sep = ""))
  LOG(paste("Total correlations: ", nrow(chaps)*nrow(prots), ". ", sep = ""))
  
  return(binari_tbl)
}


# ------------- run --------------

# essumptions: expression tables with all samples are available. 
#              list of mito prot and mito chaps are available

LOG('START RUN LOG', appnd = F)
LOG('====================\n')

# read the proteins metadata
prots_meta <- read.table("Mito_genes.tab", sep="\t", header=TRUE,
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
prots_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]

# run the validation for a single tissue
bi_tbl <- get_valid_corrs_per_tissue(tissue, chaps_meta, prots_meta)

LOG(paste("Validation done. Output data type:", typeof(bi_tbl), class(bi_tbl)))

# save table
output_pth <- paste(OUTPUT_FOLDER, tissue,".csv", sep = "")
LOG(paste('Saving results to csv:', output_pth))

write.csv(bi_tbl, output_pth, row.names = TRUE)
LOG('Binari table was saved as csv file')

LOG('\n====================')
LOG('FINISH RUN LOG')

