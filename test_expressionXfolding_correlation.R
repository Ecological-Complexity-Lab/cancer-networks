#--------------- test_expressionXfolding_correlation.R ---------------
# read the *expression* and *folding percent*
# run mantel test on it
#------------------------------

#-------- includes --------
library(readxl)
library(ggplot2)
library(tibble)
library(tidyverse)
library(magrittr)


#-------- read data -------------
# read folding data
folding_df <- read.csv("output/chap_realized_niche.csv", row.names = 1)

# read expression data
expr_df <- read.csv("output/chap_median_expressions.csv", row.names = 1)

# make sure the cols and rows are ordered is the same way
folding_df <- folding_df[,order(colnames(folding_df))]
folding_df <- folding_df[order(rownames(folding_df)),]
expr_df <- expr_df[,order(colnames(expr_df))]
expr_df <- expr_df[order(rownames(expr_df)),]

#-------- test correlation ׳with log10 ------------
exp_vec <- as.vector(as.matrix(expr_df))
fold_vec <- as.vector(as.matrix(folding_df))
log_exp_vec <- log10(exp_vec)
can_vec_2 <- rep(colnames(expr_df), each = length(rownames(expr_df)))
ch_vec_2 <- rep(rownames(expr_df), times=length(colnames(expr_df)))

tbl_all_2 <- tibble(expression=log_exp_vec, fold=fold_vec, cancer=can_vec_2, chap=ch_vec_2)

write_csv(tbl_all_2, "output/data/chap_rn_and_log_exp.csv")

g <- ggplot(tbl_all_2, aes(x=expression, y=fold, color=chap))+
  geom_point(size=2)+
  #ggtitle("Realized Niche over Log10 Median Expression levels")+
  ylab("Realized Niche (%)")+
  xlab("Log10 median expression level") + 
  labs(color="Chaperon") 
ggsave("output/paper_figures/log10_rn_vs_exp.pdf", g)


st <- cor.test(log_exp_vec, fold_vec, method="spearman", exact=FALSE)
obs_pval <- st$p.value
obs_rval <- st$estimate

#-------- test correlation ׳with log10 for each chap ------------
currs <- tibble(chaperon=character(),
                p_value=numeric(), 
                r_value=numeric())
for (chp in rownames(expr_df)) {
  data <- tbl_all_2[tbl_all_2$chap == chp,]
  
  st <- cor.test(data$expression, data$fold, method="spearman", exact=FALSE)
  obs_pval <- st$p.value
  obs_rval <- st$estimate
  
  currs %<>% add_row(chaperon=chp, p_value=obs_pval, r_value=obs_rval)
  
  ggplot(data, aes(x=expression, y=fold))+
    geom_point(size=3)+
    ggtitle(paste(chp ,"realized Niche over Log10 Median Expression levels"))+
    ylab("Realized Niche (%)")+
    xlab("Log10 median expression level")
} 

# the correlations used in the paper:
currs # non show a significant correlation but YME1L1 with p=0.00824 r=-0.720
# YME1l1 also have the smallest potential

ggplot(tbl_all_2, aes(x=expression, y=fold, color=chap)) +
     geom_point() + 
     ylab("Realized Niche (%)") + xlab("Log10 median expression level")+
     ggtitle("Realized Niche over Median Expression levels") +
    theme(
      legend.position="none",
      panel.spacing = unit(0.5, "lines"),
      strip.text.x = element_text(size = 10),
    ) +
     facet_wrap(~ chap)


#-------- run permutations to validate correlation ------------
# this validates the correlation of all the expressions together, 
# and not per chaperon.

# shuffle expression, fix fold, re-correlate
SIM_NUM <- 1000
set.seed(42)

pvals <- vector(mode='numeric',length=SIM_NUM)
rvals <- vector(mode='numeric',length=SIM_NUM)
bigger_then_obs <- 0
smaller_then_obs <- 0
for (i in 1:SIM_NUM) {
  print(i)
  perm <- sample(log_exp_vec, length(log_exp_vec), replace = FALSE)
  corr_st <- cor.test(perm, fold_vec, method="spearman", exact=FALSE)
  pvals[i] <- corr_st$p.value
  rvals[i] <- corr_st$estimate
  
  smaller_then_obs <- smaller_then_obs + (corr_st$p.value < obs_pval)
}

print("permutetion test p-value is:")
p <- smaller_then_obs/SIM_NUM
p

#save histogram plot for permutation test
pdf(file="output/expression_folding_correlation_results.pdf")

my_hist=hist(pvals , breaks=40  , plot=F)
plot(my_hist, col=rgb(0.2,0.2,0.2,0.2) , border=F , 
     main="Permutations p-values" , xlab="p-value", ylim=c(0,30) )
abline(v=obs_pval,col="blue",lwd=2, lty=6)
text(0.8, 30, paste("Permutation test p-value:", p))
text(0.8, 28.5, paste("Observed p-value:", round(obs_pval, 4)))

dev.off()
