# ------- produce_figures.r ------
# this script is meant to pipline the making of figures for the paper.
# just read the data and make a figure. 
# as a general rule only paper figures are produced here.
# --------

# -------- includes --------
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(readr)
library(tidyverse)

source("functions.r")

# ----- consts -----
output_folder <- "output/paper_figures/"

# potential and nestedness --------
# potential bars:
puf <- read.csv("output/folding_potential.csv")

g <- ggplot(puf, aes(x=reorder(name, -potential), y=potential)) +
  geom_bar(stat="identity", fill="steelblue", width=0.5)+
  labs(x ="chaperons",
       y = "folding potential (%)") + coord_flip() +
  paper_figs_theme

ggsave("output/paper_figures/chap_fold_potential.pdf", g)
ggsave("~/Dropbox/Apps/Overleaf/Cancer_ecology/Fig/chap_fold_potential.pdf", g)

# realized niche:
folding_percent <- read.csv(file = "output/chap_realized_niche.csv", row.names = 1)

# reorder heat map by row sum and col sum:
chp_sum <- rowSums(folding_percent)
chap_ordered <- t(folding_percent[order(chp_sum,decreasing=T),])
cncr_sum <- rowSums(chap_ordered)
chap_ordered <- t(chap_ordered[order(cncr_sum,decreasing=T),])

pheatmap(chap_ordered, cluster_rows = F, cluster_cols = F,
         filename = "output/paper_figures/folding_percentage_nestedness_heatmap.pdf",
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100),
         angle_col = 45 )

# generalizm: 
folding_percent_all <- 
  read.csv(file = "output/chap_folding_percent_of_all.csv", row.names = 1)

# reorder heat map by row sum and col sum:
chp_sum <- rowSums(folding_percent_all)
chap_ordered <- t(folding_percent_all[order(chp_sum,decreasing=T),])
cncr_sum <- rowSums(chap_ordered)
chap_ordered <- t(chap_ordered[order(cncr_sum,decreasing=T),])

pheatmap(chap_ordered, cluster_rows = F, cluster_cols = F,
         filename = "output/paper_figures/folding_percentage_of_total_nestedness_heatmap.pdf",
         color = colorRampPalette(brewer.pal(n = 7, name = "YlGnBu"))(100),
         angle_col = 45)

# eigen values - rn:
ev_both <- read.csv("output/data/rn_vs_shuff_eigenvalues.csv")
ev_obs <- ev_both[1,1]
ev_shuffled <- ev_both[2:nrow(ev_both), ]

pp <- ggplot(as.data.frame(ev_shuffled), aes(x = ev_shuffled)) +
  geom_histogram(colour = "darkblue", aes(fill = ..count..), bins = 30) +
  scale_x_continuous(name = "Simulated eigenvalue") +
  scale_y_continuous(name = "Count") +
  #ggtitle("Histogram of simulated eigenvalue\nfor chaperones realized niche") +
  geom_vline(xintercept = ev_obs, linetype="dashed", color = "red", size=1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  paper_figs_theme
pp
ggsave("output/paper_figures/nestedness_shuffled_realized_niche.pdf",
       plot = pp, width = 3.2, height = 2.8, units = "in")

# eigen values - generalizm:
ev_both <- read.csv("output/data/generalizm_vs_shuff_eigenvalues.csv")
ev_all_obs <- ev_both[1,1]
ev_shuffled <- ev_both[2:nrow(ev_both), ]

p <- ggplot(as.data.frame(ev_shuffled), aes(x = ev_shuffled)) +
  geom_histogram(colour = "darkblue", aes(fill = ..count..), bins = 30) +
  scale_x_continuous(name = "Simulated eigenvalue") +
  scale_y_continuous(name = "Count") +
  #ggtitle("Histogram of simulated eigenvalue\nfor chaperones' generalism") +
  geom_vline(xintercept = ev_all_obs, linetype="dashed", color = "red", size=1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  paper_figs_theme
p
ggsave("output/paper_figures/nestedness_shuffled_generalism.pdf", 
       plot = p, width = 3.2, height = 2.8, units = "in")

# expression ------
chap_expr <- read.csv("output/chap_median_expressions.csv", row.names = 1)
chap_expr_long <- melt(t(chap_expr))
chap_expr_long$log <- log10(chap_expr_long$value)

# chap boxplot
g1 <- ggplot(chap_expr_long, aes(x=reorder(Var2, log, FUN = median, ),y=log))+
      geom_boxplot(outlier.colour="black", outlier.shape=16,
                   outlier.size=2, notch=FALSE) +
      labs(y="Log10 expression levels") +
      theme_minimal()+ 
      theme(axis.text.x=element_text(angle=45, hjust=1),
            axis.title.x = element_blank()) +
      paper_figs_theme
ggsave("output/paper_figures/chap_exp_boxplot.pdf",
       plot = g1, width = 5, height = 5)

# chap histogram
p<-ggplot(chap_expr_long, aes(x=log)) + 
  geom_histogram(bins=30, fill="coral", alpha=0.7) + 
  labs(x="Log10 expression levels")+
  paper_figs_theme 
ggsave("output/paper_figures/chap_exp_hist.pdf",
       plot = p, width = 5, height = 5)

# chap heatmap
# reorder heat map by row sum and col sum:
exp_median <- read.csv("output/chap_median_expressions.csv", row.names = 1)
chp_sum <- rowSums(exp_median)
chap_ordered <- t(exp_median[order(chp_sum,decreasing=T),])
cncr_sum <- rowSums(chap_ordered)
chap_ordered <- t(chap_ordered[order(cncr_sum,decreasing=T),])


pheatmap(log10(chap_ordered), cluster_rows = F, cluster_cols = F,
         color=colorRampPalette(c("white", "orange"))(100),
         filename = "output/Median_log10_chap_expression_heatmap.pdf",
         #main = "Log10 of Median Chaperon Expression", angle_col = 315,
         display_numbers = TRUE)



# protein expression
exp_med_prot <- read.csv("output/prot_median_expressions.csv", row.names = 1)
prot_exp <- melt(t(log10(exp_med_prot)))

g <- ggplot(prot_exp, aes(x=value, fill=Var2, color=Var2)) +
  geom_histogram(position="identity", alpha=0.5)+
  xlab("Log10 of protein expression") + ylab("Count")+
  paper_figs_theme
ggsave("output/figures/prot_med_exp_per_cancer.pdf", g)

# similarity ------

# similarity boxplot - cancer
all_simlrs <- read.csv("output/jaccard_values_per_cancer.csv")

# melt the data to a long format
mlt_sim <- as.data.frame(melt(all_simlrs))
g1 <- ggplot(mlt_sim, aes(x=variable,y=value))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  labs(x=element_blank(), y="Jaccard similarity index") +
  paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
g1
ggsave("output/paper_figures/cancer_jaccard_boxplot.pdf", g1)


# similarity per cancer (between chaps) + shuff
combine_dfs <- read_csv("output/data/cancer_jaccard_with_shuff.csv")

p3 <- ggplot(combine_dfs%>% group_by(kind), aes(x=value, fill=kind)) + 
  geom_histogram(aes(y = stat(density)),
                 alpha=0.5, position = 'identity',
                 bins=30) + 
  labs(x="Jaccard similarity index",
       y="Density", fill="Population") +
  paper_figs_theme
p3
ggsave("output/paper_figures/shuffled_jaccard_per_cancer.pdf", p3)


# similarity boxplot - chaps
all_simlrs <- read_csv("output/jaccard_values_per_chap.csv")
mlt_sim <- as.data.frame(melt(all_simlrs))
g1 <- ggplot(mlt_sim, aes(x=variable,y=value))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  labs(x=element_blank(), y="Jaccard similarity index") +
  paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
g1
ggsave("output/paper_figures/chap_jaccard_boxplot.pdf", g1)


# similarity per chaps (between cancers) + shuff
combine_dfs <- read_csv("output/data/chap_jaccard_with_shuff.csv")

p3 <- ggplot(combine_dfs%>% group_by(kind), aes(x=value, fill=kind)) + 
  geom_histogram(aes(y = stat(density)),
                 alpha=0.5, position = 'identity',
                 bins=30) + 
  labs(x="Jaccard similarity index",
       y="Density", fill="Population") +
  paper_figs_theme
p3
ggsave("output/paper_figures/shuffled_jaccard_per_chap.pdf")

# infomap ------
concluting_table <- read_csv('output/multilayer_relaxed_scan_20_trials.csv')

concluting_table %>% 
  filter(relax_param==0.15) %>% 
  filter(type=='chaperone') %>%
  group_by(symbol, module) %>% 
  summarise(n=n_distinct(cancer)) %>% 
  ggplot(aes(symbol, module, fill=n, label=n))+geom_tile(color='navy')+
  geom_text()+ xlab("Chaperone") + ylab("Module number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(angle = 45, hjust=1))
ggsave("output/paper_figures/multilayer_moduls_per_chap.pdf")

# mix correlations -----
# realized niche vs similarity per chap
all_stats <- read.csv("output/data/rn_similarity_expr_stats.csv", row.names = 1)

# ggplot
g <- ggplot(all_stats, aes(x=sim_med, y=fold_med))+
      geom_point(size=2.5) +
      geom_errorbar(aes(ymax = fold_q3, ymin = fold_q1), width = 0.005, alpha=.4) + 
      geom_errorbarh(aes(xmax = sim_q3, xmin = sim_q1), height = 0.005, alpha=.4) + 
      labs(x="Similarity", y="Realized Niche (%)") + 
      paper_figs_theme
ggsave("output/paper_figures/similarity_realized_scatter.pdf", g)


# rn vs chap expression
tbl_all_2 <- as.tibble(read.csv("output/data/chap_rn_and_log_exp.csv"))

g <- ggplot(tbl_all_2, aes(x=expression, y=fold, color=chap))+
  geom_point(size=2)+
  #ggtitle("Realized Niche over Log10 Median Expression levels")+
  ylab("Realized Niche (%)")+
  xlab("Log10 median expression level") + 
  labs(color="Chaperon") +
  paper_figs_theme
ggsave("output/paper_figures/log10_rn_vs_exp.pdf", g)

# rn vs chap expression per chap
g <- ggplot(tbl_all_2, aes(x=expression, y=fold, color=chap)) +
  geom_point() + 
  ylab("Realized Niche (%)") + xlab("Log10 median expression level")+
  #ggtitle("Realized Niche over Median Expression levels") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.5, "lines"),
    strip.text.x = element_text(size = 10),
  ) +
  paper_figs_theme_no_legend + 
  facet_wrap(~ chap)
ggsave("output/paper_figures/log10_rn_vs_exp_per_chap.pdf", g)
