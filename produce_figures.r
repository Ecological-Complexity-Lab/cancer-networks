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
library(cowplot)
library(tidyverse)

source("functions.r")

# ----- consts -----
output_folder <- "output/paper_figures/"
drop_box <- "~/Dropbox/Apps/Overleaf/Cancer_ecology/Fig/"

# ----- Make figures -----

# potential and nestedness --------
# potential bars:
puf <- read.csv("output/folding_potential.csv")

pot <- ggplot(puf, aes(x=reorder(name, -potential), y=potential)) +
  geom_bar(stat="identity", fill="steelblue", width=0.5)+
  labs(x ="chaperons",
       y = "folding potential (%)") + coord_flip() +
  paper_figs_theme
pot
#ggsave("output/paper_figures/chap_fold_potential.pdf", pot)

# realized niche:
folding_percent <- read.csv(file = "output/chap_realized_niche.csv", row.names = 1)

# reorder heat map by row sum and col sum:
data1 <- melt(as.matrix(folding_percent))
colnames(data1) <- c("chap", "cancer", "value")
chap_order1 <- rownames(folding_percent)[order(rowSums(folding_percent), decreasing = FALSE)]
cancer_order1 <- colnames(folding_percent)[order(colSums(folding_percent), decreasing = TRUE)]

# ggplot Heatmap 
hm1 <- ggplot(data1, aes(x=fct_relevel(cancer, cancer_order1), 
                 y=fct_relevel(chap, chap_order1), fill=value)) + 
        geom_tile() + 
        scale_fill_gradient(low = "lightyellow", high = "maroon") +
        geom_rect(aes(xmin = 0.5, xmax = 4.5, ymin = 11.45, ymax = 15.5), 
                  color="black", fill = NA, size = 1.5) +  
        paper_figs_theme + 
        theme(axis.text.x=element_text(angle=45, hjust=1),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.title = element_blank(),
              panel.border = element_blank())
hm1

# generalizm: 
folding_percent_all <- 
  read.csv(file = "output/chap_folding_percent_of_all.csv", row.names = 1)

# reorder heat map by row sum and col sum:
data2 <- melt(as.matrix(folding_percent_all))
colnames(data2) <- c("chap", "cancer", "value")
chap_order2 <- rownames(folding_percent_all)[order(rowSums(folding_percent_all), decreasing = FALSE)]
cancer_order2 <- colnames(folding_percent_all)[order(colSums(folding_percent_all), decreasing = TRUE)]

# ggplot Heatmap 
hm2 <- ggplot(data2, aes(x=fct_relevel(cancer, cancer_order2), 
                        y=fct_relevel(chap, chap_order2), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient(low = "lightyellow", high = "navyblue") +
  # scale_fill_gradient(low = "lightyellow", high = "maroon") +
  geom_rect(aes(xmin = 0.5, xmax = 4.5, ymin = 11.45, ymax = 15.5), 
            color="black", fill = NA ,size = 1.5) +  
  paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(), 
        panel.border = element_blank())
hm2


# eigen values - rn:
ev_both <- read.csv("output/data/rn_vs_shuff_eigenvalues.csv")
ev_obs <- ev_both[1,1]
ev_shuffled <- ev_both[2:nrow(ev_both), ]

ev1 <- ggplot(as.data.frame(ev_shuffled), aes(x = ev_shuffled)) +
  geom_histogram(colour = "darkblue", aes(fill = ..count..), bins = 30) +
  scale_x_continuous(name = "Largest eigenvalue") +
  scale_y_continuous(name = "Count") +
  #ggtitle("Histogram of simulated eigenvalue\nfor chaperones realized niche") +
  geom_vline(xintercept = ev_obs, linetype="dashed", color = "red", size=1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  paper_figs_theme
ev1
#ggsave("output/paper_figures/nestedness_shuffled_realized_niche.pdf",
#       plot = ev1, width = 3.2, height = 2.8, units = "in")

# eigen values - generalizm:
ev_both <- read.csv("output/data/generalizm_vs_shuff_eigenvalues.csv")
ev_all_obs <- ev_both[1,1]
ev_shuffled <- ev_both[2:nrow(ev_both), ]

ev2 <- ggplot(as.data.frame(ev_shuffled), aes(x = ev_shuffled)) +
  geom_histogram(colour = "darkblue", aes(fill = ..count..), bins = 30) +
  scale_x_continuous(name = "Largest eigenvalue") +
  scale_y_continuous(name = "Count") +
  #ggtitle("Histogram of simulated eigenvalue\nfor chaperones' generalism") +
  geom_vline(xintercept = ev_all_obs, linetype="dashed", color = "red", size=1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  paper_figs_theme
ev2
#ggsave("output/paper_figures/nestedness_shuffled_generalism.pdf", 
#       plot = ev2, width = 3.2, height = 2.8, units = "in")

# expression ------
chap_expr <- read.csv("output/chap_median_expressions.csv", row.names = 1)
chap_expr_long <- melt(t(chap_expr))
chap_expr_long$log <- log10(chap_expr_long$value)

# chap boxplot
e1 <- ggplot(chap_expr_long, aes(x=reorder(Var2, log, FUN = median, ),y=log))+
      geom_boxplot(outlier.colour="black", outlier.shape=16,
                   outlier.size=2, notch=FALSE) +
      labs(y="Log10 expression levels") +
      theme_minimal() + paper_figs_theme +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            axis.title.x = element_blank())
 e1     
#ggsave("output/paper_figures/chap_exp_boxplot.pdf",
#       plot = e1, width = 5, height = 5)

# chap histogram
e2<-ggplot(chap_expr_long, aes(x=log)) + 
  geom_histogram(bins=30, fill="coral", alpha=0.7) + 
  labs(x="Log10 expression levels")+
  paper_figs_theme 
e2
#ggsave("output/paper_figures/chap_exp_hist.pdf",
#       plot = e2, width = 5, height = 5)

# chap heatmap
# reorder heat map by row sum and col sum:
exp_median <- log10(read.csv("output/chap_median_expressions.csv", row.names = 1))

# reorder heat map by row sum and col sum:
data <- melt(as.matrix(exp_median))
colnames(data) <- c("chap", "cancer", "value")
chap_order <- rownames(exp_median)[order(rowSums(exp_median), decreasing = FALSE)]
cancer_order <- colnames(exp_median)[order(colSums(exp_median), decreasing = TRUE)]

# ggplot Heatmap 
e3 <- ggplot(data, aes(x=fct_relevel(cancer, cancer_order), 
                        y=fct_relevel(chap, chap_order), fill=value)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "orange") +
  paper_figs_theme + 
  labs(fill="log10") + 
  geom_text(aes(label=round(value, digits = 2)),size=3) +  
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank())
e3

# protein expression
exp_med_prot <- read.csv("output/prot_median_expressions.csv", row.names = 1)
prot_exp <- melt(t(log10(exp_med_prot)))

e4 <- ggplot(prot_exp, aes(x=value, fill=Var2, color=Var2)) +
  geom_histogram(position="identity", alpha=0.6)+
  labs(x="Log10 of protein expression", y="Count") +
  paper_figs_theme + 
  theme(legend.title = element_blank())
e4
#ggsave("output/figures/prot_med_exp_per_cancer.pdf", e4)

# similarity ------

# similarity boxplot - cancer
all_simlrs <- read.csv("output/jaccard_values_per_cancer.csv")

# melt the data to a long format
mlt_sim <- as.data.frame(melt(all_simlrs))
sim1 <- ggplot(mlt_sim, aes(x=variable,y=value))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  labs(x=element_blank(), y="Jaccard similarity index") +
  paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
sim1
#ggsave("output/paper_figures/cancer_jaccard_boxplot.pdf", sim1)


# similarity per cancer (between chaps) + shuff
combine_dfs <- read_csv("output/data/cancer_jaccard_with_shuff.csv")

sim2 <- ggplot(combine_dfs%>% group_by(kind), aes(x=value, fill=kind)) + 
  geom_histogram(aes(y = stat(density)),
                 alpha=0.5, position = 'identity',
                 bins=30) + 
  labs(x="Jaccard similarity index",
       y="Density", fill="Population") +
  paper_figs_theme + 
  theme(legend.position = c(0.82,0.87))
sim2
#ggsave("output/paper_figures/shuffled_jaccard_per_cancer.pdf", sim2)


# similarity boxplot - chaps
all_simlrs <- read_csv("output/jaccard_values_per_chap.csv")
mlt_sim <- as.data.frame(melt(all_simlrs))
sim3 <- ggplot(mlt_sim, aes(x=variable,y=value))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  labs(x=element_blank(), y="Jaccard similarity index") +
  paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
sim3
#ggsave("output/paper_figures/chap_jaccard_boxplot.pdf", sim3)


# similarity per chaps (between cancers) + shuff
combine_dfs <- read_csv("output/data/chap_jaccard_with_shuff.csv")
combine_dfs %<>% 
  mutate(new_kind=case_when(kind=='obs' ~ 'Observed',
                            kind=='shuff' ~ 'Shuffled')) %>% 
  group_by(kind)

sim4 <- ggplot(combine_dfs, aes(x=value, fill=new_kind)) + 
  geom_histogram(aes(y = stat(density)),
                 alpha=0.5, position = 'identity',
                 bins=30) + 
  labs(x="Jaccard similarity index",
       y="Density", fill="Population") +
  paper_figs_theme + 
  theme(legend.position = c(0.82,0.87))
sim4
#ggsave("output/paper_figures/shuffled_jaccard_per_chap.pdf")

# infomap ------
concluting_table <- read_csv('output/multilayer_relaxed_scan_20_trials.csv')

inf <- concluting_table %>% 
        filter(relax_param==0.15) %>% 
        filter(type=='chaperone') %>%
        group_by(symbol, module) %>% 
        summarise(n=n_distinct(cancer)) %>% 
        ggplot(aes(x=fct_relevel(symbol, chap_module_order), 
                   module, fill=n, label=n))+geom_tile(color='navy') +
        geom_text()+ 
        labs(x=element_blank(), y="Module number", fill="Cancer\nnumber")+
        paper_figs_theme + 
        theme(axis.text.x = element_text(angle = 45, hjust=1))
#ggsave("output/paper_figures/multilayer_moduls_per_chap.pdf")

# mix correlations -----
# realized niche vs similarity per chap
all_stats <- read.csv("output/data/rn_similarity_expr_stats.csv", row.names = 1)

# ggplot
srn <- ggplot(all_stats, aes(x=sim_med, y=fold_med))+
      geom_point(size=2.5) +
      geom_errorbar(aes(ymax = fold_q3, ymin = fold_q1), width = 0.005, alpha=.4) + 
      geom_errorbarh(aes(xmax = sim_q3, xmin = sim_q1), height = 0.005, alpha=.4) + 
      labs(x="Similarity", y="Realized Niche (%)") + 
      paper_figs_theme
#ggsave("output/paper_figures/similarity_realized_scatter.pdf", srn)


# rn vs chap expression
tbl_all_2 <- as.tibble(read.csv("output/data/chap_rn_and_log_exp.csv"))

# rn vs chap expression per chap
ern2 <- ggplot(tbl_all_2, aes(x=expression, y=fold, color=chap)) +
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
#ggsave("output/paper_figures/log10_rn_vs_exp_per_chap.pdf", ern2)


# robustness -----
# TODO ?? add cancer lots and..? network weakening index? 
# (something to be per chap)


# ----- print all -----
# fig 1 - hm1 + hm2
pdf(paste(drop_box,'nestedness.pdf', sep = ""), 10, 4)
plot_grid(hm1 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          hm2 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()


# fig 2 - e1 + ern2
pdf(paste(drop_box,'expression.pdf', sep = ""), 10, 6)
plot_grid(e1, ern2 + theme(plot.margin = unit(c(0.2,1,1,1), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(0.4,0.6))
dev.off()

# fig 3 - sim2 + inf
pdf(paste(drop_box,'niche_separation.pdf', sep = ""), 10, 5)
plot_grid(sim2 + theme(plot.margin = unit(c(0.2,0.2,1.1,0.5), "cm")), inf, 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()

# fig 4 - sim4 + srn
pdf(paste(drop_box,'cross_cancer_similarity.pdf', sep = ""), 10, 5)
plot_grid(sim4, srn, 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()



# fig S1 - pot
pdf(paste(drop_box,'chap_fold_potential.pdf', sep = ""), 10, 6)
pot
dev.off()

# fig S2 - ev2 + ev1 
pdf(paste(drop_box,'nestedness_shuff.pdf', sep = ""), 10, 4)
plot_grid(ev2 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          ev1 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()

# fig S3 - e3 + e4
pdf(paste(drop_box,'expression_distribution.pdf', sep = ""), 10, 4)
plot_grid(e4 + theme(plot.margin = unit(c(0.2,0.2,0.75,0.2), "cm")), 
          e3 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(0.4,0.6))
dev.off()

# fig S4 - sim1 + sim3
pdf(paste(drop_box,'similarity_boxplots.pdf', sep = ""), 10, 5)
plot_grid(sim1 + theme(plot.margin = unit(c(0.2,0.2,0.75,0.5), "cm")), 
          sim3 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()
