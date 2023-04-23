# ------- produce_figures.r ------
# this script is meant to streamline the making of figures for the paper.
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
library(ggpmisc)
library(ggpubr)

source("functions.r")

# ----- consts -----
output_folder <- "output/paper_figures/"
drop_box <- "~/Dropbox/Apps/Overleaf/Cancer_ecology/"

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
  theme(legend.position = c(0.82,0.9),
        legend.title = element_blank())
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
  theme(legend.position = c(0.82,0.9),
        legend.title = element_blank())
sim4
#ggsave("output/paper_figures/shuffled_jaccard_per_chap.pdf")

# infomap ------
concluting_table <- read_csv('output/multilayer_relaxed_scan_20_trials.csv')
temp <- concluting_table %>% 
        filter(relax_param==0.15) %>% 
        filter(type=='chaperone') %>%
        group_by(symbol, module) %>% 
        summarise(n=n_distinct(cancer)) %>% 
        mutate(colour=case_when(module==1 ~ 'lightblue',
                                module==2 ~ 'pink',
                                module==3 ~ '#D1F3C5',
                                module>3  ~ 'lightgrey'))
inf <-  temp%>%
        ggplot(aes(x=fct_relevel(symbol, chap_module_order), 
                   y=factor(module), label=n)) + 
        geom_tile(color='navy', fill = temp$colour) + geom_text()+ 
        labs(x=element_blank(), y="Module number", fill="Cancer\nnumber")+
        paper_figs_theme_no_legend + 
        theme(axis.text.x = element_text(angle = 45, hjust=1))

# mix correlations -----
# realized niche vs similarity per chap
all_stats <- read.csv("output/data/rn_similarity_expr_stats.csv", row.names = 1)

# ggplot
srn <- ggplot(all_stats, aes(x=sim_med, y=fold_med, label=rownames(all_stats)))+
  geom_point(size=2.5) +
  geom_errorbar(aes(ymax = fold_q3, ymin = fold_q1), width = 0.005, alpha=.3) + 
  geom_errorbarh(aes(xmax = sim_q3, xmin = sim_q1), height = 0.005, alpha=.3) + 
  labs(x="Similarity", y="Realized Niche (%)") + 
  annotate("text", x=0.12, y=0.62, label= "R=0.64\nP=0.0097", size = 5) +
  annotate("text", x=0.395, y=0.46, label= "HSPE1", size = 4) + 
  annotate("text", x=0.48, y=0.5, label= "CLPP", size = 4) + 
  paper_figs_theme
srn

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
  theme(axis.text = element_text(size=10, color='black'))+ 
  facet_wrap(~ chap)


# robustness -----
# stb1 - network collapse by removal
colps <- read_csv("output/data/collapse_data_by_module_for_paper.csv")

stb1 <- colps %>%
  mutate(cancer=factor(cancer, levels=cancer_nestedness_order)) %>%
  mutate(y_txt= case_when(run_name=='by_module_123' ~  0.50,
                          run_name=='by_module_213' ~ 0.35,
                          run_name=='high_to_low' ~ 0.19, 
                          run_name=='random' ~ 0.05)) %>%
  ggplot(aes(prop_removed, prop_remain, color=run_name))+
  geom_point(size=2)+
  geom_line(size=1)+
  labs(x="% of chaperons removed", 
       y="% of proteins remained",
       color="Removal type") +
  # Add the R to each cancer and each removal plot in the white-space using: 
  geom_text(aes(x=0.09, y=y_txt, label=round(under_curve, digits = 3)), stat = "unique") +
  
  facet_wrap(vars(cancer), ncol = 4) +
  paper_figs_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# stb2 - correlate each removal to realized niche 
sms <- read_csv("output/data/corr_stbility_vs_rn_by_module_for_paper.csv")

my.f <- y ~ x
stb2 <- sms %>%
  ggplot(aes(x=rn_sum, y=under_curve, color=run_name)) + 
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  labs(x="Cancer realized niche sum", y="Area under extinction curve") +
  facet_wrap(vars(run_name), nrow = 2, ncol = 2) + paper_figs_theme_no_legend + 
  stat_cor(aes(label = ..r.label..), method = "spearman", 
           label.y = 0.665, label.x = 6.2, size = 3) + 
  stat_cor(aes(label = ..p.label..), method = "spearman", 
           label.y = 0.65, label.x = 6.2, size = 3)
stb2


# fig 6 - stb2
pdf(paste(drop_box,'correlation.pdf', sep = ""), 5, 5)
stb2
dev.off()

# Affirm chaperon co-expression sets-----
# using STRING DB:
perc_tibble <- as_tibble(read.csv(file = "output/data/STRING_affirm_percentage.csv")) %>% 
  filter(evidence == "experiments")
perc_tibble$type <- "obs"

# read random data
perc <- read.csv(file = "HPC/DATA/STRINGdb_rand_percentage.csv", header = F) # string
names(perc) <- c("Chap", "evidence", "value", "type")
exp_perc <- perc %>% filter(evidence == "experiments")

# plot random distribution per chap + observed as a point.
all_exp <- rbind(perc_tibble, exp_perc)

dbaff <- ggplot() +
  geom_boxplot(data=all_exp, aes(x=factor(Chap), y=value)) + 
  geom_point(data=all_exp[all_exp$type == "obs",], 
             aes(x=factor(Chap), y=value), color="red", size = 3) + 
  paper_figs_theme + ylab("% of interactions affirmed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.title.x=element_blank())

# using papers:
pprs_data <- read.csv(file = "output/data/affirmation_sims_papers.csv")
long_dist <- pprs_data %>% filter(type=="shuff")
obs <- pprs_data %>% filter(type=="obs")

# plot the distributions
ppraff <- ggplot(long_dist, aes(x=affirm_percentage)) + 
  geom_histogram() + facet_grid(Chap~.) +
  xlab("% of interactions affirmed") +
  geom_vline(data = obs, color="red",
             aes(xintercept = affirm_percentage)) + paper_figs_theme



# ----- print all -----
# fig 1 - hm1 + hm2
pdf(paste(drop_box,'nestedness.pdf', sep = ""), 5, 4)
hm1
dev.off()

# fig 2 - e1 + ern2
pdf(paste(drop_box,'expression.pdf', sep = ""), 10, 6)
plot_grid(e1, ern2 + theme(plot.margin = unit(c(0.2,1,1,1), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(0.4,0.6))
dev.off()

# fig 3 - sim4 + srn
pdf(paste(drop_box,'similarity.pdf', sep = ""), 10, 5)
plot_grid(sim4, srn, 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()

# fig 4 - sim2 + inf
pdf(paste(drop_box,'niche_separation.pdf', sep = ""), 10, 5)
plot_grid(sim2 + theme(plot.margin = unit(c(0.2,0.2,1.1,0.5), "cm")), inf, 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()



# fig 5 - stb1
pdf(paste(drop_box,'robustness.pdf', sep = ""), 10, 5)
stb1
dev.off()

# fig 6 - stb2
pdf(paste(drop_box,'correlation.pdf', sep = ""), 5, 5)
stb2
dev.off()



# fig S1 - pot
pdf(paste(drop_box,'SI_specialization.pdf', sep = ""), 10, 4)
plot_grid(pot + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
           hm2 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
           labels = c('(A)', '(B)'), 
           rel_widths = c(1,1))
dev.off()

# fig S2 - ev2 + ev1 
pdf(paste(drop_box,'SI_nestedness.pdf', sep = ""), 10, 4)
plot_grid(ev2 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          ev1 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()

# fig S3 - e3 + e4
pdf(paste(drop_box,'SI_expression.pdf', sep = ""), 10, 4)
plot_grid(e4 + theme(plot.margin = unit(c(0.2,0.2,0.75,0.2), "cm")), 
          e3 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(0.4,0.6))
dev.off()

# fig S4 - sim3 + sim1
pdf(paste(drop_box,'SI_similarity.pdf', sep = ""), 10, 5)
plot_grid(sim3 + theme(plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm")), 
          sim1 + theme(plot.margin = unit(c(0.2,0.2,0.75,0.5), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()


# revision figures:
# fig affirmation - dbaff + ppraff
pdf(paste(drop_box,'SI_affirmation.pdf', sep = ""), 10, 5)
plot_grid(dbaff + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")), 
          ppraff + theme(plot.margin = unit(c(0.25,0.25,0.5,0.5), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()

