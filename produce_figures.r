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
library(igraph)


source("functions.r")

# ----- consts -----
output_folder <- "output/paper_figures/"
drop_box <- "~/Dropbox/Apps/Overleaf/Cancer_ecology/Nature Comm/Revision/"

# ----- Make figures -----

## potential and nestedness --------
# potential bars:
puf <- read.csv("output/folding_potential.csv")

pot <- ggplot(puf, aes(x=reorder(name, -potential), y=potential)) +
  geom_bar(stat="identity", fill="steelblue", width=0.5)+
  labs(x ="chaperons",
       y = "folding potential (%)") + coord_flip() +
  paper_figs_theme
pot

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

## expression ------
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

# chap histogram
e2<-ggplot(chap_expr_long, aes(x=log)) + 
  geom_histogram(bins=30, fill="coral", alpha=0.7) + 
  labs(x="Log10 expression levels")+
  paper_figs_theme 
e2

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

## similarity ------

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

# similarity per cancer (between chaps) + shuff
combine_dfs <- read_csv("output/data/cancer_jaccard_with_shuff.csv")

sim2 <- ggplot(combine_dfs%>% group_by(kind), aes(x=value, fill=kind)) + 
  geom_histogram(aes(y = stat(density)),
                 alpha=0.5, position = 'identity',
                 bins=30) + 
  labs(x="Jaccard similarity index",
       y="Density", fill="Population") +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  paper_figs_theme + 
  theme(legend.position = c(0.82,0.88),
        legend.title = element_blank())
sim2


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

## infomap ------
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

## mix correlations -----
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
tbl_all_2 <- as_tibble(read.csv("output/data/chap_rn_and_log_exp.csv"))

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


## robustness -----
# stb1 - network collapse by removal
colps <- read_csv("output/data/collapse_data_by_module_for_paper.csv")

stb1 <- colps %>%
  mutate(cancer=factor(cancer, levels=cancer_nestedness_order)) %>%
  mutate(y_txt= case_when(run_name=='high_to_low' ~ 0.50, 
                          run_name=='by_module_21' ~ 0.35,
                          run_name=='by_module_12' ~  0.19,
                          run_name=='random' ~ 0.05)) %>%
  mutate(labels=case_when(run_name=='high_to_low' ~ 'High to low', 
                          run_name=='by_module_21' ~ 'By module 2 then 1',
                          run_name=='by_module_12' ~  'By module 1 then 2',
                          run_name=='random' ~ 'Random removal')) %>%
  ggplot(aes(prop_removed, prop_remain, color=labels))+
  geom_point(size=2)+
  geom_line(size=1)+
  labs(x="% of chaperons removed", 
       y="% of proteins remained",
       color="Removal type") +
  # rearrange the legend items
  scale_color_discrete(breaks=c('High to low', 'By module 2 then 1', 'By module 1 then 2', 'Random removal')) +
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


## Affirm chaperon co-expression sets-----
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

# pvalues taken from interaction_evidance.r
dat_text <- data.frame(label = c("p=0.034", "p<0.01", "p<0.01"),
                       Chap = c("CLPP", "HSPD1", "TRAP1"))

# plot the distributions
ppraff <- ggplot(long_dist, aes(x=affirm_percentage)) + 
  geom_histogram() + 
  facet_grid(Chap~.) +
  xlab("% of interactions affirmed") + ylab("Count") +
  geom_vline(data = obs, color="red",
             aes(xintercept = affirm_percentage)) + 
          paper_figs_theme_no_legend +
  theme(axis.title.y = element_text(margin=margin(r=10))) + 
  geom_text(data = dat_text,
             mapping = aes(x = -Inf, y = -Inf, label = label, size=2),
             hjust   = -0.1, vjust   = -8.5)
ppraff


## Inter-cancer relationships --------
# plot jaccard data
cj_data <- read.csv(file = "output/data/cancer_edges_jaccard.csv")
cncr_order <- cj_data[1:12, 2]

cj <- ggplot(cj_data, aes(x=cancer1, y=fct_relevel(cancer2, rev(cncr_order)), fill=jaccard)) + 
  geom_tile() + 
  scale_fill_gradient(low = "lightyellow", high = "#009640", limits=c(0,1)) +
  paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank())
cj


# plot xei's work to have a standard look in the paper:
# plot link prediction data
lp_data <- read.csv(file = "output/data/link_prediction_plot_ready.csv")

lp <- ggplot(lp_data, aes(x=From, y=fct_relevel(To, rev(cncr_order)), fill=AUC)) + 
  geom_tile() + 
  scale_fill_gradient2(midpoint = 0.5, low = "#832424", mid="white", 
                       high = "#3A3A98", limits=c(0,1)) +
  paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank())
lp

# physical node membership - probabilities (TODO maybe to be removed)
mem_data <- read.csv(file = "output/data/chap_membership_matrix.csv") %>%
  select(Chaperon=Name, X1, X2)
mem_data <- melt(mem_data)
mem_data <- mem_data %>% mutate(module=case_when(variable=="X1" ~ 1,
                                                 variable=="X2" ~ 2))

mprb <- ggplot(mem_data, aes(y=fct_relevel(Chaperon, rev(chap_module_order)), 
                             x=factor(module), fill=value)) + 
          geom_tile() + 
          scale_fill_gradient(low = "lightyellow", high = "navyblue") +
          paper_figs_theme + labs(x="module ID", fill = "memb.\nprob.") +
          theme(axis.title.y = element_blank(),
                panel.border = element_blank())

# physical node membership - hard membership
bi_file <- "output/data/bipartite_membership.csv"
membership_vercors <- read.table(bi_file, sep=",", header=TRUE,
                                 stringsAsFactors=FALSE, quote="", fill=FALSE)
dta <- membership_vercors[,c(1,3)] %>% select(Chaperon=X, module=community_2)
dta <- dcast(dta, Chaperon ~ module,fill = 0)
dta <- melt(dta) %>% mutate(module=case_when(value>0 ~ 1,
                                             value==0 ~ 0))
mprb <- ggplot(dta, aes(y=fct_relevel(Chaperon, rev(chap_module_order)), 
                        x=factor(variable), fill=module)) + 
  geom_tile() + 
  scale_fill_gradient(low = "lightyellow", high = "#CC79A7") +
  paper_figs_theme_no_legend + 
  labs(x="module ID", fill = "memb.\nprob.") +
  theme(axis.title.y = element_blank(),
        panel.border = element_blank())
mprb


# 



# Print figures -----
# print the figures according to needed in the paper, in pdf format.
# note: plot.margin order of element is: t -> r -> b -> l


## main paper ----
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


## supplementary ----
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


## revision ----
# fig affirmation - dbaff + ppraff
pdf(paste(drop_box,'SI_affirmation.pdf', sep = ""), 10, 5)
plot_grid(dbaff + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")), 
          ppraff + theme(plot.margin = unit(c(0.25,0.25,0.5,0.5), "cm")), 
          labels = c('(A)', '(B)'), 
          rel_widths = c(1,1))
dev.off()

# cancer comparison- remastering figure 4 - sim2 + mprb + cj + lp
pdf(paste(drop_box,'niche_separation.pdf', sep = ""), 10, 8)
plot_grid(sim2 + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")), 
          mprb + theme(plot.margin = unit(c(0.40,2.5,0.5,0.7), "cm")),
          lp + theme(plot.margin = unit(c(0.75,0.25,0.5,0.5), "cm")), 
          cj + theme(plot.margin = unit(c(0.75,0.25,0.5,0.5), "cm")), 
          labels = c('(A)', '(B)', '(C)', "(D)"), 
          rel_widths = c(1,1))
dev.off()


# Side Figures: -----------------
# plotting the cancer as separate networks
# ----- plot cancer layers
# metadata:
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
clients_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]

# read the layers:
aj_file <- "output/data/adjacency_edgelist.csv"
aj_data <- read.csv(aj_file, header=TRUE, stringsAsFactors=FALSE, fill=FALSE)

# function for one cancer network
plot_cancer_net <- function(cancer_name, nets, meta_chap, meta_client){
  # prepare cancer:
  net <- nets[, c("node1","node2", cancer_name)]
  names(net) <- c("Chap", "Prot", "is_signfcnt")
  net <- net %>% filter(is_signfcnt > 0)
  
  #create vertices:
  vert_ch <- meta_chap %>% add_column(shape = "rectangle", 
                                      color = "green",
                                      label= meta_chap$Symbol,
                                      size = 15) %>% 
    select(Symbol, shape, color, size, label)
  vert_pr <- meta_client %>% add_column(shape = "circle", 
                                        color = "orange",
                                        label= "",
                                        size = 3) %>% 
    select(Symbol, shape, color, size, label)
  vert_pr <- vert_pr[vert_pr$Symbol %in% net$Prot,]
  vert <- rbind(vert_ch, vert_pr)
  
  # plot the network:
  rrr <- igraph::graph_from_data_frame(d = net, vertices = vert ,directed = FALSE)
  plot.igraph(rrr,  axes = FALSE, #vertex.frame.color = NA,
              vertex.label = vert$label, vertex.label.cex=0.5,
              vertex.size = vert$size,
              vertex.color = vert$color,
              vertex.shape = vert$shape,
              main = cancer_name,
              layout=layout_with_mds)
}

#for all cancers
pdf(paste(drop_box,'SI_networks_as_graphs.pdf', sep = ""))
for (c in sort(cancer_sample_size_order)) {
  plot_cancer_net(c, aj_data, chaps_meta, clients_meta)
}
dev.off()



# Conjugation rate was estimated as:
#   Ψ * ln(1 + (T / R)(N / D))*(1 / (N – N0)
