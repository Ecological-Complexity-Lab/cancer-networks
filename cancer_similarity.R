#------------------------------
# use similarity indexes to find the differences 
# by comparing between the chaperons in each cancer
#------------------------------

#-------- includes --------
library(readxl)
library(vegan)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(igraph)

source("functions.r")

#-------- functions --------
long_format_from_dist_triangle <- function(similr, cancer_name) {
  # this function gets a distances triangle and turns it into a long format
  x <- as.matrix(similr)
  r_ind <- cbind(as.data.frame(rownames(x)), 1:15)
  names(r_ind) <- c("chap_row", "row")
  c_ind <- cbind(as.data.frame(rownames(x)), 1:15)
  names(c_ind) <- c("chap_col", "col")
  
  ind <- which(upper.tri(x, diag = FALSE), arr.ind = TRUE)
  sims_cancer <- as.data.frame(cbind(ind, value=x[ind])) %>% 
    left_join(y=r_ind, by=c("row")) %>%
    left_join(y=c_ind, by=c("col")) %>% 
    select(row=chap_row, col=chap_col, similarity=value) %>% 
    add_column(cancer=cancer_name)
  return(sims_cancer)
}

#only one cancer at a time
print_cancer_igraph <- function(cancer_name, c_att, inters) {
  # prepare vertices and edges
  intr_vals <- inters %>% filter(cancer == cancer_name) %>% select(from, to, prots)
  cancertail <- in_tail %>% filter(cancer == cancer_name) %>% 
    left_join(intr_vals, by=c("row"="from","col"="to")) %>% 
    left_join(intr_vals, by=c("row"="to","col"="from"))
  cancertail[is.na(cancertail)] <- 0
  cancertail <- cancertail %>% mutate(prots=prots.x+prots.y) %>% 
    select(chap1 = row, chap2 = col, cancer, prots)
  cancertail
  
  # plot the graph
  cancer_igraph <- graph_from_data_frame(d = cancertail, vertices = c_att, directed = FALSE)
  plot.igraph(cancer_igraph,  axes = FALSE, vertex.frame.color = NA,
              #vertex.label = V(g)$name, vertex.label.color = "gray20",
              vertex.size = 40, vertex.size2 = 30,
              vertex.color = chap_attrib$colour, #"gray90", vertex.frame.color = "gray20",
              vertex.shape = chap_attrib$shape,
              #edge.arrow.size=0.5, edge.color=col, 
              edge.width = cancertail$prots / 30,
              #edge.curved = T,
              main = cancer_name,
              margin = c(-0.05,-0.05,-0.07,-0.1),
              layout = as.matrix(chap_attrib[2:3]))
}


#-------- load the networks from an excel file, and calculate similarity --------
excel_path <- "HPC/binari_validated_corrs.xlsx"

# do the convert for every cancer
sheet_names <- excel_sheets(excel_path)

networks <- list()
sim_triangles <- list()
all_simlrs <- matrix(0, nrow = 105, ncol = 0)
all_simlrs_long <- NULL
for (name in sheet_names) {
  x <- as.data.frame(read_excel(excel_path, sheet = name, col_names = TRUE))
  x2 <- x[,-1]
  rownames(x2) <- x[,1]
  networks[[name]] <- x2
  
  # generate similarity between the chaperons in this cancer
  res <- vegdist(x2, method="jaccard")
  simlr <- 1-res
  
  all_simlrs <- cbind(all_simlrs, simlr)
  sim_triangles[[name]] <- simlr
  
  long <- long_format_from_dist_triangle(simlr, name)
  all_simlrs_long <- rbind(all_simlrs_long, long)
}

colnames(all_simlrs) <- sheet_names
write.csv(all_simlrs, file = "output/jaccard_values_per_cancer.csv", row.names = FALSE)
write.csv(all_simlrs_long, file = "output/jaccard_values_per_cancer_long_format.csv", row.names = FALSE)

#-------- visualize the results --------
all_simlrs <- read.csv("output/jaccard_values_per_cancer.csv")

# melt the data to a long format
mlt_sim <- as.data.frame(melt(all_simlrs))
g1 <- ggplot(mlt_sim, aes(x=variable,y=value))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  xlab("") + ylab("Similarity index")
g1

# distribute similarities
p<-ggplot(mlt_sim, aes(x=value)) + 
  geom_histogram(bins=40, color="black") + 
  labs(x="Jaccard similarity index")
p

# distribute similarities by cancer
p2<-ggplot(mlt_sim, aes(x=value, fill= variable)) + 
  geom_histogram(bins=20, color="black", position="identity", alpha=0.5) + 
  labs(x="Jaccard similarity index", fill = "Cancer")
p2

ggsave("output/paper_figures/cancer_jaccard_boxplot.pdf", g1)

#-------- investigate who is in the upper tail --------
all_simlrs_long <- read.csv("output/jaccard_values_per_cancer_long_format.csv")
in_tail <- all_simlrs_long %>% filter(similarity > 0.4)


chap_attrib <- read.csv("output/data/chap_attributes.csv") %>%
                arrange(v_order) %>% 
                mutate(shape=case_when(function.=='protease' ~ 'rectangle',
                                       function.=='folding' ~ 'circle')) %>%
                mutate(colour=case_when(module==1 ~ 'lightblue',
                                        module==2 ~ 'pink',
                                        module==3 ~ 'lightgreen'))

intersections_df <- read.csv(file = "output/data/chapchap_intersections.csv") %>%
                    select(from=Var1, to=Var2, value, cancer) %>% 
                    mutate(prots=value*1143)



excel_path <- "HPC/binari_validated_corrs.xlsx"
sheet_names <- excel_sheets(excel_path)
pdf("output/cancer_similarity_tail_chapchap_network.pdf")
for (cancer_nm in sheet_names) {
  print_cancer_igraph(cancer_nm, chap_attrib, intersections_df)
}
dev.off()

# ------ prepare networks for figure ---
all_simlrs_long <- read.csv("output/jaccard_values_per_cancer_long_format.csv")
in_tail <- all_simlrs_long %>% filter(similarity > 0.3)

chap_attrib <- read.csv("output/data/chap_attributes.csv") %>%
  arrange(v_order) %>% 
  mutate(shape='circle') %>%
  mutate(colour=case_when(module==1 ~ 'lightblue',
                          module==2 ~ 'pink',
                          module==3 ~ 'lightgreen'))

intersections_df <- read.csv(file = "output/data/chapchap_intersections.csv") %>%
  select(from=Var1, to=Var2, value, cancer) %>% 
  mutate(prots=value*1143)

pdf("output/cancer_similarity_samp_networks.pdf")
for (cancer_nm in c("PRAD","LUSC","HNSC","LIHC")) {
  print_cancer_igraph(cancer_nm, chap_attrib, intersections_df)
}
dev.off()

# all cancers togather ---- 
# distribution in tail  by cancer
p2<-ggplot(in_tail, aes(x=similarity, fill= cancer)) + 
  geom_histogram(bins=20, color="black", position="stack", alpha=0.5) + 
  labs(x="Jaccard similarity index", fill = "Cancer")
p2

n_distinct(c(in_tail$col, in_tail$row)) # all 15 chaps are in the tail.
table(c(in_tail$col, in_tail$row)) # how do they distribute in the chaps?

# igraph for all cancers in the same plot
rrr <- igraph::graph_from_data_frame(d = in_tail, vertices = chap_attrib, directed = FALSE)
plot.igraph(rrr,  axes = FALSE, vertex.frame.color = NA,
            #vertex.label = V(g)$name, vertex.label.color = "gray20",
            vertex.size = 25, #ideg*25 + 40, vertex.size2 = 30,
            vertex.color = chap_attrib$colour, #"gray90", vertex.frame.color = "gray20",
            vertex.shape = chap_attrib$shape,
            #edge.arrow.size=0.5, edge.color=col, edge.width = E(g)$weight / 10,
            #edge.curved = T,
            main = "all",
            margin = c(-0.4,-0.4,-0.15,-0.2),
            layout = as.matrix(chap_attrib[2:3]))


# ------ in modules vs between modules -----
# investigate for each cancer what is the percent of connection 
# inside each module, the percent between modules
modules_connections_by_filter <- function(cutoff=0.4) {

  all_simlrs_long <- read.csv("output/jaccard_values_per_cancer_long_format.csv")
  in_tail <- all_simlrs_long %>% filter(similarity > cutoff)
  head(in_tail)
  
  chap_attrib <- read.csv("output/data/chap_attributes.csv") %>%
                      select(chap=X, prot_type=function., module)
  head(chap_attrib)
  with_modules <- in_tail %>% 
      left_join(chap_attrib, by = c("row" = "chap")) %>%
      left_join(chap_attrib, by = c("col" = "chap")) %>%
      select(from=row, to=col, cancer, from_module=module.x, to_module=module.y)
  head(with_modules)
  
  # calculate percentages
  res <- with_modules %>%
         group_by(cancer) %>%
         summarise(in_module_num=length(cancer[from_module==to_module]), 
                   between_module_num=length(cancer[from_module!=to_module]), 
                   all=length(cancer)) %>% 
         mutate(in_percent=in_module_num/all, 
                out_percent=between_module_num/all)
  return(res)
}

modules_connections_by_filter(0.1)
modules_connections_by_filter(0.2)
modules_connections_by_filter(0.3)
info <- modules_connections_by_filter(0.4)

info %>% ggplot(aes(x=cancer, y= in_percent)) + 
  geom_boxplot() + 
  labs(x=element_blank(), y="Percent of connections within modules") + 
  paper_figs_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
