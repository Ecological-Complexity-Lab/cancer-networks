#------------------------------
# chaperones' couple therapy
# ron over all the chaperones that have a common protein, generate a table that 
# shows how consistant each chaperone-chaperone-client triangle is.
#
# This scripts uses Infomap and assumes it is installed in its working directory.
# For more information:
# https://ecological-complexity-lab.github.io/infomap_ecology_package
#------------------------------

#-------- includes --------
library(readxl)
library(infomapecology)

#-------- prepare protein meta data -------- 
prots_meta <- read.table("HPC/Mito_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
chaps_meta <- read.table("HPC/Mito_ch_genes.tab", sep="\t", header=TRUE, 
                         stringsAsFactors=FALSE, quote="", fill=FALSE)
prots_meta <- prots_meta[!prots_meta$ENSID %in% chaps_meta$ENSID, ]


#-------- load the networks from an excel file --------
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

#-------- build infomap format to be run --------
# build node tibble
chap_nodes <- cbind.data.frame(1:nrow(chaps_meta),"chaperone",
                    chaps_meta["Symbol"], chaps_meta["ENSID"])
colnames(chap_nodes) <- c("node_id", "type", "symbol", "ENSID")
prot_start_id <- nrow(chap_nodes) + 1
prot_nodes <- cbind.data.frame(prot_start_id:(nrow(prots_meta)+prot_start_id-1),
                           "protein", prots_meta["Symbol"], prots_meta["ENSID"])
colnames(prot_nodes) <- c("node_id", "type", "symbol", "ENSID")
all_nodes <- rbind.data.frame(chap_nodes, prot_nodes)

# make symbols as row number
rownames(chap_nodes) <- chap_nodes[,"symbol"]
rownames(prot_nodes) <- prot_nodes[,"ENSID"]

all_link <- data.frame(layer_from=numeric(),
                       node_from=numeric(),
                       layer_to=numeric(),
                       node_to=numeric(),
                       weight=numeric())

# build bridges dateframe
cncr_nms <- names(networks)
for (cancer_id in 1:length(networks)) {
  net <- networks[[cancer_id]]
  
  for (chap in chap_nodes$symbol) {
    chap_id <- as.numeric(chap_nodes[chap_nodes$symbol==chap, "node_id"])
    
    for (prot in prot_nodes$ENSID) {
      prot_id <- as.numeric(prot_nodes[prot_nodes$ENSID==prot, "node_id"])
      link_value <- net[chap, prot]
      
      # if there's an interaction to report
      if (link_value==1) {
        all_link <- rbind(all_link, c(cancer_id, chap_id, cancer_id, prot_id, 1))
      }
    }
  }
}
names(all_link) <- c("layer_from","node_from","layer_to","node_to","weight")

# remove nodes that have no links in any cancer
active_nodes <- all_nodes %>% 
  filter(node_id %in% all_link$node_to | (node_id %in% all_link$node_from))

# create a table for cancer ids
cncr_nms <- names(networks)
all_layers <- tibble(layer_id=1:length(cncr_nms), cancer=cncr_nms)

# turn them into tibbles
tbl_nodes <- as_tibble(active_nodes[,-4])
tbl_links <- as_tibble(all_link)

#-------- run multilayer infomap across relaxed values -------
# prepare multilayer object
net_obj <- create_multilayer_object(extended = tbl_links, 
                                    nodes = tbl_nodes, 
                                    intra_output_extended = F, 
                                    inter_output_extended = F)

# here run the multi-layer on the r=0
r_is_0 <- run_infomap_multilayer(M = net_obj, 
                                 flow_model = 'undirected', 
                                 silent = TRUE, 
                                 trials = 20, # change to 100 if this is not slow
                                 seed = 1234,
                                 relax = TRUE, 
                                 multilayer_relax_rate = 0,
                                 multilayer_relax_limit = -1, 
                                 temporal_network = FALSE)
# see the number of nodes in each module
table(r_is_0$modules$module)

A <- r_is_0$modules %>% select(module, node_id, layer_id)
all_modules <- A
all_modules["relax_param"] <- 0 

I_or <- rbind(NULL, tibble(r=0, I=1)) # save I between observed and relaxed
for (r in seq(0.05,1,0.05)){
  print(r)
  mods_relax <- run_infomap_multilayer(M = net_obj, 
                                       flow_model = 'undirected', 
                                       silent = TRUE, 
                                       trials = 20, 
                                       seed = 1234,
                                       relax = TRUE, 
                                       multilayer_relax_rate = r,
                                       multilayer_relax_limit = -1, 
                                       temporal_network = FALSE)
  B <- mods_relax$modules %>% select(module, node_id, layer_id)
  all_modules <- rbind(all_modules, data.frame(B, relax_param=r))
  
  N <- inner_join(A,B,by=c('node_id','layer_id')) %>% 
       group_by(module.y) %>% 
       select(module.y, module.x) %>% table()
  I_or <- rbind(I_or, tibble(r=r, I=NMI(N))) 
}

concluting_table <- left_join(all_modules, all_nodes, by="node_id") %>% 
                    left_join(all_layers, by="layer_id") %>% 
                    select("relax_param", "cancer", "symbol", "module", "type","ENSID")

# get the number of modules fer relax plot
modules_per_cancer <- concluting_table %>% select(relax_param, module) %>%
                      group_by(relax_param) %>% summarise(module = max(module))
plot(modules_per_cancer)

# Save modules to a file
write.csv(concluting_table, file = "output/multilayer_relaxed_scan_20_trials.csv")

View(I_or)
# NMI of AxA results in 1.
# NMI of Axf(r) when 0.20<r<0.55 is 0 for some reason: 
# all the nodes are in the same cluster(??)

concluting_table <- read_csv('output/multilayer_relaxed_scan_20_trials.csv')

ggplot(modules_per_cancer, aes(relax_param, module))+geom_point()+geom_line()+
  theme(axis.text = element_text(size=20))


concluting_table %>% 
  filter(relax_param==0.15) %>% 
  filter(type=='chaperone') %>%
  group_by(cancer, module) %>% 
  summarise(n=n_distinct(symbol)) %>% 
  ggplot(aes(cancer, module, fill=n, label=n))+geom_tile()+
  geom_text()+
  theme(axis.text = element_text(size=20))

concluting_table %>% 
  filter(relax_param==0.15) %>% 
  filter(type=='chaperone') %>%
  group_by(symbol, module) %>% 
  summarise(n=n_distinct(cancer)) %>% 
  ggplot(aes(symbol, module, fill=n, label=n))+geom_tile(color='red')+
  geom_text()+
  theme(axis.text = element_text(size=20))

concluting_table %>% 
  filter(relax_param==0.15) %>% 
  filter(type=='chaperone') %>%
  filter(module==1) %>% 
  group_by(symbol, cancer) %>% 
  summarise(n=n_distinct(type)) %>% 
  ggplot(aes(symbol, cancer, fill=n, label=n))+geom_tile(color='red')+
  geom_text()+
  theme(axis.text = element_text(size=20))
