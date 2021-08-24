#------------------------------
# chaperones' couple therapy
# ron over all the chaperones that have a common protein, generate a table that 
# shows how consistant each chaperone-chaperone-client triangle is.
#------------------------------

#-------- includes --------
library(readxl)
library(infomapecology)



#-------- prepare protein meta data -------- TODO: make sure this is needed
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


#-------- run multilayer infomap across relaxed values -------
# https://ecological-complexity-lab.github.io/infomap_ecology_package/
# TODO not done at all!
# TODO this is copied from Shai's example. need to be adjusted.

# Create a multilayer object
NEE2017 <- create_multilayer_object(extended = siberia1982_7_links, 
                                    nodes = siberia1982_7_nodes, 
                                    intra_output_extended = F, 
                                    inter_output_extended = F)

# Ignore interlayer edges
NEE2017$inter <- NULL

#Run Infomap and return the network's modular structure at increasing relax-rates.
relaxrate_modules <- NULL
for (r in seq(0,1,0.05)){
  print(r)
  modules_relax_rate <- run_infomap_multilayer(NEE2017, relax = T, silent = T, 
                                               trials = 50, seed = 497294, 
                                               multilayer_relax_rate = r, 
                                               multilayer_relax_limit_up = 1, 
                                               multilayer_relax_limit_down = 0, 
                                               temporal_network = T)
  modules_relax_rate$modules$relax_rate <- r
  relaxrate_modules <- rbind(relaxrate_modules, modules_relax_rate$modules)
}




# Other:
# here run the multilayer on the r=0


I_or <- NULL # I between observed and relaxed
for (r in seq(0.05,1,0.05)){
  print(r)
  mods_relax <- run_infomap_multilayer(M = mln_relax, flow_model = 'undirected', silent = T, 
                                       trials = 20, seed = 12345,
                                       relax = T, 
                                       multilayer_relax_rate = r,
                                       multilayer_relax_limit = -1, 
                                       temporal_network = F)
  A <- mods$modules %>% select(module, node_id, layer_id)
  B <- mods_relax$modules %>% select(module, node_id, layer_id)
  
  N <- inner_join(A,B,by=c('node_id','layer_id')) %>% 
    group_by(module.y) %>% 
    select(module.x) %>% table()
  I_or <- rbind(I_or, tibble(r=r, I=NMI(N))) 
}


# TODO: download the package in linux, 
# make sure you can run infomap through the code
# prepare data of in the format to be run
# sanity check to the code shai sent me
# run the multilayer infomap on my cancer networks

