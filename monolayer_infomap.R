#------------------------------
# This scripts uses Infomap and assumes it is installed in its working directory.
# For more information:
# https://ecological-complexity-lab.github.io/infomap_ecology_package
#
# the script runs only on one cancer at a time
#------------------------------

#-------- includes --------
library(readxl)
library(infomapecology)
library(tidyr)

source("functions.r")

#-------- load the networks from an excel file --------
networks <- load_cancer_mats()

#-------- run infomap on each cancer ---- 
all_mono_modules <- NULL
plots <- list()
for (cncr in sheet_names) {
  net <- networks[[cncr]]
  
  # run multilayer infomap across relaxed values:
  net_obj <- create_monolayer_object(x = t(as.matrix(net)), 
                                    directed = FALSE, 
                                    bipartite = TRUE, 
                                    group_names = c("Chaperones","Proteins"))
  
  # here run the infomap
  res <- run_infomap_monolayer(x = net_obj, 
                                   flow_model = 'undirected',
                                   silent = TRUE, 
                                   trials = 100, # change to 100 if this is not slow
                                   seed = 1234)
  # see the number of nodes in each module
  table(res$modules$module_level1)
  
  # save results
  mono_modules <- res$modules %>% 
                  select(id = node_id, name = node_name, 
                         module = module_level1, type = node_group) %>% 
                  add_column(cancer = cncr)
  
  all_mono_modules <- rbind(all_mono_modules, 
                            mono_modules %>% filter(!is.na(module)))
}

write.csv(all_mono_modules, file = "output/monolayer_infomap_per_cancer.csv")

g <- all_mono_modules %>% filter(type=='Chaperones') %>%
  ggplot( aes(x=name, y=module, group=cancer, fill=cancer)) +
  geom_tile() + ylab("Module ID") +
  ggtitle("Chaperone distribution across modules per cancer") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.5, "lines"),
    strip.text.x = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_blank()
  ) +
  facet_wrap(~cancer)

ggsave("output/monolayer_infomap_per_cancer.pdf")
dev.off()
