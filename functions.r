# ------------ "functions.r" ------------------
# a file containing functions that repeat and are common 
# to multiple script files
# ------------------------------

# ----- global consts ------
cancer_sample_size_order <- c("KIRP", # 288
                              "LIHC", # 371
                              "STAD", # 375
                              "COAD", # 478
                              "PRAD", # 498
                              "HNSC", # 500
                              "LUSC", # 502
                              "THCA", # 502
                              "LUAD", # 533
                              "KIRC", # 538
                              "UCEC", # 551
                              "BRCA") # 1102
cancer_percent_of_valid_interactions <- c("BRCA", # 22%
                                          "LUSC", # 31%
                                          "LUAD", # 33%
                                          "HNSC", # 48%
                                          "KIRC", # 49%
                                          "UCEC", # 49%
                                          "COAD", # 50%
                                          "LIHC", # 58%
                                          "PRAD", # 58%
                                          "STAD", # 61%
                                          "THCA", # 61%
                                          "KIRP") # 100%
chap_module_order <- 
  c("HSPD1", "HSPA9", "CLPX", "TRAP1", "AFG3L2", "DNAJA3", "GRPEL2", "YME1L1",
    "HSPE1", "LONP1", "CLPP", "HSCB", "DNAJC19", "HTRA2", 
    "SPG7")
# ----- global functions ------

# load the cancer network matrices
load_cancer_mats <- function(excel_path="HPC/binari_validated_corrs.xlsx") {
  # do the convert for every cancer
  sheet_names <- excel_sheets(excel_path)
  
  networks <- list()
  for (name in sheet_names) {
    x <- as.data.frame(read_excel(excel_path, sheet = name, col_names = TRUE))
    x2 <- x[,-1]
    rownames(x2) <- x[,1]
    networks[[name]] <- x2
  }
  
  return(networks)
}

# build a cancer union matrix
get_cancer_union <- function(networks) { 
  # assumption: all tables have the same structure
  union_table <- networks[[1]]
  union_table[,] <- 0 
  
  for (i in 1:length(networks)) {
    net <- networks[[i]]
    union_table <- union_table + net
  }

  return(union_table)
}

# ---------- ploting themes -----

paper_figs_theme <- 
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank(), 
        legend.text=element_text(size=12, color='black'),
        legend.title=element_text(size=14, color='black'))

paper_figs_theme_no_legend <- 
  paper_figs_theme +
  theme(legend.position = 'none')