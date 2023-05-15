
# Ecological network analysis reveals cancer-dependent chaperone-client interaction structure and robustnes
Published in Nature Communications 2023, by Geut Galai, [Xie He](https://sites.google.com/view/xie-he/about-me), [Barak Rotblat](https://barakrotblat.wixsite.com/rotblatlab),
and [Shai Pilosof](https://lifewp.bgu.ac.il/wp/pilos/).

**Please cite the paper when using the data or code.**

This repository contains all the data and source code used to produce the results
presented in this paper:

> TODO: add references of the paper and/or preprint when available.

|Item  | Location |
|:-|:-----|
| Preprint | [https://doi.org/10.1101/2022.11.24.517852](https://doi.org/10.1101/2022.11.24.517852) |
| Data | Deposited on Figshare. DOI: 10.6084/m9.figshare.22779755 |

## About

This paper is a result of a collaboration between the [Ecological Complexity Lab](https://lifewp.bgu.ac.il/wp/pilos/) and the [Rotblat Lab](https://barakrotblat.wixsite.com/rotblatlab). We set out to understand how chaperone-client interaction networks are structured, how structure varies across cancer environments and what are the consequences of this variation to the networks' robustness to chaperon targeting. Such knowledge is crucial for understanding the molecular mechanisms underlying metabolic reprogramming, with potential applications for cancer-specific drug development. To address these questions we adopt theory and methodology from the discipline of Network Ecology. Analogous studies in ecology aim to discover non-random variation in species interactions (e.g., plant-pollinator) across environments.


## Abstract

Cancer cells alter the expression levels of metabolic enzymes to fuel proliferation. The mitochondrion is a central hub of metabolic reprogramming, where chaperones service hundreds of clients, forming chaperone-client interaction networks (CCINs). How network structure affects its robustness to chaperone targeting is key to developing cancer-specific drug therapy. However, how structure and robustness vary across different cancer tissues remains unknown. We revealed a non-random, hierarchical pattern whereby the cancer type modulates the chaperones' ability to realize their potential client interactions. Despite the low similarity between the CCINs, we highly accurately predicted links in one cancer type based on another. Moreover, we identify groups of chaperones that interact with similar clients. Simulations of network robustness show that this group structure affects cancer-specific response to chaperone removal. Our results open the door for new hypotheses regarding the ecology and evolution of CCINs and can inform caner-specific drug development strategies.


# System requirements
R Programming language: 4.1.0 \
- 'vegan' package: 2.5-7 \
- 'bipartite' package: 2.16

Python: 2.7.18 

Infomap (MacOS version): Version 1.7.1

tested on: \
- R code - MacOS BigSur 11.6.6 \
- Python Link Prediction Analysis - Linux 6.2, Ubuntu 23.04 (Detailed requirement could be found in the req.txt)\
- HPC OS version: Oracle Linux Server 8.7 

## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/Ecological-Complexity-Lab/cancer-networks.git

## Folder organization

The main repository folder contains the code performing the analysis presented 
in the paper that are done locally (not on the HPC). \

Key script in the main folder: \
`calc_folding_efficiency_per_cancer.R`: Implementation of the realized niche analysis. \
`cancer_similarity.R`: Calculate Jaccard similarity between each pair of chaperones for a given cancer.  \
`chap_expression_analysis.R`: Process and analyse chaperone expressions\
`chaps_similerity_index.R`: Calculate Jaccard similarity between each pair of cancers for a given chaperone.  \
`functions.R`: Functions etc used in more then one script to prevent code duplication.\
`interaction_evidance.R`: Implementation of the co-expression affirmation using external resources. \
`multilayer_infomap.R`: Community detection analysis using Infomap. \
`shuffle_cancer_networks.R`: Implements the shuffling of cancer networks used for analysis validations.\
`stability_check.R`: The robustness analysis implementation (including correlations)\
`test_expressionXfolding_correlation.R`: Checks the correlation between chaperone expression and realized niche.


Nested folders: \
`code_for_link_prediction&community_detection`: As the name suggests, contains the code for the link prediction analysis and the community detection, using SBM. \
`data_for_lp&cd`: This folder holds the data used by the code above. \
`external_data`: Data taken from other papers and STRING database for the affirmation analysis (fig. S6). \
`HPC`: This folder contains the code, the input and the output of the analysis run on the HPC.  \
`output`: All output files from the local code is directed to this folder. \
`prediction_data_80_k2`: All the original output files from the link prediction results. The cleaned up version of this could be found as a csv file in `code_for_link_prediction&community_detection` named as `link_prediction_80_k2.csv` \
`Source Data For Figures`: contains the zip code that holds the information detailing the what is shown in the figures of the papers. \


## Reproducing the results

The software in this study can be divided into 3 parts -

1. R scripts compatible with running on an HPC server as independent jobs - this includes:
* Building the initial binary tables from the expression data (`build_binari_table.R`) for each cancer type.
* Running the analysis to validate the co-expression method using external data from STRING database (`chap_stringdb_affirmation.r`). (fig. S6A in the paper)
Everything needed is available in the `HPC` folder (including the results used in the paper).

2. Python scripts to be run on linux server with this includes:
* Running SBM analysis itself. Access code_for_link_prediction&community_detection, run `Best_Community_Arrangement.py`. Note that you could choose either to run on the projected network or the original bipartite network based on your need. 
* Running link prediction analysis. Access code_for_link_prediction&community_detection, run `run_link_prediction.py` with Python: 2.7.18. To change the number of predefined community, change the number of "K" in that code, detailed comments could be found in the code.

3. R scripts meant to be ran on a local computer, including:
* Producing the format needed for the SBM analysis and processing the results (`SBM_analysis.r`).
* Preforming the rest of the scripts presented in the paper.
All the results will be found under the `output` folder.


## License

All source code is made available under the CC BY-NC-SA license. This license lets you remix, tweak, and build upon this work non-commercially, as long as you credit the authors and license the new creations under the identical terms.
