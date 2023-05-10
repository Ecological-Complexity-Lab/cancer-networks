
# Ecological network analysis reveals cancer-dependent chaperone-client interaction structure and robustnes
Published in Nature Communications,  
by
Geut Galai,
Xie He,
[Barak Rotblat](https://barakrotblat.wixsite.com/rotblatlab),
and [Shai Pilosof](https://lifewp.bgu.ac.il/wp/pilos/)

This repository contains all the data and source code used to produce the results
presented in this paper:

> TODO: add references of the paper and/or preprint when available.

|  | Info |
|-:|:-----|
| Version of record | https://doi.org/JOURNAL_DOI |
| Open-access version on biorxiv | https://doi.org/BIORXIV_DOI |
| Archive of this repository | https://doi.org/FIGSHARE_ZENODO_DOI |

## About

This paper is a result of a collaboration between the [Ecological Complexity Lab](https://lifewp.bgu.ac.il/wp/pilos/) and the [Rotblat Lab](https://barakrotblat.wixsite.com/rotblatlab). It sets out to view cancer types as different ecological environments and as such regard the actors in them as ecological interactors thus revealing how different environments effect those interactions. Here we chose to study mitochondrial chaperone-clients interaction networks (similar to plant-pollinator networks in ecological studies) using co-expression correlations.


## Abstract

Cancer cells alter the expression levels of metabolic enzymes to fuel proliferation. The mitochondrion is a central hub of metabolic reprogramming, where chaperones service hundreds of clients, forming chaperone-client interaction networks (CCINs). How network structure affects its robustness to chaperone targeting is key to developing cancer-specific drug therapy. However, how structure and robustness vary across different cancer tissues remains unknown. We revealed a non-random, hierarchical pattern whereby the cancer type modulates the chaperones' ability to realize their potential client interactions. Despite the low similarity between the CCINs, we highly accurately predicted links in one cancer type based on another. Moreover, we identify groups of chaperones that interact with similar clients. Simulations of network robustness show that this group structure affects cancer-specific response to chaperone removal. Our results open the door for new hypotheses regarding the ecology and evolution of CCINs and can inform caner-specific drug development strategies.


# System requirements
R Programming language: 4.1.0

Python: 2.7.18

Link Prediction: tested on  Linux 6.2, Ubuntu 23.04

Infomap (MacOS version): Version 1.7.1

tested on: MacOS BigSur 11.6.6

HPC OS version: Oracle Linux Server 8.7


## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/Ecological-Complexity-Lab/cancer-networks.git


## Reproducing the results

The software in this study can be divided into 3 parts -

1. R scripts compatible with running on an HPC server as independent jobs - this includes:
* Building the initial binary tables from the expression data (`build_binari_table.R`) for each cancer type.
* Running the analysis to validate the co-expression method using external data from STRING database (`chap_stringdb_affirmation.r`). (fig. S6A in the paper)
Everything needed is available in the `HPC` folder (including the results used in the paper).

2. Python scripts to be run on linux server with this includes:
* Running SBM analysis itself. Access code_for_link_prediction&community_detection, run `Best_Community_Arrangement.py`. Note that you could choose either to run on the projected network or the original bipartite network based on your need. 
* Running link prediction analysis. Access code_for_link_prediction&community_detection, run `run_link_prediction.py` with Python: 2.7.18. To change the number of predefined community, change the number of "K" in that code, detailed comments could be found in the code.

3. R scripts meant to be ran on a local computer (on our case a Mac), including:
* Producing the format needed for the SBM analysis and processing the results (`SBM_analysis.r`).
* Preforming the rest of the scripts presented in the paper.
All the results will be found under the `output` folder.


## License

> TODO: what license does our code operate under? how should I adjust this section to be relevant to us?

All Python source code (including `.py` and `.ipynb` files) is made available
under the MIT license. You can freely use and modify the code, without
warranty, so long as you provide attribution to the authors. See
`LICENSE-MIT.txt` for the full license text.

The manuscript text (including all LaTeX files), figures, and data/models
produced as part of this research are available under the [Creative Commons
Attribution 4.0 License (CC-BY)][cc-by]. See `LICENSE-CC-BY.txt` for the full
license text.

[cc-by]: https://creativecommons.org/licenses/by/4.0/
