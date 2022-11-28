# Meta-analysis of proteasomal gene expression in T cells from human cancers

This repository contains scripts to reproduce bioinformatics analyses from the work by de De Blas et al. (manuscript in preparation). Briefly, we re-analyse the large collection of scRNA-seq dataset collated by [Zheng et al. (2021) Science](https://www.science.org/doi/10.1126/science.abe6474) and analyse the association of proteasomal gene expression to CD8+ T cell exhaustion.

The annotated notebook [classify_cancer_samples.Rmd](./classify_cancer_samples.Rmd) downloads and processes a collection of CD8+ T cells from multiple studies and cancer types from the study by Zheng et al. For this meta-analysis, we focused on droplet-based technologies (25 datasets covering 19 cancer types), to ensure that expression profiles are comparable across studies.

High-quality T cells are filtered using [scGate](https://github.com/carmonalab/scGate), removing potential contaminants and excluding very small samples (at least 200 cells). The [ProjecTILs](https://github.com/carmonalab/ProjecTILs) algorithm and a mouse TIL reference map are then used to annotate the datasets into CD8+ T cell functional subtypes. 

Focusing on T exhausted cells (Tex) and T effector memory cells (Tem), we carried out differential gene expression analysis. Besides the expected markers for these subtypes, this also revealed several genes encoding for proteasomal subunits. Gene signatures for proteasomal genes, exhaustion and ROS detoxification were evaluated using [UCell](https://bioconductor.org/packages/release/bioc/html/UCell.html). Additional signatures were obtained from the [MSigDB database](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) from the Broad Institute.
 
