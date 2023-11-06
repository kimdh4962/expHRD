# expHRD
An Individualized, Transcriptome-based Prediction Model for Homologous Recombination Deficiency (HRD) Assessment in Cancer


## Description
This is a source code for expHRD, An Individualized, Transcriptome-based Prediction Model for Homologous Recombination Deficiency Assessment in Cancer.


## Packages required
  - R (v 4.2.0)
  - GSVA (v 1.46.0)
  - ggplot2 (v 3.4.0)
  - ggpubr (v 0.5.0)

## Installation

    git clone https://github.com/kimdh4962/expHRD



## Examples
  - Input: RNA-seq data of triple negative breast cancer (TNBC) patients from The Cancer Genome Atlas Program (TCGA)
  - Output: Predicted HRD scores of the patients

    cd expHRD/code
    Rscript CAL_exp_HRD_github.R
