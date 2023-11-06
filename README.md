# expHRD
An Individualized, Transcriptome-based Prediction Model for Homologous Recombination Deficiency (HRD) Assessment in Cancer


## Description
The expHRD repository contains the source code for predicting expHRD scores in cancer patients. These scores are calculated using an algorithm that relies on transcriptome-based n-of-1 style HRD scoring.


## Packages required
To run this code, you will need the following packages:
  - R (v 4.2.0)
  - GSVA (v 1.46.0)
  - ggplot2 (v 3.4.0)
  - ggpubr (v 0.5.0)

## Installation
To get started, clone this repository using the following command:
  ```bash
  git clone https://github.com/kimdh4962/expHRD
  ```

## Examples
Here's how you can use the expHRD code with an example:

- Input: RNA-seq data of triple negative breast cancer (TNBC) patients from The Cancer Genome Atlas Program (TCGA).
- Output: Predicted HRD scores for the patients.
To run the code, follow these steps:

1. Navigate to the expHRD directory:
```bash
cd expHRD/code
```
2. Execute the R script for HRD score prediction:
```bash
Rscript CAL_exp_HRD_github.R
```
3. Check the output file for the prediction:
```bash
vi ./../output/OUT_TABLE.csv
```
<br><hr/>
Feel free to replace the input data with your own if needed. If you encounter any issues or have questions, please refer to the documentation or open an issue on this repository.
