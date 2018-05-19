# InDePTH: influential gene detection in perturbed transcriptome hierarchical network

## Introduction

InDePTH is a novel algorithm to detect hubs of influential genes from reconstructed upstream and downstream relationships among differentially expressed genes (DEGs) (user-prepared, query DEGs), by referring to the rank matrix of Z-scores from a database of comprehensive genetic perturbations, such as the LINCS L1000 dataset (publicly available, reference data).

## Installation

```r
devtools::install_github("koido/InDePTH")
```

## Usage

Currently, InDePTH is composed of 4 Rscript files:
* init_fun.R
* cscore_calc.R
* make_threshold.R
* influence_calc.R

Please see each help file.

## Reference
* M. Koido, Y. Tani, S. Tsukahara, Y. Okamoto, A. Tomida, “InDePTH: Detection of hub genes for developing gene expression networks under anticancer drug treatment”, *Oncotarget*, accepted
