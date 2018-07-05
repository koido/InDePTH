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

First of all, please collect public data descrived in the help file of init_fun.R. We cannot redistribute these public files.

Then, run:
```r
# directory is only required.
init(input.dir = "input", out.dir = "init")
```

This make all files required for InDePTH.

### Example of InDePTH with arbitrary threshold

```r
# cscore calculation
res1 <- cscore_LINCS(up.sig = c( "217761_at", "217398_x_at", "218744_s_at" ),
                     dn.sig = c( "221856_s_at", "218193_s_at", "218457_s_at" ),
                     input.dir = "init",
                     output.dir = "out",
                     write.name = "test1.tsv")
# influence calculation
res2 <- influence_calc(up.sig = c( "217761_at", "217398_x_at", "218744_s_at" ),
                      dn.sig = c( "221856_s_at", "218193_s_at", "218457_s_at" ),
                      up.ratio = c( 2, 4, 6 ),
                      dn.ratio = c( 1/2, 1/3, 1/5 ),
                      tot.thr = 0.5,
                      input.dir = "init",
                      output.dir = "out",
                      cscore = res1,
                      write.name = "test1")
```

### Example of InDePTH with data-driven threshold)

```r
# threshold calculation using res1
thr_res <- make_threshold( pert_id_vec = c( "TRCN0000010389", "TRCN0000010390", "TRCN0000010391"),
                           cscore = res1,
                           cell = "HT29",
                           input.dir = "init",
                           output.dir = "out" )
# influence calculation
res2 <- influence_calc(up.sig = c( "217761_at", "217398_x_at", "218744_s_at" ),
                      dn.sig = c( "221856_s_at", "218193_s_at", "218457_s_at" ),
                      up.ratio = c( 2, 4, 6 ),
                      dn.ratio = c( 1/2, 1/3, 1/5 ),
                      tot.thr = thr_res$cutoff,
                      input.dir = "init",
                      output.dir = "out",
                      cscore = res1,
                      write.name = "test1_thr")
```

Now, res2 contains information of hub score. Using test1.rds or test1_thr.rds, the gene regulatory network can be visualized by the R {igraph} package.

The all analysis conditions which were described here are dummy.
Please see each help file for the detail usage.

## Reference
* M. Koido, Y. Tani, S. Tsukahara, Y. Okamoto, A. Tomida, “InDePTH: Detection of hub genes for developing gene expression networks under anticancer drug treatment”, *Oncotarget*, 9(49), 29097-29111, 2018

## License

GPL-3.
