# IDclust
IDclust is an unsupervised method for clustering single-cell datasets (scRNA, scEpigenome). The cells are iteratively re-processed and re-clustered. Subclusters are created provided they are significantly different from each other.

# Installation

Install package

```
if (!require("devtools"))
    {
      install.packages("devtools", dep=TRUE)
        if(!require("devtools")) stop("Package not found")
    }

devtools::install_github("vallotlab/IDclust")
```

## For scRNA use 
```
install.packages("Seurat")
```

## For scEpigenomics use

```
devtools::install_github("vallotlab/ChromSCape")
```

# Usage

## scRNA 

```
library(Seurat)
iterative_differential_clustering_scRNA(
    Seu,
    method = "Seurat",
    output_dir = "./",
    plotting = TRUE,
    saving = TRUE,
    nPCA = 50,
    nfeatures = 2000,
    logFC.th = log2(2),
    qval.th = 0.01,
    min.pct = 0.1,
    min.pct.cell_assigned = 0.25,
    limit = 10,
    k = 100,
    resolution = 0.5,
    biological_replicate_col = NULL
)
   
```

## scEpigenomics 

```
library(ChromSCape)
iterative_differential_clustering_scEpigenomics <- function(
    scExp,
    output_dir = "./",
    plotting = TRUE,
    saving = TRUE,
    nPCA = 10,
    percent_feature = 1,
    quantile.activation = 0.7,
    min.pct.cell_assigned = 0.1,
    FC.th = 2,
    qval.th = 0.1,
    limit = 5,
    k = 100,
    starting.resolution = 0.1,
    resolution = 0.8,
    verbose = TRUE
)
   
```

![alt text](https://github.com/vallotlab/IDclust/blob/master/inst/www/network.png?raw=true)


