# IDclust
IDclust is an unsupervised method for clustering single-cell datasets (scRNA, scEpigenome). The cells are iteratively re-processed and re-clustered. Subclusters are created provided they are significantly different from each other.

![alt text](https://github.com/vallotlab/IDclust/blob/master/inst/www/scheme.png?raw=true)

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

## Running Iterative Differential Clustering 

```
library(IDclust)
library(Seurat)

data("Seu")
Seu <- iterative_differential_clustering(
    Seu,
    output_dir = "./",
    plotting = TRUE,
    saving = TRUE,
    n_dims = 50,
    dim_red = "pca",
    vizualization_dim_red = "umap",
    processing_function = processing_Seurat,
    differential_function = differential_edgeR_pseudobulk_LRT,
    logFC.th = log2(1.5),
    qval.th = 0.01
    )
```

## Plotting network
## scEpigenomics 

```
plot_cluster_network(Seu)
```

![alt text](https://github.com/vallotlab/IDclust/blob/master/inst/www/network.png?raw=true)


