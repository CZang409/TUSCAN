# TUSCAN

`TUSCAN` is an R package that aims to perform tumor segmentation and classification analysis by constructing a spatial copy number profile. TUSCAN takes both the gene information from spatial transcriptomics data and the histology image as input to improve the accuracy of tumor region identification. Below is the overall workflow of TUSCAN, which consists of three main steps. The first step is to find a subset of normal spots. TUSCAN performs a preliminary
clustering of all spots based on the gene information dimensionally reduced by an autoencoder, along with RGB values extracted from the histology image. Then, it will automatically pick a cluster with the highest confidence of representing normal cells. The second step is to infer the copy number profile of all spots applying the selected cluster as the reference baseline. The third step is to segment all spots into tumors and non-tumors via a consensus clustering with the copy number profile as the input.

![Untitled (3)](https://github.com/CZang409/TUSCAN/assets/166551317/cdccc0c6-6feb-47ce-9782-f36a044eae2e)

## Software dependencies
- R version >= 4.0.0.  
- R packages:  magick, Seurat, cluster, ConsensusClusterPlus, dplyr, ggplot2, igraph, keras, tensorflow, dlm, RColorBrewer, ComplexHeatmap, grid, stats
  
## Installation
The R package can be installed from github:
```R
# Install devtools, if necessary
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github('CZang409/TUSCAN')
```

Note that TUSCAN relies on the keras and tensorflow R package, which can be installed through:
```R
install.packages("keras")
library(keras)
install_keras()
# install_keras() will install a version of TensorFlow that is compatible with the current version of the keras package.
```

TUSCAN depends on R package `ConsensusClusterPlus`, if you encouter this error when trying to install this package: "package 'ConsensusClusterPlus' is not available for this version of R", please install it through [Bioconductor](https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html):
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")
```

## Tutorial
Details in Tutorial
