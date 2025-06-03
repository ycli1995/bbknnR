# bbknnR
Use batch balanced KNN (BBKNN) in R

## Introduction
BBKNN is a fast and intuitive batch effect removal tool for single-cell data. It is originally used in the [scanpy](https://scanpy.readthedocs.io/en/stable/) workflow, and now can be used with [Seurat](https://satijalab.org/seurat/) seamlessly.

### System requirements
`bbknnR` has been tested on R versions >= 4.1. Please consult the `DESCRIPTION` file for more details on required R packages. bbknnR has been tested on Linux platforms

To use the full features of `bbknnR`, you also need to install the [bbknn](https://bbknn.readthedocs.io/en/latest/) python package:
```
pip install bbknn
```

### Installation
`bbknnR` has been released to CRAN:
```
install.packages("bbknnR")
```
or can be installed from github:
```
devtools::install_github("ycli1995/bbknnR")
```

### Quick start
```
library(bbknnR)
library(Seurat)
data("panc8_small")
panc8_small <- RunBBKNN(panc8_small, batch_key = "tech")
```

## Release
### 2.0.0
* Remove `reticulate` dependency. Now use kNN algorithms provided by `RcppAnnoy` and `rnndescent`
* Add `return.umap.model` for `RunBBKNN.Seurat`
* Improvements for `testthat`

### 1.1.0
* Compatibility with Seurat v5
* Improvements for documentation and verbose.

### 1.0.2
* Explicit import of `get_dummies.()` from tidytable
* Fix a bug when pass only one `batch_key` to `RidgeRegression()`

### 1.0.1
* Import public function `similarity_graph()` from `uwot==0.1.14` in `compute_connectivities_umap()` to follow the CRAN policy

### 1.0.0
* Initially released to CRAN

## Citation
Please cite this implementation R in if you use it:
```
Yuchen Li (2022). bbknnR: Use batch balanced KNN (BBKNN) in R.
package version 0.1.0 https://github.com/ycli1995/bbknnR
```

Please also cite the original publication of this algorithm.
```
Polanski, Krzysztof, et al. "BBKNN: fast batch alignment of single cell transcriptomes." Bioinformatics 36.3 (2020): 964-965.
```
