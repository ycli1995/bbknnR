#' A small example version of the pancreas scRNA-seq dataset
#'
#' A subsetted version of the pancreas scRNA-seq dataset to test BBKNN
#'
#' @format A Seurat object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains one assay ("RNA" - scRNA-seq expression data)
#'   \item{counts - Raw expression data}
#'   \item{data - Normalized expression data}
#'   \item{scale.data - Scaled expression data}
#'   \item{var.features - names of the current features selected as variable}
#'   \item{meta.features - Assay level metadata such as mean and variance}
#'    }}
#'   \item{meta.data}{Cell level metadata}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{graphs}{Empty}
#'   \item{reductions}{Dimensional reductions: currently PCA}
#'   \item{version}{Seurat version used to create the object}
#'   \item{commands}{Command history}
#' }
#' @source SeuratData \url{https://github.com/satijalab/seurat-data}
#'
"panc8_small"