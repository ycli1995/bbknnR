
#' Perform batch balanced KNN
#' 
#' Batch balanced KNN, altering the KNN procedure to identify each cell’s top 
#' neighbours in each batch separately instead of the entire cell pool with no 
#' accounting for batch. The nearest neighbours for each batch are then merged 
#' to create a final list of neighbours for the cell. Aligns batches in a quick 
#' and lightweight manner.
#' 
#' @param object An object
#' @param ... Arguments passed to other methods
#' 
#' @references 
#' Polański, Krzysztof, et al. "BBKNN: fast batch alignment of single cell 
#' transcriptomes." Bioinformatics 36.3 (2020): 964-965.
#' 
#' @export RunBBKNN
RunBBKNN <- function(object, ...) {
  UseMethod(generic = 'RunBBKNN', object = object)
}

#' Perform ridge regression on scaled expression data
#' 
#' Perform ridge regression on scaled expression data, accepting both technical 
#' and biological categorical variables. The effect of the technical variables 
#' is removed while the effect of the biological variables is retained. This is 
#' a preprocessing step that can aid BBKNN integration.
#' 
#' @param object An object
#' @param ... Arguments passed to other methods
#' 
#' @references Park, Jong-Eun, et al. "A cell atlas of human thymic development 
#' defines T cell repertoire formation." Science 367.6480 (2020): eaay3224.
#' 
#' @export RidgeRegression
RidgeRegression <- function(object, ...) {
  UseMethod(generic = 'RidgeRegression', object = object)
}
