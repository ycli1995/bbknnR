
annoy_knn <- function(
    data, 
    query = data, 
    k = 5,
    metric = "euclidean", 
    n_trees = 50, 
    seed = 42,
    ...
) {
  k <- min(nrow(query), k)
  ann <- annoy_build(data, metric = metric, n_trees = n_trees, seed = seed)
  annoy_search(ann, query, k = k, ...)
}

#' @importFrom RcppAnnoy AnnoyAngular AnnoyEuclidean AnnoyHamming AnnoyManhattan
annoy_build <- function(data, metric = "euclidean", n_trees = 50, seed = 42) {
  ann <- switch(
    EXPR = metric,
    "euclidean" = new(AnnoyEuclidean, ncol(data)),
    "angular" = new(AnnoyAngular, ncol(data)),
    "manhattan" = new(AnnoyManhattan, ncol(data)),
    "hamming" = new(AnnoyHamming, ncol(data)),
    stop("Invalid metric '", metric, "'.")
  )
  ann$setSeed(as.integer(seed))
  for (i in seq_len(nrow(data))) {
    ann$addItem(i - 1, data[i, ])
  }
  ann$build(n_trees)
  ann
}

#' @importFrom methods is
annoy_search <- function(index, query, k, search.k = -1, return.dist = TRUE) {
  idx <- matrix(0, nrow(query), k)
  dist <- matrix(0, nrow(query), k)
  for (i in seq_len(nrow(query))) {
    res <- index$getNNsByVectorList(query[i, ], k, search.k, return.dist)
    idx[i, ] <- res$item + 1L
    if (return.dist) {
      if (is(index, "Rcpp_AnnoyAngular")) {
        dist[i, ] <- 0.5 * (res$distance) ^ 2
      } else {
        dist[i, ] <- res$distance
      }
    }
  }
  return(list(idx = idx, dist = dist))
}

#' @importFrom rnndescent rnnd_build rnnd_query
nndescent_knn <- function(
    data, 
    query = data, 
    k = 5,
    metric = "euclidean", 
    n_trees = 50, 
    seed = 42,
    n_threads = 0,
    ...
) {
  metric <- match.arg(metric, .nnd_metrics)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  ann <- rnnd_build(
    data = data, 
    k = k, 
    metric = metric, 
    n_trees = n_trees, 
    n_threads = n_threads,
    ...
  )
  rnnd_query(index = ann, query = query, k = k, n_threads = n_threads)
}

.annoy_metrics <- c("euclidean", "angular", "manhattan", "hamming")

.nnd_metrics <- c(
  'euclidean', 
  'manhattan', 
  'sqeuclidean',
  'chebyshev', 
  'canberra', 
  'braycurtis', 
  'cosine', 
  'correlation', 
  'hamming', 
  'jaccard',
  'dice',
  'matching',
  'russellrao',
  'kulsinski',
  'rogerstanimoto',
  'sokalmichener',
  'sokalsneath',
  'tsss',
  'yule',
  'hellinger'
)

