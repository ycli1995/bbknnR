
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title 
#' Perform batch balanced KNN
#' 
#' @description 
#' Batch balanced KNN, altering the KNN procedure to identify each cell’s top neighbours 
#' in each batch separately instead of the entire cell pool with no accounting for batch. 
#' The nearest neighbours for each batch are then merged to create a final list of neighbours 
#' for the cell. Aligns batches in a quick and lightweight manner.
#' 
#' @param object An object
#' @param ... Arguments passed to other methods
#' 
#' @references 
#' Polański, Krzysztof, et al. "BBKNN: fast batch alignment of single cell 
#' transcriptomes." Bioinformatics 36.3 (2020): 964-965.
#' 
#' @name RunBBKNN
#' 
#' @export
RunBBKNN <- function(object, ...) {
  UseMethod(generic = 'RunBBKNN', object = object)
}

#' @param batch_list A character vector with the same length as nrow(pca)
#' @param n_pcs Number of dimensions to use. Default is 50.
#' @param neighbors_within_batch How many top neighbours to report for each batch; 
#' total number of neighbours in the initial k-nearest-neighbours computation will 
#' be this number times the number of batches. This then serves as the basis for 
#' the construction of a symmetrical matrix of connectivities.
#' @param trim Trim the neighbours of each cell to these many top connectivities. 
#' May help with population independence and improve the tidiness of clustering. 
#' The lower the value the more independent the individual populations, at the 
#' cost of more conserved batch effect. Default is 10 times neighbors_within_batch 
#' times the number of batches. Set to 0 to skip.
#' @param approx If TRUE, use approximate neighbour finding - RcppAnnoy or pyNNDescent. 
#' This results in a quicker run time for large datasets while also potentially 
#' increasing the degree of batch correction.
#' @param use_annoy Only used when approx = TRUE. If TRUE, will use RcppAnnoy for 
#' neighbour finding. If FALSE, will use pyNNDescent instead.
#' @param annoy_n_trees Only used with annoy neighbour identification. The number of 
#' trees to construct in the annoy forest. More trees give higher precision when 
#' querying, at the cost of increased run time and resource intensity.
#' @param pynndescent_n_neighbors Only used with pyNNDescent neighbour identification. 
#' The number of neighbours to include in the approximate neighbour graph. More neighbours 
#' give higher precision when querying, at the cost of increased run time and resource 
#' intensity.
#' @param pynndescent_random_state Only used with pyNNDescent neighbour identification. 
#' The RNG seed to use when creating the graph.
#' @param use_faiss If approx = FALSE and the metric is "euclidean", use the faiss 
#' package to compute nearest neighbours if installed. This improves performance at 
#' a minor cost to numerical precision as faiss operates on float32.
#' @param metric What distance metric to use. The options depend on the choice of 
#' neighbour algorithm. "euclidean", the default, is always available.
#' @param set_op_mix_ratio Pass to 'set_op_mix_ratio' parameter for \code{\link[uwot]{umap}}
#' @param local_connectivity Pass to 'local_connectivity' parameter for \code{\link[uwot]{umap}}
#' @param seed Set a random seed. By default, sets the seed to 42. Setting NULL will 
#' not set a seed.
#' 
#' @import reticulate methods Matrix
#' 
#' @rdname RunBBKNN
#' @export
RunBBKNN.default <- function(
  object,
  batch_list, 
  n_pcs = 50L,
  neighbors_within_batch = 3L,
  trim = NULL,
  approx = TRUE,
  use_annoy = TRUE, 
  annoy_n_trees = 10L,
  pynndescent_n_neighbors = 30L, 
  pynndescent_random_state = 0L, 
  use_faiss = TRUE, 
  metric = "euclidean",
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  seed = 42,
  ...
) {
  batches <- unique(x = batch_list)
  trim <- trim %||% neighbors_within_batch * length(x = batches) * 10
  
  if (approx && use_annoy) {
    cat("Running BBKNN using RcppAnnoy...\n")
    out <- bbknn_annoy(
      pca = object, 
      batch_list = batch_list,
      neighbors_within_batch = neighbors_within_batch, 
      n_pcs = n_pcs, 
      trim = trim,
      metric = metric,
      annoy_n_trees = annoy_n_trees,
      set_op_mix_ratio = set_op_mix_ratio,
      local_connectivity = local_connectivity
    )
  } else {
    cat("Checking python modules...\n")
    if (!py_module_available(module = "bbknn")){
      stop("Cannot find bbknn, please install through pip (e.g. pip install bbknn).")
    }
    cat("Running BBKNN via python...\n")
    bbknn <- import(module = "bbknn", convert = FALSE)
    bbknn_out <- bbknn$bbknn_matrix(
      object, 
      batch_list,
      neighbors_within_batch = as.integer(x = neighbors_within_batch), 
      n_pcs = as.integer(x = n_pcs), 
      trim = as.integer(x = trim),
      approx = approx, 
      annoy_n_trees = as.integer(x = annoy_n_trees), 
      pynndescent_n_neighbors = as.integer(x = pynndescent_n_neighbors), 
      pynndescent_random_state = pynndescent_random_state, 
      use_annoy = FALSE, 
      use_faiss = use_faiss, 
      metric = metric,
      set_op_mix_ratio = set_op_mix_ratio,
      local_connectivity = local_connectivity
    )
    cat("Getting BBKNN Graph and distances...\n")
    bbknn_cnts <- new(
      Class = "dgCMatrix",
      i = as.vector(py_to_r(bbknn_out[[1]]$indices)),
      p = as.vector(py_to_r(bbknn_out[[1]]$indptr)),
      x = as.vector(py_to_r(bbknn_out[[1]]$data)),
      Dim = unlist(x = py_to_r(bbknn_out[[1]]$shape))
    )
    bbknn_dist <- new(
      Class = "dgCMatrix",
      i = as.vector(py_to_r(bbknn_out[[0]]$indices)),
      p = as.vector(py_to_r(bbknn_out[[0]]$indptr)),
      x = as.vector(py_to_r(bbknn_out[[0]]$data)),
      Dim = unlist(x = py_to_r(bbknn_out[[0]]$shape))
    )
    out <- list(cnts = bbknn_cnts, dist = bbknn_dist)
  }
  if (!is.null(x = rownames(x = object))) {
    rownames(x = out$cnts) <- colnames(x = out$cnts) <- rownames(x = object)
    rownames(x = out$dist) <- colnames(x = out$dist) <- rownames(x = object)
  }
  cat("BBKNN finish!\n")
  return(out)
}

#' @param batch_key Column name in meta.data discriminating between your batches.
#' @param assay used to construct Graph.
#' @param reduction Which dimensional reduction to use for the BBKNN input. 
#' Default is PCA
#' @param graph_name Name of the generated BBKNN graph. Default is bbknn.
#' @param run_TSNE Whether or not to run t-SNE based on BBKNN results.
#' @param TSNE_name Name to store t-SNE dimensional reduction.
#' @param TSNE_key Specifies the string before the number of the t-SNE dimension 
#' names. tSNE by default.
#' @param run_UMAP Whether or not to run UMAP based on BBKNN results.
#' @param UMAP_name Name to store UMAP dimensional reduction.
#' @param UMAP_key Specifies the string before the number of the UMAP dimension 
#' names. tSNE by default.
#' @param min_dist Pass to 'min_dist' parameter for \code{\link[uwot]{umap}}
#' @param spread Pass to 'spread' parameter for \code{\link[uwot]{umap}}
#' 
#' @return Returns a Seurat object containing a new BBKNN Graph. If run t-SNE or UMAP, will 
#' also return corresponded reduction objects.
#' 
#' @import Seurat SeuratObject uwot Rtsne
#' 
#' @rdname RunBBKNN
#' @export
RunBBKNN.Seurat <- function(
  object,
  batch_key,
  assay = NULL,
  reduction = "pca",
  n_pcs = 50L,
  graph_name = "bbknn",
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  run_TSNE = TRUE,
  TSNE_name = "tsne",
  TSNE_key = "tSNE_",
  run_UMAP = TRUE,
  UMAP_name = "umap",
  UMAP_key = "UMAP_",
  min_dist = 0.3,
  spread = 1,
  seed = 42,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  embeddings <- object[[reduction]]@cell.embeddings
  batch_list <- as.character(x = object@meta.data[, batch_key])
  bbknn_out <- RunBBKNN(
    embeddings,
    batch_list,
    n_pcs = n_pcs,
    seed = seed,
    ...
  )
  bbknn_out$cnts <- as.Graph(x = bbknn_out$cnts)
  slot(object = bbknn_out$cnts, name = "assay.used") <- assay
  object[[graph_name]] <- bbknn_out$cnts
  if (run_TSNE || run_UMAP) {
    bbknn <- retrive_knn(dist = bbknn_out$dist)
    if (!is.null(x = seed)) {
      set.seed(seed = seed)
    }
    if (run_UMAP) {
      cat("Running UMAP...\n")
      umap <- umap(
        X = embeddings,
        n_neighbors = ncol(x = bbknn$idx),
        nn_method = bbknn,
        set_op_mix_ratio = set_op_mix_ratio,
        local_connectivity = local_connectivity,
        min_dist = min_dist,
        spread = spread
      )
      colnames(x = umap) <- paste0(UMAP_key, c(1, 2))
      rownames(x = umap) <- rownames(x = embeddings)
      object[[UMAP_name]] <- CreateDimReducObject(
        embeddings = umap,
        assay = assay,
        key = UMAP_key
      )
    }
    if (run_TSNE) {
      cat("Running tSNE...\n")
      perplexity <- ncol(x = bbknn$idx) - 1
      tsne <- Rtsne_neighbors(
        index = bbknn$idx,
        distance = bbknn$dist,
        perplexity = perplexity
      )$Y
      colnames(x = tsne) <- paste0(TSNE_key, c(1, 2))
      rownames(x = tsne) <- rownames(x = embeddings)
      object[[TSNE_name]] <- CreateDimReducObject(
        embeddings = tsne,
        assay = assay,
        key = TSNE_key
      )
    }
  }
  cat("All done!\n")
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Perform batch balanced KNN (BBKNN) using RcppAnnoy
#
# @param pca Matrix of input dimensional reduction
# @param batch_list A character vector with the same length as nrow(pca)
# @param neighbors_within_batch How many top neighbours to report for each batch
# @param n_pcs Number of dimensions to use. Default is 50
# @param trim Trim the neighbours of each cell to these many top connectivities.
# @param metric Distance metric for annoy. Options include: euclidean, cosine, 
# manhattan, and hamming
# @param annoy_n_trees More trees gives higher precision when using annoy approximate 
# nearest neighbor search
# @param set_op_mix_ratio Interpolate between (fuzzy) union and intersection as the set 
# operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial 
# sets. Both fuzzy set operations use the product t-norm. The value of this parameter should 
# be between 0.0 and 1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use 
# a pure fuzzy intersection.
# @param local_connectivity The local connectivity required – i.e. the number of nearest 
# neighbors that should be assumed to be connected at a local level. The higher this value 
# the more connected the manifold becomes locally. In practice this should be not more than 
# the local intrinsic dimension of the manifold.
# @param seed Set a random seed. By default, sets the seed to 42. Setting NULL will not set 
# a seed.
# 
# @return Returns a list containing sparse matrices for distances and connectivities. 
# The connectivities are the actual neighbourhood graph.
# 
#' @import Rcpp RcppAnnoy
#'
bbknn_annoy <- function(
  pca,
  batch_list,
  neighbors_within_batch = 3, 
  n_pcs = 50, 
  trim = NULL,
  metric = c("euclidean", "angular", "manhattan", "hamming"),
  annoy_n_trees = 10L,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  seed = 42
) {
  metric <- match.arg(arg = metric)
  batch_list <- as.character(x = batch_list)
  batches <- unique(x = batch_list)
  k <- neighbors_within_batch * length(x = batches)
  trim <- trim %||% 10 * k
  if (min(table(batch_list)) < neighbors_within_batch) {
    stop("Not all batches have at least 'neighbors_within_batch' cells in them.")
  }
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  knn_dist <- matrix(0, nrow = nrow(x = pca), ncol = k)
  knn_index <- matrix(0L, nrow = nrow(x = pca), ncol = k)
  for (to_ind in seq(length(x = batches))) {
    batch_to <- batches[to_ind]
    mask_to <- batch_list == batch_to
    ind_to <- seq(length(x = batch_list))[mask_to]
    data <- pca[mask_to, 1:n_pcs]
    ckd <- switch(
      EXPR = metric,
      "euclidean" = new(AnnoyEuclidean, ncol(x = data)),
      "angular" = new(AnnoyAngular, ncol(x = data)),
      "manhattan" = new(AnnoyManhattan, ncol(x = data)),
      "hamming" = new(AnnoyHamming, ncol(x = data)),
    )
    for (i in seq(nrow(x = data))) {
      ckd$addItem(i - 1, data[i, ])
    }
    ckd$build(annoy_n_trees)
    for (from_ind in seq(length(x = batches))) {
      batch_from <- batches[from_ind]
      mask_from <- batch_list == batch_from
      ind_from <- seq(length(x = batch_list))[mask_from]
      ckdout <- query_annoy_tree(
        query = pca[mask_from, 1:n_pcs],
        ckd = ckd,
        k = neighbors_within_batch
      )
      ckdout$index <- GetRawIndex(ckdout$index, ind_to)
      col_range <- seq((to_ind - 1) * neighbors_within_batch + 1, to_ind * neighbors_within_batch)
      knn_index <- UpdateIntMat(knn_index, ckdout$index, ind_from - 1L, col_range - 1L)
      knn_dist <- UpdateNumMat(knn_dist, ckdout$dist, ind_from - 1L, col_range - 1L)
    }
    ckd$unload()
  }
  for (i in seq(nrow(x = knn_index))) {
    knn_index[i, ] <- knn_index[i, ][order(knn_dist[i, ])]
    knn_dist[i, ] <- knn_dist[i, ][order(knn_dist[i, ])]
  }
  out <- compute_connectivities_umap(
    knn_index = knn_index,
    knn_dist = knn_dist,
    n_obs = nrow(x = knn_index),
    n_neighbors = ncol(x = knn_index),
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    seed = seed
  )
  if (trim > 0) {
    out$cnts <- trimming(out$cnts, trim = trim)
  }
  return(out)
}

# Query k nearest neighbors from AnnoyIndex
# 
# @param query Query data
# @param ckd Pre-computed AnnoyIndex
# @param k k nearest neighbors
# 
#' @import RcppAnnoy
#' 
query_annoy_tree <- function(query, ckd, k = 3) {
  index <- matrix(nrow = nrow(x = query), ncol = k)
  dist <- matrix(nrow = nrow(x = query), ncol = k)
  for (i in seq(nrow(x = query))) {
    holder <- ckd$getNNsByVectorList(
      query[i, ], 
      n = k,
      search_k = -1, 
      include_distances = TRUE
    )
    index[i, ] <- holder$item
    dist[i, ] <- holder$distance
  }
  return(list(index = index + 1L, dist = dist))
}

# fuzzy_simplicial_set
#
# Export fuzzy_simplicial_set() function from \code{uwot}
#
#' @import uwot
#' @importFrom utils getFromNamespace
#' 
fuzzy_simplicial_set <- getFromNamespace(x = "fuzzy_simplicial_set", ns = "uwot")

# Compute connectivities using UMAP
# 
# @param knn_index KNN index matrix
# @param knn_dist KNN distances matrix
# @param set_op_mix_ratio Interpolate between (fuzzy) union and intersection as the 
# set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy 
# simplicial sets. Both fuzzy set operations use the product t-norm.
# @param local_connectivity The local connectivity required – i.e. the number of 
# nearest neighbors that should be assumed to be connected at a local level. 
# @param seed Set a random seed. 
#
#' @import Matrix methods uwot dplyr
#' @importFrom Rcpp evalCpp
#'
compute_connectivities_umap <- function(
  knn_index,
  knn_dist,
  n_obs,
  n_neighbors,
  set_op_mix_ratio = 1.0,
  local_connectivity = 1.0,
  seed = 42
) {
  X <- matrix(0, nrow = n_obs, ncol = 1)
  if (!is.null(x = seed)) {
    set.seed(seed)
  }
  cnts <- fuzzy_simplicial_set(
    nn = list(idx = knn_index, dist = knn_dist),
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity
  )
  dist <- do.call(
    sparseMatrix, 
    c(GetSparseDist(knn_index, knn_dist, n_obs, n_neighbors), repr = "T", index1 = FALSE)
  ) %>%
    drop0() %>%
    as(Class = "dgCMatrix")
  return(list(dist = dist, cnts = cnts))
}

# Trim the connectivities
# 
# @param cnts A weighted adjacency dgCMatrix of the neighborhood graph of data points.
# @param trim Trim the neighbours of each cell to these many top connectivities. 
#
#' @importFrom Matrix drop0
#' @importFrom Rcpp evalCpp
#'
trimming <- function(cnts, trim) {
  Trimming(cnts@x, cnts@i, cnts@p, trim)
  cnts <- drop0(cnts)
  return(cnts)
}

# Retrive knn index and distances from BBKNN output
#  
# @param dist A dgCMatrix of pairwise distances
#
retrive_knn <- function(dist) {
  ncol <- dist@p[2] - dist@p[1]
  nrow <- ncol(x = dist)
  knn_index <- matrix(0L, nrow = nrow, ncol = ncol)
  knn_dist <- matrix(0, nrow = nrow, ncol = ncol)
  for (i in 1:nrow) {
    knn_dist[i, ] <- dist@x[((i - 1) * ncol + 1):(i * ncol)][order(dist@x[((i - 1) * ncol + 1):(i * ncol)])]
    knn_index[i, ] <- dist@i[((i - 1) * ncol + 1):(i * ncol)][order(dist@x[((i - 1) * ncol + 1):(i * ncol)])] + 1L
  }
  knn_dist <- cbind(rep(x = 0, nrow), knn_dist)
  knn_index <- cbind(1:nrow, knn_index)
  return(list(idx = knn_index, dist = knn_dist))
}

