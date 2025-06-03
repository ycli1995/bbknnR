
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
#' @name RunBBKNN
#' @export RunBBKNN
RunBBKNN <- function(object, ...) UseMethod('RunBBKNN', object)

#' @param batch_list A character vector with the same length as nrow(pca)
#' @param neighbors_within_batch How many top neighbours to report for each 
#' batch; total number of neighbours in the initial k-nearest-neighbours 
#' computation will be this number times the number of batches. This then serves 
#' as the basis for the construction of a symmetrical matrix of connectivities.
#' @param n_pcs Number of dimensions to use. Default is 50.
#' @param method Method to find k nearest neighbors (kNNs). One of "annoy" and 
#' "nndescent".
#' @param metric Metric to calculate the distances. The options depend on the 
#' choice of kNN \code{method}. The following metrics are supported in both 
#' \code{annoy} and \code{nndescent}:
#' \itemize{
#' \item 'euclidean' (the default)
#' \item 'manhattan'
#' \item 'hamming'
#' }
#' The following metrics are only supported in \code{nndescent}:
#' \itemize{
#' \item 'sqeuclidean'
#' \item 'chebyshev'
#' \item 'canberra'
#' \item 'braycurtis' 
#' \item 'cosine' 
#' \item 'correlation'
#' \item 'jaccard' 
#' \item 'dice'
#' \item 'matching'
#' \item 'russellrao' 
#' \item 'kulsinski' 
#' \item 'rogerstanimoto' 
#' \item 'sokalmichener'
#' \item 'sokalsneath'
#' \item 'tsss'
#' \item 'yule'
#' \item 'hellinger'
#' }
#' @param n_trees The number of trees to use in the random projection forest. 
#' More trees give higher precision when querying, at the cost of increased run 
#' time and resource intensity.
#' @param trim Trim the neighbours of each cell to these many top 
#' connectivities. May help with population independence and improve the 
#' tidiness of clustering. The lower the value the more independent the 
#' individual populations, at the cost of more conserved batch effect. Default 
#' is 10 times neighbors_within_batch times the number of batches. Set to 0 to 
#' skip.
#' @param set_op_mix_ratio Pass to 'set_op_mix_ratio' parameter for 
#' \code{\link[uwot]{umap}}
#' @param local_connectivity Pass to 'local_connectivity' parameter for 
#' \code{\link[uwot]{umap}}
#' @param seed Set a random seed. By default, sets the seed to 42. Setting 
#' \code{NULL} will not set a seed.
#' @param verbose Whether or not to print output to the console
#' 
#' @rdname RunBBKNN
#' @export
#' @method RunBBKNN matrix
RunBBKNN.matrix <- function(
    object,
    batch_list,
    neighbors_within_batch = 3, 
    n_pcs = 50, 
    method = c("annoy", "nndescent"),
    metric = "euclidean",
    n_trees = 10L,
    trim = NULL,
    set_op_mix_ratio = 1,
    local_connectivity = 1,
    seed = 42,
    verbose = TRUE,
    ...
) {
  method <- match.arg(method)
  run_bbknn(
    pca = object,
    batch_list = batch_list,
    neighbors_within_batch = neighbors_within_batch, 
    n_pcs = n_pcs, 
    method = method,
    metric = metric,
    n_trees = n_trees,
    trim = trim,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    seed = seed,
    verbose = verbose,
    ...
  )
}

#' @param batch_key Column name in meta.data discriminating between your 
#' batches.
#' @param assay Used to construct Graph.
#' @param reduction Which dimensional reduction to use for the BBKNN input. 
#' Default is PCA
#' @param graph_name Name of the generated BBKNN graph. Default is "bbknn".
#' @param run_TSNE Whether or not to run t-SNE based on BBKNN results.
#' @param TSNE_name Name to store t-SNE dimensional reduction.
#' @param TSNE_key Specifies the string before the number of the t-SNE dimension 
#' names. tSNE by default.
#' @param run_UMAP Whether or not to run UMAP based on BBKNN results.
#' @param UMAP_name Name to store UMAP dimensional reduction.
#' @param UMAP_key Specifies the string before the number of the UMAP dimension 
#' names. tSNE by default.
#' @param return.umap.model Whether UMAP will return the uwot model.
#' @param min_dist Pass to 'min_dist' parameter for \code{\link[uwot]{umap}}
#' @param spread Pass to 'spread' parameter for \code{\link[uwot]{umap}}
#' 
#' @return Returns a Seurat object containing a new BBKNN Graph and Neighbor 
#' data. If run t-SNE or UMAP, will also return corresponded reduction objects.
#' 
#' @examples
#' data("panc8_small")
#' panc8_small <- RunBBKNN(panc8_small, "tech")
#' 
#' @importFrom SeuratObject as.Graph Cells CreateDimReducObject DefaultAssay 
#' Embeddings FetchData Misc Misc<-
#' @importFrom uwot optimize_graph_layout umap
#' @importFrom Rtsne Rtsne_neighbors
#' @importFrom future nbrOfWorkers
#' @importFrom methods new slot<-
#' 
#' @importClassesFrom SeuratObject Graph Neighbor Seurat
#' 
#' @rdname RunBBKNN
#' @export
#' @method RunBBKNN Seurat
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
    return.umap.model = FALSE,
    min_dist = 0.3,
    spread = 1,
    seed = 42,
    verbose = TRUE,
    ...
) {
  assay <- assay %||% DefaultAssay(object)
  embeddings <- Embeddings(object, reduction = reduction)
  batch_list <- FetchData(object, vars = batch_key[1])[, 1]
  bbknn_out <- RunBBKNN(
    embeddings,
    batch_list,
    n_pcs = n_pcs,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    seed = seed,
    verbose = verbose,
    ...
  )
  graph <- as.Graph(bbknn_out$connectivities)
  slot(graph, name = "assay.used") <- assay
  object[[paste0(assay, "_", graph_name)]] <- graph
  
  object[[paste0(assay, ".", graph_name)]] <- new(
    "Neighbor",
    nn.idx = bbknn_out$nn.idx,
    nn.dist = bbknn_out$nn.dist,
    cell.names = Cells(object)
  )
  if (!run_TSNE & !run_UMAP) {
    return(object)
  }
  if (run_UMAP) {
    if (verbose) {
      message("Running UMAP...")
    }
    if (!is.null(seed)) {
      set.seed(seed)
    }
    if (return.umap.model) {
      warning(
        "Trimmed connectivity graph will be ignored ", 
        "since 'return.umap.model = TRUE'",
        call. = FALSE, immediate. = TRUE
      )
      umap <- umap(
        X = embeddings,
        n_threads = nbrOfWorkers(),
        n_neighbors = ncol(bbknn_out$nn.idx),
        nn_method = list(idx = bbknn_out$nn.idx, dist = bbknn_out$nn.dist),
        set_op_mix_ratio = set_op_mix_ratio,
        local_connectivity = local_connectivity,
        min_dist = min_dist,
        spread = spread,
        ret_model = TRUE
      )
      umap.model <- umap
      umap <- umap$embedding
    } else {
      umap <- optimize_graph_layout(
        X = embeddings,
        graph = bbknn_out$connectivities,
        n_components = 2,
        min_dist = min_dist, 
        spread = spread
      )
      umap.model <- NULL
    }
    colnames(umap) <- paste0(UMAP_key, seq_len(ncol(umap)))
    rownames(umap) <- rownames(embeddings)
    object[[UMAP_name]] <- CreateDimReducObject(
      embeddings = umap,
      assay = assay,
      key = UMAP_key
    )
    if (return.umap.model) {
      Misc(object[[UMAP_name]], slot = "model") <- umap.model
    }
  }
  if (run_TSNE) {
    if (verbose) {
      message("Running tSNE...")
    }
    if (!is.null(seed)) {
      set.seed(seed)
    }
    perplexity <- ncol(bbknn_out$nn.idx) - 1
    tsne <- Rtsne_neighbors(
      index = bbknn_out$nn.idx,
      distance = bbknn_out$nn.dist,
      perplexity = perplexity
    )$Y
    colnames(tsne) <- paste0(TSNE_key, seq_len(ncol(umap)))
    rownames(tsne) <- rownames(embeddings)
    object[[TSNE_name]] <- CreateDimReducObject(
      embeddings = tsne,
      assay = assay,
      key = TSNE_key
    )
  }
  return(object)
}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Internal #####################################################################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
run_bbknn <- function(
    pca,
    batch_list,
    neighbors_within_batch = 3, 
    n_pcs = 50, 
    method = c("annoy", "nndescent"),
    metric = "euclidean",
    n_trees = 10L,
    trim = NULL,
    set_op_mix_ratio = 1,
    local_connectivity = 1,
    seed = 42,
    verbose = verbose,
    ...
) {
  batch_list <- as.character(batch_list)
  stopifnot(length(batch_list) == nrow(pca))
  
  if (length(rownames(pca)) == 0) {
    rownames(pca) <- seq_len(nrow(pca))
  }
  
  batch_counts <- table(batch_list)
  batches <- names(batch_counts)
  min_batch <- min(batch_counts)
  if (min_batch < neighbors_within_batch) {
    b <- batches[which.min(batch_counts)]
    stop(
      "Not all batches have at least 'neighbors_within_batch' cells in them: ",
      "\n neighbors_within_batch: ", neighbors_within_batch,
      "\n cells of ", b, ": ", min_batch
    )
  }
  
  method <- match.arg(method)
  if (n_pcs < ncol(pca)) {
    pca <- pca[, seq_len(n_pcs)]
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  metric <- metric[1]
  if (any(!metric %in% c(.annoy_metrics, .nnd_metrics))) {
    stop("Invalid kNN metric: ", metric)
  }
  if (any(!metric %in% .annoy_metrics) & method == "annoy") {
    stop("Invalid kNN metric for 'annoy': ", metric)
  }
  if (any(!metric %in% .nnd_metrics) & method == "nndescent") {
    stop("Invalid kNN metric for 'nndescent': ", metric)
  }
  knn_func <- switch(method, annoy = annoy_knn, nndescent = nndescent_knn)
  
  my.lapply <- lapply
  if (nbrOfWorkers() > 1) {
    my.lapply <- future_lapply
  }
  if (verbose) {
    message(
      "Find BBKNN for each batch: \n",
      "  neighbors_within_batch: ", neighbors_within_batch, "\n",
      "  n_pcs: ", n_pcs, "\n",
      "  method: '", method, "'\n",
      "  metric: '", metric, "'"
    )
  }
  all.nn <- my.lapply(batches, function(b) {
    data.idx <- which(batch_list == b)
    data <- pca[data.idx, ]
    nn <- knn_func(
      data = data, 
      query = pca, 
      k = neighbors_within_batch, 
      metric = metric, 
      n_trees = n_trees, 
      ...
    )
    nn$idx <- correct_idx_cpp(nn$idx, data.idx)
    nn
  })
  nn.idx <- Reduce(cbind, lapply(all.nn, `[[`, 1))
  nn.dist <- Reduce(cbind, lapply(all.nn, `[[`, 2))
  for (i in seq_len(nrow(nn.idx))) {
    nn.idx[i, ] <- nn.idx[i, ][order(nn.dist[i, ])]
    nn.dist[i, ] <- nn.dist[i, ][order(nn.dist[i, ])]
  }
  trim <- trim %||% 10 * ncol(nn.idx)
  if (verbose) {
    message("Compute connectivity graph with 'trim = ", trim, "'")
  }
  cnts <- get_connectivities_umap(
    nn.idx = nn.idx,
    nn.dist = nn.dist,
    set_op_mix_ratio = 1.0,
    local_connectivity = 1.0,
    trim = trim,
    seed = seed
  )
  rownames(cnts) <- colnames(cnts) <- rownames(pca)
  rownames(nn.idx) <- rownames(nn.dist) <- rownames(pca)
  if (verbose) {
    message("BBKNN finish!")
  }
  list(nn.idx = nn.idx, nn.dist = nn.dist, connectivities = cnts)
}

#' @importFrom uwot similarity_graph
get_connectivities_umap <- function(
    nn.idx,
    nn.dist,
    set_op_mix_ratio = 1.0,
    local_connectivity = 1.0,
    trim = NULL,
    seed = 42
) {
  trim <- trim %||% 10 * nrow(nn.idx)
  X <- matrix(0, nrow = nrow(nn.idx), ncol = 2)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  cnts <- similarity_graph(
    X = X,
    nn_method = list(idx = nn.idx, dist = nn.dist),
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity
  )
  if (trim > 0) {
    cnts <- trim_graph(cnts, as.integer(trim))
  }
  cnts
}


