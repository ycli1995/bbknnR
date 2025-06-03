
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
#' @name RidgeRegression
#' @export RidgeRegression
RidgeRegression <- function(object, ...) UseMethod('RidgeRegression', object)

#' @param latent_data Extra data to regress out, should be cells x latent data
#' @param batch_key Variables to regress out as technical effects. Must be 
#' included in column names of latent_data
#' @param confounder_key Variables to to retain as biological effects. Must be 
#' included in column names of latent_data
#' @param lambda A user supplied lambda sequence. pass to 
#' \code{\link[glmnet]{glmnet}}
#' @param seed Set a random seed. By default, sets the seed to 42. Setting NULL 
#' will not set a seed.
#' @param verbose Whether or not to print output to the console
#' 
#' @importFrom glmnet glmnet coef.glmnet
#' @importFrom tidytable get_dummies
#' 
#' @rdname RidgeRegression
#' @export
#' @method RidgeRegression default
RidgeRegression.default <- function(
  object,
  latent_data,
  batch_key,
  confounder_key,
  lambda = 1,
  seed = 42,
  verbose = TRUE,
  ...
) {
  latent_data <- latent_data[, c(batch_key, confounder_key), drop = FALSE]
  if (verbose) {
    message("Running ridge regression...")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  mod <- get_dummies(.df = latent_data)
  mod <- as.matrix(mod[, (ncol(latent_data) + 1):ncol(mod)])
  batch_idx <- sapply(
    X = batch_key,
    FUN = function(x) startsWith(x = colnames(x = mod), prefix = paste0(x, "_"))
  )
  batch_idx <- Reduce(xor, batch_idx)
  mod_tech <- mod[, batch_idx]
  resid_data <- matrix(0, nrow = nrow(object), ncol = ncol(object))
  for (i in 1:nrow(x = object)) {
    fit <- glmnet(
      x = mod, 
      y = object[i, ], 
      alpha = 0, 
      intercept = FALSE, 
      standardize = FALSE,
      lambda = lambda,
      ...
    )
    cf <- as.matrix(coef.glmnet(fit))[-1, ]
    cf <- cf[batch_idx]
    explained <- mod_tech %*% cf
    regression.mat <- object[i, ] - explained
    resid_data[i, ] <- regression.mat
  }
  rownames(resid_data) <- rownames(object)
  colnames(resid_data) <- colnames(object)
  return(resid_data)
}

#' @param assay Name of Assay ridge regression is being run on
#' @param features Features to compute ridge regression on. If features=NULL, 
#' ridge regression will be run using the variable features for the Assay. 
#' @param run_pca Whether or not to run pca with regressed expression data 
#' (TRUE by default)
#' @param npcs Total Number of PCs to compute and store (50 by default)
#' @param reduction.name Dimensional reduction name (pca by default)
#' @param reduction.key Dimensional reduction key, specifies the string before 
#' the number for the dimension names (PC by default)
#' @param replace Whether or not to replace original scale.data with regressed 
#' expression data (TRUE by default)
#' 
#' @return Returns a Seurat object.
#' 
#' @examples
#' data("panc8_small")
#' panc8_small <- RidgeRegression(panc8_small, "tech", c("nCount_RNA"))
#' 
#' @importFrom Seurat RunPCA
#' @importFrom SeuratObject DefaultAssay FetchData GetAssayData SetAssayData 
#' VariableFeatures
#' @importFrom utils packageVersion
#' @importClassesFrom SeuratObject Seurat
#' 
#' @rdname RidgeRegression
#' @export
#' @method RidgeRegression Seurat
RidgeRegression.Seurat <- function(
  object,
  batch_key,
  confounder_key,
  assay = NULL,
  features = NULL,
  lambda = 1,
  run_pca = TRUE,
  npcs = 50,
  reduction.name = "pca",
  reduction.key = "PC_",
  replace = FALSE,
  seed = 42,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  features <- features %||% VariableFeatures(object, assay = assay)
  if (!any(c(run_pca, replace))) {
    warning(
      "At least one of 'run_pca' or 'replace' should be set up.\n", 
      "Return the original object", 
      call. = FALSE, immediate. = TRUE
    )
    return(object)
  }
  if (packageVersion("Seurat") >= package_version("5.0.0")) {
    data.expr <- GetAssayData(
      object = object, 
      assay = assay, 
      layer = "scale.data"
    )
  } else {
    data.expr <- GetAssayData(
      object = object, 
      assay = assay, 
      slot = "scale.data"
    )
  }
  features <- intersect(features, rownames(data.expr))
  data.expr <- data.expr[features, , drop = FALSE]
  latent_data <- FetchData(object, vars = c(batch_key, confounder_key))
  data.resid <- RidgeRegression(
    object = data.expr,
    latent_data = latent_data,
    batch_key = batch_key,
    confounder_key = confounder_key,
    lambda = lambda,
    seed = seed,
    verbose = verbose,
    ...
  )
  if (replace) {
    if (verbose) {
      message("Replace original scale.data...")
    }
    if (packageVersion("Seurat") >= package_version("5.0.0")) {
      object <- SetAssayData(
        object = object, 
        assay = assay, 
        layer = "scale.data", 
        new.data = data.resid
      )
    } else {
      object <- SetAssayData(
        object = object, 
        assay = assay, 
        slot = "scale.data", 
        new.data = data.resid
      )
    }
  }
  if (run_pca) {
    if (verbose) {
      message("Running PCA...")
    }
    pca <- RunPCA(
      object = data.resid, 
      assay = assay,
      npcs = npcs, 
      reduction.key = reduction.key,
      verbose = FALSE
    )
    object[[reduction.name]] <- pca
  }
  return(object)
}
