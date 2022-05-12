
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title 
#' Perform ridge regression on scaled expression data
#' 
#' @description 
#' Perform ridge regression on scaled expression data, accepting both technical 
#' and biological categorical variables. The effect of the technical variables 
#' is removed while the effect of the biological variables is retained. This is 
#' a preprocessing step that can aid BBKNN integration.
#' 
#' @param object An object
#' @param ... Arguments passed to other methods
#' 
#' @references Park, Jong-Eun, et al. "A cell atlas of human thymic development defines T 
#' cell repertoire formation." Science 367.6480 (2020): eaay3224.
#' 
#' @name RidgeRegression
#' 
#' @export
RidgeRegression <- function(object, ...) {
  UseMethod(generic = 'RidgeRegression', object = object)
}

#' @param latent_data Extra data to regress out, should be cells x latent data
#' @param batch_key Variables to regress out as technical effects. Must be included in column 
#' names of latent_data
#' @param confounder_key Variables to to retain as biological effects. Must be included in 
#' column names of latent_data
#' @param lambda A user supplied lambda sequence. pass to \code{\link[glmnet]{glmnet}}
#' @param seed Set a random seed. By default, sets the seed to 42. Setting NULL will not set 
#' a seed.
#' 
#' @import glmnet tidytable
#' 
#' @rdname RidgeRegression
#' @method RidgeRegression default
#' @export
RidgeRegression.default <- function(
  object,
  latent_data,
  batch_key,
  confounder_key,
  lambda = 1,
  seed = 42,
  ...
) {
  latent_data <- latent_data[, c(batch_key, confounder_key), drop = FALSE]
  cat("Running ridge regression...\n")
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  mod <- get_dummies.(.df = latent_data)
  mod <- as.matrix(x = mod[, (ncol(latent_data) + 1):ncol(mod)])
  batch_idx <- do.call(
    `|`, lapply(
      X = batch_key, 
      FUN = function(x) startsWith(x = colnames(x = mod), prefix = x)
    )
  )
  mod_tech <- mod[, batch_idx]
  resid_data <- matrix(
    nrow = nrow(x = object),
    ncol = ncol(x = object)
  )
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
    cf <- as.matrix(x = coef.glmnet(fit))[-1, ]
    cf <- cf[batch_idx]
    explained <- mod_tech %*% cf
    regression.mat <- object[i, ] - explained
    resid_data[i, ] <- regression.mat
  }
  rownames(x = resid_data) <- rownames(x = object)
  colnames(x = resid_data) <- colnames(x = object)
  return(resid_data)
}

#' @param assay Name of Assay ridge regression is being run on
#' @param features Features to compute ridge regression on. If features=NULL, ridge regression 
#' will be run using the variable features for the Assay. 
#' @param run_pca Whether or not to run pca with regressed expression data (TRUE by default)
#' @param npcs Total Number of PCs to compute and store (50 by default)
#' @param reduction.name Dimensional reduction name (pca by default)
#' @param reduction.key Dimensional reduction key, specifies the string before the number for the 
#' dimension names (PC by default)
#' @param replace Whether or not to replace original scale.data with regressed expression data (TRUE by default)
#' 
#' @return Returns a Seurat object.
#' 
#' @import Seurat SeuratObject
#' 
#' @rdname RidgeRegression
#' @method RidgeRegression Seurat
#' @export
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
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  features <- features %||% VariableFeatures(object = object, assay = assay)
  if (!any(c(run_pca, replace))) {
    warning("At least one of 'run_pca' or 'replace' should be set up.\nReturn the original object", immediate. = TRUE)
    return(object)
  }
  data.expr <- GetAssayData(object = object, assay = assay, slot = "scale.data")
  features <- intersect(x = features, y = rownames(x = data.expr))
  data.expr <- data.expr[features, ]
  latent_data <- object@meta.data[, c(batch_key, confounder_key), drop = FALSE]
  data.resid <- RidgeRegression(
    object = data.expr,
    latent_data = latent_data,
    batch_key = batch_key,
    confounder_key = confounder_key,
    lambda = lambda,
    seed = seed,
    ...
  )
  if (replace) {
    cat("Replace original scale.data...\n")
    object <- SetAssayData(object = object, assay = assay, slot = "scale.data", new.data = data.resid)
  }
  if (run_pca) {
    cat("Running PCA...\n")
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
