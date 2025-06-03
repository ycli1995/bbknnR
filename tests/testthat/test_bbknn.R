
test_that("BBKNN works for matrix", {
  data("panc8_small")
  pca <- Embeddings(panc8_small, reduction = "pca")
  b <- panc8_small$dataset
  n <- length(unique(b))
  k <- 3
  
  bbknn <- RunBBKNN(pca, b, neighbors_within_batch = k, verbose = FALSE)
  expect_equal(nrow(bbknn$nn.idx), nrow(pca))
  expect_equal(nrow(bbknn$nn.dist), nrow(pca))
  expect_equal(ncol(bbknn$nn.idx), k * n)
  
  expect_equal(nrow(bbknn$connectivities), nrow(pca))
  expect_equal(ncol(bbknn$connectivities), nrow(pca))
})

test_that("BBKNN works for nndescent", {
  data("panc8_small")
  pca <- Embeddings(panc8_small, reduction = "pca")
  b <- panc8_small$dataset
  n <- length(unique(b))
  k <- 3
  
  bbknn <- RunBBKNN(
    pca, b, 
    neighbors_within_batch = k, 
    method = "nndescent",
    verbose = FALSE
  )
  expect_equal(nrow(bbknn$nn.idx), nrow(pca))
  expect_equal(nrow(bbknn$nn.dist), nrow(pca))
  expect_equal(ncol(bbknn$nn.idx), k * n)
  
  expect_equal(nrow(bbknn$connectivities), nrow(pca))
  expect_equal(ncol(bbknn$connectivities), nrow(pca))
})

test_that("BBKNN works for Seurat", {
  data("panc8_small")
  
  k <- 3
  panc8_small <- RunBBKNN(
    panc8_small, 
    batch_key = "tech", 
    neighbors_within_batch = k, 
    method = "nndescent",
    verbose = FALSE
  )
  expect_s4_class(panc8_small[['RNA_bbknn']], "Graph")
  expect_s4_class(panc8_small[['RNA.bbknn']], "Neighbor")
  expect_s4_class(panc8_small[['umap']], "DimReduc")
  expect_s4_class(panc8_small[['tsne']], "DimReduc")
})

test_that("BBKNN returns umap model", {
  data("panc8_small")
  
  k <- 3
  expect_warning(panc8_small <- RunBBKNN(
    panc8_small, 
    batch_key = "tech", 
    neighbors_within_batch = k, 
    method = "nndescent",
    return.umap.model = TRUE,
    verbose = FALSE
  ))
  expect_true(is.list(Misc(panc8_small[['umap']], "model")))
})



