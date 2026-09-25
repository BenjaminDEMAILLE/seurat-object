set.seed(42)

# A UMAP built with return.model = TRUE keeps the embedding it was fit on in
# misc$model. Subsetting left that untouched, so the model and the reduction
# described different sets of cells, and a query projected onto the reduction
# was placed against the wrong coordinates
build_reduction <- function(ncell = 40L, ndim = 3L, with_model = TRUE) {
  cells <- paste0("c", seq_len(ncell))
  embeddings <- matrix(
    rnorm(ncell * ndim),
    nrow = ncell,
    dimnames = list(cells, paste0("umap_", seq_len(ndim)))
  )
  reduction <- CreateDimReducObject(embeddings = embeddings, key = "umap_", assay = "RNA")
  # a model of the shape uwot stores
  if (with_model) {
    Misc(reduction, slot = "model") <- list(
      embedding = embeddings,
      n_neighbors = 5L,
      metric = list(euclidean = list())
    )
  }
  reduction
}

test_that("subsetting a reduction subsets the model it carries", {
  reduction <- build_reduction()
  kept <- Cells(reduction)[seq(from = 11, to = 40)]
  subset.reduction <- expect_no_warning(subset(reduction, cells = kept))

  expect_identical(Cells(subset.reduction), kept)
  expect_identical(rownames(Misc(subset.reduction, slot = "model")$embedding), kept)
  # and the coordinates themselves are the ones the reduction has
  expect_identical(
    Misc(subset.reduction, slot = "model")$embedding,
    Embeddings(subset.reduction)
  )
})

test_that("renaming cells keeps the stored model aligned through subsetting", {
  reduction <- build_reduction()
  old.cells <- Cells(reduction)
  new.cells <- paste0("renamed_", old.cells)
  renamed <- expect_no_warning(RenameCells(reduction, new.names = new.cells))
  expect_identical(rownames(Misc(renamed, slot = "model")$embedding), new.cells)
  expect_identical(Misc(renamed, slot = "model")$embedding, Embeddings(renamed))

  kept <- new.cells[seq(from = 11, to = 40)]
  subset.reduction <- expect_no_warning(subset(renamed, cells = kept))
  expect_identical(rownames(Misc(subset.reduction, slot = "model")$embedding), kept)
  expect_identical(Misc(subset.reduction, slot = "model")$embedding, Embeddings(subset.reduction))
})

test_that("the rest of the model is left alone", {
  reduction <- build_reduction()
  subset.reduction <- subset(reduction, cells = Cells(reduction)[1:20])
  model <- Misc(subset.reduction, slot = "model")
  expect_identical(model$n_neighbors, 5L)
  expect_identical(model$metric, list(euclidean = list()))
})

test_that("a reduction with no model is unchanged", {
  reduction <- build_reduction(with_model = FALSE)
  subset.reduction <- subset(reduction, cells = Cells(reduction)[1:20])
  expect_null(Misc(subset.reduction, slot = "model"))
  expect_identical(nrow(Embeddings(subset.reduction)), 20L)

  renamed <- expect_no_warning(RenameCells(reduction, new.names = paste0("renamed_", Cells(reduction))))
  expect_null(Misc(renamed, slot = "model"))
  expect_identical(Cells(renamed), paste0("renamed_", Cells(reduction)))
})

test_that("subsetting refuses a model that cannot match every retained cell", {
  reduction <- build_reduction()
  model <- Misc(reduction, slot = "model")
  model$embedding <- model$embedding[-1, , drop = FALSE]
  slot(reduction, name = "misc")[["model"]] <- model

  expect_error(
    subset(reduction, cells = Cells(reduction)[1:20]),
    "Cannot align stored reduction model embedding"
  )
})
