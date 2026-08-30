set.seed(42)

# Build a Seurat object whose counts layer is a BPCells matrix on disk
make_ondisk_object <- function(nfeatures = 40, ncells = 25) {
  counts <- matrix(
    data = rpois(n = nfeatures * ncells, lambda = 1.5),
    nrow = nfeatures,
    ncol = ncells
  )
  dimnames(counts) <- list(
    paste0("gene", seq_len(nfeatures)),
    paste0("cell", seq_len(ncells))
  )
  store <- tempfile("bpcells-store")
  suppressWarnings(BPCells::write_matrix_dir(mat = as.sparse(counts), dir = store))
  object <- suppressWarnings(
    CreateSeuratObject(counts = BPCells::open_matrix_dir(dir = store))
  )
  list(object = object, counts = counts)
}

# Save `object` into a fresh directory, then rename that directory, mimicking
# an object being shared or reorganized after it was written
save_and_move <- function(object, relative) {
  saved <- tempfile("saved")
  dir.create(saved)
  suppressWarnings(
    SaveSeuratRds(
      object = object,
      file = file.path(saved, "object.Rds"),
      move = TRUE,
      relative = relative
    )
  )
  moved <- tempfile("moved")
  expect_true(file.rename(from = saved, to = moved))
  file.path(moved, "object.Rds")
}

test_that("LoadSeuratRds finds relative layer paths from any directory", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("fs")
  built <- make_ondisk_object()
  rds <- save_and_move(object = built$object, relative = TRUE)
  # Loading from an unrelated working directory used to resolve the relative
  # paths against the working directory and drop every layer
  elsewhere <- tempfile("elsewhere")
  dir.create(elsewhere)
  old <- setwd(elsewhere)
  on.exit(setwd(old), add = TRUE)
  loaded <- expect_no_warning(LoadSeuratRds(file = rds))
  expect_identical(Layers(loaded), "counts")
  expect_equal(
    suppressWarnings(as.matrix(LayerData(loaded, layer = "counts"))),
    built$counts
  )
})

test_that("LoadSeuratRds finds absolute layer paths after a move", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("fs")
  built <- make_ondisk_object()
  rds <- save_and_move(object = built$object, relative = FALSE)
  loaded <- expect_no_warning(LoadSeuratRds(file = rds))
  expect_identical(Layers(loaded), "counts")
  expect_equal(
    suppressWarnings(as.matrix(LayerData(loaded, layer = "counts"))),
    built$counts
  )
})

test_that("LoadSeuratRds still loads an object that has not moved", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("fs")
  built <- make_ondisk_object()
  saved <- tempfile("saved")
  dir.create(saved)
  rds <- file.path(saved, "object.Rds")
  suppressWarnings(SaveSeuratRds(object = built$object, file = rds, move = TRUE))
  loaded <- expect_no_warning(LoadSeuratRds(file = rds))
  expect_identical(Layers(loaded), "counts")
  expect_equal(
    suppressWarnings(as.matrix(LayerData(loaded, layer = "counts"))),
    built$counts
  )
})
