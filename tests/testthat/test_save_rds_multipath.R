set.seed(42)

# A joined BPCells layer references its backing store once per sub-matrix, and
# .FilePath() encodes those as one comma-separated string. Saving has to
# understand that encoding, the way LoadSeuratRds() already does.
build_joined_bpcells <- function(store, nfeat = 120, ncell = 80) {
  counts <- matrix(rpois(nfeat * ncell, lambda = 2), nrow = nfeat, ncol = ncell)
  dimnames(counts) <- list(
    paste0("gene", seq_len(nfeat)),
    paste0("cell", seq_len(ncell))
  )
  suppressWarnings(BPCells::write_matrix_dir(mat = as.sparse(counts), dir = store))
  object <- suppressWarnings(
    CreateSeuratObject(counts = BPCells::open_matrix_dir(dir = store))
  )
  object <- split(object, f = rep(c("a", "b"), length.out = ncol(object)))
  object[["RNA"]] <- JoinLayers(object[["RNA"]])
  list(object = object, counts = counts)
}

test_that("a joined on-disk layer records several paths in one entry", {
  skip_if_not_installed("BPCells")
  built <- build_joined_bpcells(store = tempfile("bpcells-store"))
  path <- SeuratObject:::.FilePath(LayerData(built$object, layer = "counts"))
  expect_length(path, 1L)
  expect_true(grepl(",", path, fixed = TRUE))
})

test_that("SaveSeuratRds moves every path of a joined on-disk layer", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("fs")
  store <- tempfile("bpcells-store")
  built <- build_joined_bpcells(store = store)
  saved <- tempfile("saved")
  dir.create(saved)
  rds <- file.path(saved, "object.Rds")
  # Previously: Error, Can't find path: 'store,store'
  expect_no_error(
    suppressWarnings(SaveSeuratRds(object = built$object, file = rds, move = TRUE))
  )
  expect_true(file.exists(rds))

  loaded <- suppressWarnings(LoadSeuratRds(file = rds))
  expect_identical(Layers(loaded), "counts")
  expect_equal(
    suppressWarnings(as.matrix(LayerData(loaded, layer = "counts"))),
    built$counts
  )
})

test_that("relative = TRUE handles a multi-path entry", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("fs")
  built <- build_joined_bpcells(store = tempfile("bpcells-store"))
  saved <- tempfile("saved")
  dir.create(saved)
  rds <- file.path(saved, "object.Rds")
  expect_no_error(suppressWarnings(
    SaveSeuratRds(object = built$object, file = rds, move = TRUE, relative = TRUE)
  ))
  cache <- Tool(object = readRDS(rds), slot = "SaveSeuratRds")
  # every component relative, and none mangled into a single comma-bearing name
  parts <- unlist(strsplit(cache$path, split = ","))
  expect_gt(length(parts), 1L)
  expect_false(any(fs::is_absolute_path(parts)))
})

test_that("the path helpers handle repeats, order and single entries", {
  skip_if_not_installed("fs")
  dir <- tempfile("multi")
  dir.create(dir)
  a <- file.path(dir, "part-a"); dir.create(a)
  b <- file.path(dir, "part-b"); dir.create(b)
  dest <- tempfile("dest"); dir.create(dest)

  # a repeated path is moved once but still reported for each position
  moved <- SeuratObject:::.FileMoveAll(paste(a, b, a, sep = ","), new_path = dest)
  parts <- unlist(strsplit(moved, split = ","))
  expect_length(parts, 3L)
  expect_identical(parts[1], parts[3])
  expect_true(all(fs::dir_exists(parts)))

  # single paths are unchanged in shape
  single <- SeuratObject:::.FileMoveAll(parts[2], new_path = dest)
  expect_false(grepl(",", single, fixed = TRUE))

  rel <- SeuratObject:::.PathRelAll(paste(a, b, sep = ","), start = dir)
  expect_identical(rel, "part-a,part-b")
  expect_identical(SeuratObject:::.PathRelAll(a, start = dir), "part-a")
})
