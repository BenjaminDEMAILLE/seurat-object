# Every list slot of a Seurat object is validated with pass.zero = TRUE, which
# allows the slot to be empty, except `assays`. That made an object with no
# assays invalid, so the class could neither be constructed by new() nor
# extended, since new() builds from the prototype.

test_that("a Seurat object can be constructed empty", {
  expect_no_error(new("Seurat"))
  expect_s4_class(new("Seurat"), "Seurat")
})

test_that("the class can be extended", {
  setClass(Class = "TestExtended", contains = "Seurat", slots = c(newslot = "list"))
  on.exit(removeClass("TestExtended"), add = TRUE)
  expect_no_error(new("TestExtended"))

  extended <- as(pbmc_small, "TestExtended")
  expect_s4_class(extended, "TestExtended")
  expect_true(is(extended, "Seurat"))
  # the data survives the coercion
  expect_identical(Assays(extended), Assays(pbmc_small))
  expect_equal(ncol(extended), ncol(pbmc_small))
  expect_identical(colnames(extended), colnames(pbmc_small))

  extended@newslot <- list(a = 1)
  expect_no_error(validObject(extended))
  expect_identical(names(extended@newslot), "a")
})

test_that("ordinary objects are still valid", {
  expect_no_error(validObject(pbmc_small))
})

test_that("an unnamed assay list is still rejected", {
  object <- new("Seurat")
  slot(object, "assays") <- list(pbmc_small[["RNA"]])
  expect_error(validObject(object), "named list")
})

test_that("assays with an empty name are still rejected", {
  object <- new("Seurat")
  slot(object, "assays") <- setNames(list(pbmc_small[["RNA"]]), "")
  expect_error(validObject(object), "named list")
})

test_that("IsNamedList treats an empty list as the callers expect", {
  expect_true(is.function(SeuratObject:::IsNamedList))
  expect_false(SeuratObject:::IsNamedList(list()))
  expect_true(SeuratObject:::IsNamedList(list(), pass.zero = TRUE))
  expect_true(SeuratObject:::IsNamedList(list(a = 1)))
})
