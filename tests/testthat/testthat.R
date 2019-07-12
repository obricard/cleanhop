library(testthat)
library(cleanhop)

sce <- SingleCellExperiment(assays = list(counts = as.matrix(cleanhop_counts)), colData=cleanhop_annotations)
cleaned_sce <- cleanhop(sce)

test_that("cleanhop test example data", {
  expect_equal(as.data.frame(counts(cleaned_sce)),as.data.frame(cleaned_cleanhop_counts))
})


#test_check("cleanhop")
