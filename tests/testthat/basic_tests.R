library(testthat)
library(ReAnnotateR)

#using existing test data from vignettes
load(system.file("data/example_bed.rda", package = "ReAnnotateR"))
load(system.file("data/example_annotation.rda", package = "ReAnnotateR"))
load(system.file("data/example_enrichment.rda", package = "ReAnnotateR"))

gr <- ReAnnotateR::read_bed(example_bed)

test_that("read_bed works", {
  expect_error(ReAnnotateR::read_bed(example_bed), NA)
})

test_that("annotate_regions works", {
  expect_error({
    gr_annotated <- ReAnnotateR::annotate_regions(gr, txdb_info = example_annotation)
  }, NA)
})

test_that("nearest_gene works", {
  gr_annotated <- ReAnnotateR::annotate_regions(gr, txdb_info = example_annotation)
  expect_error({
    gr_nearest <- ReAnnotateR::nearest_gene(gr_annotated, txdb_info = example_annotation)
  }, NA)
})

test_that("functional_terms works", {
  gr_annotated <- ReAnnotateR::annotate_regions(gr, txdb_info = example_annotation)
  gr_nearest <- ReAnnotateR::nearest_gene(gr_annotated, txdb_info = example_annotation)
  expect_error({
    enrichment <- ReAnnotateR::functional_terms(gr_nearest, OrgDb = example_enrichment)
  }, NA)
})

test_that("plot functions run", {
  gr_annotated <- ReAnnotateR::annotate_regions(gr, txdb_info = example_annotation)
  gr_nearest <- ReAnnotateR::nearest_gene(gr_annotated, txdb_info = example_annotation)
  enrichment <- ReAnnotateR::functional_terms(gr_nearest, OrgDb = example_enrichment)
  
  expect_error(ReAnnotateR::plot_feature_composition(gr_annotated), NA)
  expect_error(ReAnnotateR::plot_chromosomal_density(gr_annotated), NA)
  expect_error(ReAnnotateR::plot_enrichment(enrichment), NA)
  expect_error(ReAnnotateR::plot_ideogram(gr_annotated), NA)
})

test_that("fisher_enrichment works", {
  gr_annotated <- ReAnnotateR::annotate_regions(gr, txdb_info = example_annotation)
  expect_error({
    fisher_results <- ReAnnotateR::fisher_enrichment(gr_annotated, category_col = "feature")
  }, NA)
})

test_that("export_results works", {
  gr_annotated <- ReAnnotateR::annotate_regions(gr, txdb_info = example_annotation)
  gr_nearest <- ReAnnotateR::nearest_gene(gr_annotated, txdb_info = example_annotation)
  enrichment <- ReAnnotateR::functional_terms(gr_nearest, OrgDb = example_enrichment)
  
  expect_error({
    ReAnnotateR::export_results(
      annotated_gr = gr_nearest,
      enrichment_results = enrichment,
      output_dir = tempdir()
    )
  }, NA)
})
