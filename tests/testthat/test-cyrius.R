test_that("Cyrius Python requirements are declared on load", {
  pkgs <- reticulate::py_require()$packages
  expect_true(all(
    c("numpy>=1.16", "scipy>=1.2", "pysam>=0.15.3", "statsmodels>=0.9") %in% pkgs
  ))
})

test_that("cyrius() validates input files", {
  expect_error(cyrius("does-not-exist.bam"), "do not exist")

  tmp <- tempfile(fileext = ".txt")
  bam <- tempfile(fileext = ".BAM")
  file.create(c(tmp, bam))
  on.exit(unlink(c(tmp, bam)), add = TRUE)

  expect_error(cyrius(tmp), "\\.bam or \\.cram")
  expect_error(cyrius(bam), "not indexed")
})

test_that("Cyrius runs", {
  skip_on_cran()
  skip_if_not_installed("here")

  cram_dir <- here::here("data-raw", "sample_cram")
  fa <- here::here("data-raw", "fasta", "GRCh38_full_analysis_set_plus_decoy_hla.fa")
  skip_if_not(dir.exists(cram_dir), "sample CRAMs not available")
  skip_if_not(file.exists(fa), "reference fasta not available")

  cram <- list.files(path = cram_dir, pattern = "\\.cram$", full.names = TRUE)
  names(cram) <- gsub(".final", "", tools::file_path_sans_ext(basename(cram)))
  result <- cyrius(cram, fa)

  expect_s4_class(result, "DataFrame")
  expect_equal(colnames(result), c("Sample", "Genotype", "Filter"))
  expect_equal(result$Sample, names(cram))
})
