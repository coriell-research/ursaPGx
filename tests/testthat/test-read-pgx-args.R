make_small_vcf <- function(...) {
  write_test_vcf(gene_lines("CYP2C9", gt = "0|1"), contigs = std_contigs(), ...)
}

test_that("build argument defaults to GRCh38 and rejects invalid values", {
  f <- make_small_vcf()
  expect_identical(pgxBuild(readPGx(f, "CYP2C9")), "GRCh38")
  expect_identical(pgxBuild(readPGx(f, "CYP2C9", build = c("GRCh38", "GRCh37"))), "GRCh38")

  for (b in list("T2T", "grch38", "hg38", NA, NA_character_, character(0),
                 c("GRCh37", "GRCh38"), c("GRCh38", "GRCh38"))) {
    expect_error(readPGx(f, "CYP2C9", build = b), "Genome must be one of", info = paste(b, collapse = ","))
  }
})

test_that("gene argument is validated", {
  f <- make_small_vcf()
  expect_error(readPGx(f, gene = "PTEN"), "Gene not in gene list")
  expect_error(readPGx(f, gene = "cyp2c9"), "availableGenes")
  expect_error(readPGx(f, gene = NA), "Exactly one gene")
  expect_error(readPGx(f, gene = character(0)), "Exactly one gene")
  expect_error(readPGx(f, gene = c("CYP2C9", "DPYD")), "Exactly one gene")
  expect_error(readPGx(f, gene = "CYP2D6"), "cyrius")
  expect_error(readPGx(f), "gene")
})

test_that("file argument is validated", {
  f <- make_small_vcf()
  expect_error(readPGx(c(f, f), "CYP2C9"), "Only a single file")
  expect_error(readPGx(character(0), "CYP2C9"), "Only a single file")
  expect_error(readPGx(NA_character_, "CYP2C9"), "single path")
  expect_error(readPGx(NA, "CYP2C9"), "single path")
  expect_error(readPGx("", "CYP2C9"), "single path")
  expect_error(readPGx(list(f), "CYP2C9"), "single path")
  expect_error(readPGx("/no/such/file.vcf.gz", "CYP2C9"), "not found")
  expect_error(readPGx("sample.vcf", "CYP2C9"), "indexed VCF")
  expect_error(readPGx("sample.VCF", "CYP2C9"), "indexed VCF")
  expect_error(readPGx("sample.bcf", "CYP2C9"), "BCF")
})

test_that("a named file vector is accepted", {
  f <- make_small_vcf()
  p <- readPGx(c(chr10 = f), "CYP2C9")
  expect_s4_class(p, "PGx")
})

test_that("missing index gives an informative error", {
  f <- make_small_vcf(index = "none")
  expect_error(readPGx(f, "CYP2C9"), "No index found")
})

test_that("a CSI-only index is accepted", {
  skip_if_not(has_bcftools(), "bcftools not available")
  f <- make_small_vcf(index = "csi")
  expect_false(file.exists(paste0(f, ".tbi")))
  p <- readPGx(f, "CYP2C9")
  expect_s4_class(p, "PGx")
  expect_equal(nrow(p), length(availableGeneRanges("CYP2C9")))
})

test_that("a non-BGZF file with a stale index errors instead of crashing", {
  good <- make_small_vcf()
  dir <- tempfile("ursaPGx-gzip-")
  dir.create(dir)
  bad <- file.path(dir, "plain.vcf.gz")
  con <- gzfile(bad, "w")
  writeLines(readLines(good), con)
  close(con)
  file.copy(paste0(good, ".tbi"), paste0(bad, ".tbi"))
  expect_error(readPGx(bad, "CYP2C9"), "bgzip")

  # An uncompressed file with a .gz name
  plain <- file.path(dir, "text.vcf.gz")
  writeLines(readLines(good), plain)
  file.copy(paste0(good, ".tbi"), paste0(plain, ".tbi"))
  expect_error(readPGx(plain, "CYP2C9"), "bgzip")
})

test_that("a sites-only VCF is rejected", {
  ref <- availableGeneRanges("CYP2C9")
  lines <- vrec("chr10", GenomicRanges::start(ref)[1], strrep("N", GenomicRanges::width(ref)[1]), "A")
  f <- write_test_vcf(lines, contigs = std_contigs(), sites_only = TRUE)
  expect_error(readPGx(f, "CYP2C9"), "no samples")
})
