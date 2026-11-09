test_that("readPGx returns a well-formed PGx object", {
  f <- write_test_vcf(gene_lines("CYP2C9", gt = c("0|1", "1|0", "0|0")), contigs = std_contigs(),
                      samples = c("S1", "S2", "S3"))
  p <- readPGx(f, "CYP2C9")

  expect_s4_class(p, "PGx")
  expect_s4_class(p, "CollapsedVCF")
  expect_true(validObject(p))
  expect_identical(pgxGene(p), "CYP2C9")
  expect_identical(pgxBuild(p), "GRCh38")
  expect_identical(pgxCallableAlleles(p), "")
  expect_false(p@hasCallableAlleles)
  expect_false(p@hasReferenceDataFrame)
  expect_s4_class(pgxReferenceDataFrame(p), "DFrame")
  expect_equal(dim(pgxReferenceDataFrame(p)), c(0, 0))
  expect_equal(dim(pgxGenotypeMatrix(p)), dim(p))
  expect_identical(colnames(p), c("S1", "S2", "S3"))

  out <- capture.output(show(p))
  expect_true(any(grepl("PGx gene: CYP2C9", out)))
  expect_true(any(grepl("Number of samples: 3", out)))
})

test_that("GRCh37 objects are valid and accept GRCh37 callable alleles", {
  f <- write_test_vcf(gene_lines("CYP2C9", "GRCh37", gt = "0|1"), contigs = std_contigs("GRCh37", "NCBI"))
  p <- readPGx(f, "CYP2C9", build = "GRCh37")
  expect_identical(pgxBuild(p), "GRCh37")
  expect_true(validObject(p))
  p <- determineCallableAlleles(p)
  expect_gt(length(pgxCallableAlleles(p)), 0)
  expect_true(validObject(p))
})

test_that("single samples and unusual sample names are preserved", {
  f <- write_test_vcf(gene_lines("CYP2C9", gt = "0|1"), contigs = std_contigs(), samples = "HG00096")
  p <- readPGx(f, "CYP2C9")
  expect_identical(colnames(p), "HG00096")
  expect_equal(ncol(pgxGenotypeMatrix(p)), 1)

  nms <- c("NA-12878.b", "sample 3")
  f <- write_test_vcf(gene_lines("CYP2C9", gt = c("0|1", "1|1")), contigs = std_contigs(), samples = nms)
  expect_identical(colnames(readPGx(f, "CYP2C9")), nms)
})

test_that("unphased genotypes produce a warning", {
  f <- write_test_vcf(gene_lines("CYP2C9", gt = c("0/1", "0|0")), contigs = std_contigs(), samples = c("A", "B"))
  expect_warning(readPGx(f, "CYP2C9"), "unphased")

  f <- write_test_vcf(gene_lines("CYP2C9", gt = c("0|1", "0|0")), contigs = std_contigs(), samples = c("A", "B"))
  expect_no_warning(readPGx(f, "CYP2C9"))
})

test_that("missing and haploid genotypes produce a warning", {
  f <- write_test_vcf(gene_lines("CYP2C9", gt = c(".|.", "0|0")), contigs = std_contigs(), samples = c("A", "B"))
  expect_warning(readPGx(f, "CYP2C9"), "missing")

  f <- write_test_vcf(gene_lines("CYP2C9", gt = c("1", "0|0")), contigs = std_contigs(), samples = c("A", "B"))
  expect_warning(readPGx(f, "CYP2C9"), "haploid")
})
