test_that("exactly the reference positions are returned for every gene and build", {
  for (b in c("GRCh38", "GRCh37")) {
    for (g in setdiff(availableGenes(b), "CYP2D6")) {
      ref <- availableGeneRanges(g, build = b)
      lines <- gene_lines(g, b, gt = "0|1", offtarget = TRUE)
      # Some off-target records overlap PGx ranges and reach the exact-match filter
      expect_gt(length(lines), length(ref) + 2)
      f <- write_test_vcf(lines)
      p <- suppressMessages(readPGx(f, g, build = b))
      rr <- SummarizedExperiment::rowRanges(p)
      label <- paste(g, b)
      expect_equal(nrow(p), length(ref), info = label)
      expect_true(all(IRanges::overlapsAny(rr, ref, type = "equal")), info = label)
      expect_true(all(IRanges::overlapsAny(ref, rr, type = "equal")), info = label)
      expect_false(any(duplicated(GenomicRanges::granges(rr))), info = label)
      expect_equal(anyDuplicated(rownames(p)), 0, info = label)
      expect_false(is.unsorted(rr), info = label)
    }
  }
})

test_that("overlapping reference ranges do not duplicate records", {
  # These genes have PharmVar ranges that overlap one another (e.g. an SNV
  # inside a deletion), which used to make tabix return records twice.
  for (b in c("GRCh38", "GRCh37")) {
    for (g in c("CYP2A6", "CYP2C9", "DPYD", "NUDT15")) {
      ref <- availableGeneRanges(g, build = b)
      expect_true(any(IRanges::countOverlaps(ref, ref) > 1), info = paste(g, b))
      f <- write_test_vcf(gene_lines(g, b, gt = "0|1"))
      p <- suppressMessages(readPGx(f, g, build = b))
      expect_equal(nrow(p), length(ref), info = paste(g, b))
      expect_equal(anyDuplicated(rownames(p)), 0, info = paste(g, b))
    }
  }
})

test_that("a record inside the CYP2C9*25 deletion range is returned once", {
  lines <- c(
    vrec("chr10", 94942212, "AAGAAATGGAA", "A", "0|1"),
    vrec("chr10", 94942216, "A", "G", "1|0"),
    vrec("chr10", 94942290, "C", "T", "0|0")
  )
  p <- readPGx(write_test_vcf(lines, contigs = std_contigs()), "CYP2C9")
  expect_equal(nrow(p), 3)
  expect_false(any(duplicated(GenomicRanges::granges(SummarizedExperiment::rowRanges(p)))))
})

test_that("off-target, partially overlapping and symbolic records are excluded", {
  lines <- c(
    vrec("chr10", 94942205, "AAAAAAAAAAAA", "A", "0|1"),           # deletion spanning *25 start, not equal
    vrec("chr10", 94942212, "A", "G", "0|1"),                      # SNV at the start of the *25 deletion range
    vrec("chr10", 94942280, "N", "<DEL>", "0|1", info = "END=94942300;SVTYPE=DEL"), # symbolic SV spanning *2
    vrec("chr10", 94942289, "G", "A", "0|1"),                      # adjacent off-target SNV
    vrec("chr10", 94942290, "C", "T", "0|1"),                      # *2, kept
    vrec("chr10", 95000000, "G", "A", "0|1")                       # far off-target
  )
  f <- write_test_vcf(lines, contigs = std_contigs())
  p <- suppressMessages(readPGx(f, "CYP2C9"))
  expect_identical(rownames(p), "chr10:94942290_C/T")
})

test_that("non-normalized indels overlapping PGx positions are reported", {
  # CYP2C9*25 is AAGAAATGGAA>A at 94942212; a right-shifted representation
  # of a deletion does not match and is dropped with a message.
  lines <- c(
    vrec("chr10", 94942213, "AGAAATGGAAG", "A", "0|1"),
    vrec("chr10", 94942290, "C", "T", "0|1")
  )
  f <- write_test_vcf(lines, contigs = std_contigs())
  expect_message(p <- readPGx(f, "CYP2C9"), "normalized")
  expect_equal(nrow(p), 1)
})

test_that("a VCF without the gene's chromosome gives an informative error", {
  f <- write_test_vcf(vrec("chr22", 100000, "A", "G", "0|1"), contigs = std_contigs())
  err <- expect_error(readPGx(f, "CYP2C9"), "chr10")
  expect_match(conditionMessage(err), "CYP2C9")
  expect_match(conditionMessage(err), "chr22")
  expect_no_match(conditionMessage(err), "scanTabix")
})

test_that("a header-only VCF gives an informative error", {
  f <- write_test_vcf(character(), contigs = std_contigs())
  err <- expect_error(readPGx(f, "CYP2C9"), "CYP2C9")
  expect_match(conditionMessage(err), "none")
})

test_that("zero matching positions on the right chromosome warns and returns a valid object", {
  f <- write_test_vcf(vrec("chr10", 95000000, "G", "A", "0|1"), contigs = std_contigs())
  expect_warning(p <- readPGx(f, "CYP2C9"), "No overlapping ranges")
  expect_equal(nrow(p), 0)
  expect_true(validObject(p))
  expect_error(determineCallableAlleles(p), "no callable alleles")
})

test_that("FILTER values are not used to exclude records", {
  lines <- vrec("chr10", 94942290, "C", "T", "0|1", filter = "LowQual")
  f <- write_test_vcf(lines, contigs = std_contigs(), extra_header = '##FILTER=<ID=LowQual,Description="Low quality">')
  p <- readPGx(f, "CYP2C9")
  expect_equal(nrow(p), 1)
  expect_identical(as.character(VariantAnnotation::filt(p)), "LowQual")
})

test_that("rows are sorted and rownames are canonical and unique", {
  ref <- availableGeneRanges("CYP2C9")
  lines <- gene_lines("CYP2C9", chrom = "10", gt = "0|1")
  f <- write_test_vcf(lines, contigs = "10")
  p <- suppressMessages(readPGx(f, "CYP2C9"))
  expect_false(is.unsorted(SummarizedExperiment::rowRanges(p)))
  expect_true(all(startsWith(rownames(p), "chr10:")))

  # Records sharing an rsID (split multiallelic sites) get unique names
  lines <- c(
    vrec("chr10", 94942249, "C", "G", "0|1", id = "rs1"),
    vrec("chr10", 94942249, "C", "T", "0|0", id = "rs1"),
    vrec("chr10", 94942290, "C", "T", "0|1", id = "rs1799853")
  )
  p <- readPGx(write_test_vcf(lines, contigs = std_contigs()), "CYP2C9")
  expect_equal(nrow(p), 3)
  expect_equal(anyDuplicated(rownames(p)), 0)
  expect_true("rs1799853" %in% rownames(p))
})

test_that("records without an ID are named CHROM:POS_REF/ALT", {
  p <- readPGx(write_test_vcf(vrec("chr10", 94942290, "C", "T", "0|1"), contigs = std_contigs()), "CYP2C9")
  expect_identical(rownames(p), "chr10:94942290_C/T")
})

test_that("split multiallelic records are both kept", {
  lines <- c(
    vrec("chr10", 94942249, "C", "G", "0|1"),
    vrec("chr10", 94942249, "C", "T", "1|0")
  )
  p <- readPGx(write_test_vcf(lines, contigs = std_contigs()), "CYP2C9")
  expect_equal(nrow(p), 2)
  expect_setequal(as.character(unlist(VariantAnnotation::alt(p))), c("G", "T"))
})

test_that("genuinely duplicated input records are kept and reported", {
  lines <- c(
    vrec("chr10", 94942290, "C", "T", "0|1"),
    vrec("chr10", 94942290, "C", "T", "1|0")
  )
  f <- write_test_vcf(lines, contigs = std_contigs())
  expect_warning(p <- readPGx(f, "CYP2C9"), "duplicated record")
  expect_equal(nrow(p), 2)
})

test_that("a multi-chromosome VCF returns only the gene's chromosome", {
  lines <- c(
    gene_lines("CYP2C9", gt = "0|1"),
    gene_lines("CYP4F2", gt = "0|1")
  )
  f <- write_test_vcf(lines, contigs = std_contigs())
  p <- suppressMessages(readPGx(f, "CYP2C9"))
  expect_identical(GenomeInfoDb::seqlevelsInUse(p), "chr10")
  expect_equal(nrow(p), length(availableGeneRanges("CYP2C9")))
  p2 <- suppressMessages(readPGx(f, "CYP4F2"))
  expect_identical(GenomeInfoDb::seqlevelsInUse(p2), "chr19")
})

test_that("a record spanning several query blocks is returned once", {
  # INFO/END makes tabix treat this record as spanning later PGx positions,
  # which belong to separate reduced query blocks.
  ref <- availableGeneRanges("CYP2C9")
  st <- sort(GenomicRanges::start(ref))
  lines <- c(
    vrec("chr10", st[1], "N", "A", "0|1", info = paste0("END=", st[4])),
    gene_lines("CYP2C9", gt = "0|1")[-1]
  )
  f <- write_test_vcf(lines, contigs = std_contigs())
  expect_no_warning(p <- suppressMessages(readPGx(f, "CYP2C9")))
  expect_equal(nrow(p), length(ref))
  expect_equal(sum(GenomicRanges::start(SummarizedExperiment::rowRanges(p)) == st[1]), 1)
  expect_false("paramRangeID" %in% names(S4Vectors::mcols(SummarizedExperiment::rowRanges(p))))

  # A long deletion covering several PGx positions is reported once
  del <- vrec("chr10", st[1] - 2, strrep("N", st[4] - st[1] + 4), "N", "0|1")
  f2 <- write_test_vcf(c(del, gene_lines("CYP2C9", gt = "0|1")), contigs = std_contigs())
  expect_message(p2 <- readPGx(f2, "CYP2C9"), "^1 record\\(s\\) overlap")
  expect_equal(nrow(p2), length(ref))
})

test_that("records dropped by the exact-match filter do not trigger genotype warnings", {
  # An unphased SNV inside the CYP2C9*25 deletion span, but not a PGx position
  lines <- c(
    vrec("chr10", 94942212, "AAGAAATGGAA", "A", "0|1"),
    vrec("chr10", 94942215, "C", "T", "0/1"),
    vrec("chr10", 94942290, "C", "T", "0|1")
  )
  f <- write_test_vcf(lines, contigs = std_contigs())
  expect_no_warning(expect_no_message(p <- readPGx(f, "CYP2C9")))
  expect_equal(nrow(p), 2)
})
