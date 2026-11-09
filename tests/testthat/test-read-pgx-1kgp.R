# Tests against the 1000 Genomes high-coverage phased panel in
# data-raw/sample_vcf. These are skipped on CRAN/CI and when the files are not
# present locally.

test_that("every 1kGP gene can be read with readPGx", {
  genes <- kgp_genes()
  for (i in seq_len(nrow(genes))) {
    gene <- genes$gene[i]
    p <- suppressMessages(readPGx(genes$file[i], gene = gene))
    expect_true(validObject(p), info = gene)
    expect_gt(nrow(p), 0)
    expect_equal(ncol(p), 3202, info = gene)
    expect_equal(anyDuplicated(rownames(p)), 0, info = gene)

    rr <- SummarizedExperiment::rowRanges(p)
    key <- paste(as.character(rr), vapply(VariantAnnotation::alt(p), function(a) paste(as.character(a), collapse = ","), ""))
    expect_equal(anyDuplicated(key), 0, info = gene)
    expect_true(all(IRanges::overlapsAny(rr, availableGeneRanges(gene), type = "equal")), info = gene)
  }
})

test_that("NUDT15*5 carriers in 1kGP are called", {
  f <- kgp_file("chr13")
  skip_if(is.na(f), "chr13 1kGP file not available")
  p <- suppressMessages(readPGx(f, "NUDT15"))
  gt <- VariantAnnotation::geno(p)$GT
  row <- which(GenomicRanges::start(SummarizedExperiment::rowRanges(p)) == 48037798 &
               as.character(unlist(VariantAnnotation::alt(p))) == "A")
  skip_if(length(row) != 1, "NUDT15*5 position not present")
  carriers <- colnames(p)[grepl("1", gt[row, ])]
  skip_if(length(carriers) == 0, "No NUDT15*5 carriers in data")

  calls <- calls_of(run_pipeline(f, "NUDT15"))
  expect_true(any(grepl("\\*5", calls[carriers])))
  expect_false(any(grepl("Amb", calls[carriers])))
})
