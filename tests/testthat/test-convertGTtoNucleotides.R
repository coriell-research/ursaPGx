test_that("convertGTtoNucleotides is able to construct a valid matrix for each sample", {
  df <- kgp_genes()

  # Read in all gene data and test if each gene has callable alleles
  for (i in seq_len(nrow(df))) {
    file <- df[i, "file", drop = TRUE]
    gene <- df[i, "gene", drop = TRUE]
    p <- suppressMessages(readPGx(file, gene = gene))

    p <- determineCallableAlleles(p)
    p <- buildReferenceDataFrame(p)
    orig <- geno(p)$GT

    p <- convertGTtoNucleotides(p)
    converted <- geno(p)$GT

    expect_equal(all(dim(orig) == dim(converted)), TRUE, label = "Converted matrix has same dimensions")
    expect_equal(all(grepl("[ACGT]", converted)), TRUE, label = "All positions have been converted to nucleotides")
  }
})
