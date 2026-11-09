test_that("callPhasedDiplotypes produces a valid output for each gene", {
  df <- kgp_genes()

  # Read in all gene data and test if each gene has callable alleles
  for (i in seq_len(nrow(df))) {
    file <- df[i, "file", drop = TRUE]
    gene <- df[i, "gene", drop = TRUE]
    p <- suppressMessages(readPGx(file, gene = gene))

    p <- determineCallableAlleles(p)
    p <- buildReferenceDataFrame(p)
    p <- convertGTtoNucleotides(p)
    result <- callPhasedDiplotypes(p)

    expect_equal(is(result, "DFrame"), TRUE, label = "Output is a DataFrame")
    }
})
