test_that("buildReferenceDataFrame constructs a valid reference for each sample", {
  df <- kgp_genes()

  # Read in all gene data and test if each gene has callable alleles
  for (i in seq_len(nrow(df))) {
    file <- df[i, "file", drop = TRUE]
    gene <- df[i, "gene", drop = TRUE]
    p <- suppressMessages(readPGx(file, gene = gene))

    p <- determineCallableAlleles(p)
    p <- buildReferenceDataFrame(p)
    DF <- pgxReferenceDataFrame(p)

    expect_equal(is(DF, "DFrame"), TRUE, label = "Reference DF is a vaild DataFrame")
  }
})
