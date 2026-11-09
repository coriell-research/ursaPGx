test_that("callableAlleles can be extracted for every sample", {
  df <- kgp_genes()

    # Read in all gene data and test if each gene has callable alleles
    for (i in seq_len(nrow(df))) {
        file <- df[i, "file", drop = TRUE]
        gene <- df[i, "gene", drop = TRUE]
        p <- suppressMessages(readPGx(file, gene = gene))
    

        p <- determineCallableAlleles(p)
        expect_type(pgxCallableAlleles(p), "character")
        expect_equal(length(pgxCallableAlleles(p)) > 0, TRUE, label = "There exists callable alleles")
    }
})
