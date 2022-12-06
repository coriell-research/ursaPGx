test_that("Stop conditions catch errors on readPGx", {
  file <- here::here("data-raw", "sample_vcf", "1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
  expect_error(readPGx(file, gene = "PTEN"), "Gene not in gene list")
  expect_error(readPGx(file, gene = "CYP2D6", build = "T2T"), "Genome not available")
  expect_error(readPGx(file, gene = c("CYP2D6", "DPYD")), "Multple genes supplied to function")
})

test_that("Every 1kg sample can be read in with readPGx", {
  chroms <- c("chr10", "chr19", "chr22", "chr1", "chr11", "chr7", "chr13", "chr12")
  d <- here::here("data-raw", "sample_vcf")
  f <- paste0("1kGP_high_coverage_Illumina.", chroms, ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
  files <- file.path(d, f)
  names(files) <- chroms

  # Mapping for chromosome to genes
  df <- data.frame(
    chrom = c("chr19", "chr19", "chr19", "chr10", "chr10", "chr10", 
              "chr22", "chr7", "chr7", "chr19", "chr13", "chr12"
    ),
    gene = c("CYP2A13", "CYP2A6", "CYP2B6", "CYP2C19", "CYP2C8", "CYP2C9", 
             "CYP2D6", "CYP3A4", "CYP3A5", "CYP4F2", "NUDT15", "SLCO1B1"
    )
  )

  # Read in all gene data and test if valid objects can be created
  for (i in 1:nrow(df)) {
    chrom <- df[i, "chrom", drop = TRUE]
    gene <- df[i, "gene", drop = TRUE]
    p <- readPGx(files[chrom], gene = gene)
    expect_true(validObject(p))

    nranges <- length(SummarizedExperiment::rowRanges(p))
    expect_gt(nranges, 0)

    samples <- ncol(p)
    expect_equal(samples, 3202)
  }
})
