test_that("Phased allele caller functions", {
  vcf <- here::here(
    "data-raw",
    "sample_vcf",
    "1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    )
  p <- ursaPGx::readPGx(vcf, gene = "CYP2C19")
  result <- callPhasedDiplotypes(p)
  result_all <- callPhasedDiplotypes(p, exact=FALSE)
  result_sub <- callPhasedDiplotypes(p, summarize=FALSE)
  result_allSub <- callPhasedDiplotypes(p, exact=FALSE, summarize=FALSE)
  expect_equal(nrow(result), 3202)
  expect_equal(nrow(result_all), 3202)
  expect_equal(nrow(result_sub), 3202)
  expect_equal(nrow(result_allSub), 3202)
})
