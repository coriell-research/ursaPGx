# End-to-end diplotype calling on synthetic VCFs carrying real PharmVar
# alleles. These guard the parts of readPGx() that change what reaches the
# caller: record duplication, split multiallelic sites and seqinfo handling.

cyp2c9_carriers <- list(
  S2 = c("CYP2C9*2", NA),
  S3 = c(NA, "CYP2C9*3"),
  S25 = c("CYP2C9*25", NA),
  S41 = c("CYP2C9*41", NA),
  S26 = c("CYP2C9*26", NA),
  S44 = c("CYP2C9*44", NA),
  REF = c(NA, NA)
)
cyp2c9_expected <- c(
  S2 = "*2|*1", S3 = "*1|*3", S25 = "*25|*1", S41 = "*41|*1",
  S26 = "*26|*1", S44 = "*44|*1", REF = "*1|*1"
)

call_cyp2c9 <- function(..., build = "GRCh38", contigs = std_contigs()) {
  lines <- pv_lines("CYP2C9", cyp2c9_carriers, build = build, ...)
  f <- write_test_vcf(lines, contigs = contigs, samples = names(cyp2c9_carriers))
  calls_of(run_pipeline(f, "CYP2C9", build = build))
}

test_that("CYP2C9 alleles are called, including those at overlapping ranges", {
  # *41 (SNV) lies inside the *25 deletion range
  expect_identical(call_cyp2c9(), cyp2c9_expected)
})

test_that("split multiallelic records are called correctly in either order", {
  # *26 (C>G) and *44 (C>T) share position 94942249
  expect_identical(call_cyp2c9(split = TRUE, reverse_split = FALSE), cyp2c9_expected)
  expect_identical(call_cyp2c9(split = TRUE, reverse_split = TRUE), cyp2c9_expected)
})

test_that("merged multi-ALT records are called correctly", {
  expect_identical(call_cyp2c9(split = FALSE), cyp2c9_expected)
})

test_that("calls do not depend on the input chromosome naming", {
  for (chrom in c("10", "NC_000010.11", "ch10")) {
    expect_identical(call_cyp2c9(chrom = chrom, contigs = chrom), cyp2c9_expected, info = chrom)
  }
  # A numeric header that also lists chrM used to produce mixed genome tags
  expect_identical(
    call_cyp2c9(chrom = "10", contigs = c(as.character(1:22), "X", "Y", "chrM")),
    cyp2c9_expected
  )
})

test_that("an unrelated split record before a PGx allele does not break calling", {
  carriers <- list(S2 = c("CYP2C9*2", NA), REF = c(NA, NA))
  lines <- pv_lines("CYP2C9", carriers)
  pos <- grep("^chr10\t94942290\t", lines)
  decoys <- c(
    vrec("chr10", 94942290, "C", "A", c("0|0", "0|0")),  # non-PharmVar ALT
    vrec("chr10", 94942290, "C", "*", c("0|0", "0|0"))   # spanning deletion
  )
  lines <- append(lines, decoys, after = pos - 1)
  f <- write_test_vcf(lines, contigs = std_contigs(), samples = names(carriers))
  calls <- calls_of(run_pipeline(f, "CYP2C9"))
  expect_identical(calls, c(S2 = "*2|*1", REF = "*1|*1"))
})

test_that("NUDT15*5 is called inside the *9 deletion range", {
  carriers <- list(S5 = c("NUDT15*5", NA), S9 = c("NUDT15*9", NA), S3 = c("NUDT15*3", "NUDT15*3"))
  f <- write_test_vcf(pv_lines("NUDT15", carriers), contigs = std_contigs(), samples = names(carriers))
  expect_identical(calls_of(run_pipeline(f, "NUDT15")), c(S5 = "*5|*1", S9 = "*9|*1", S3 = "*3|*3"))
})

test_that("GRCh37 calling works for NCBI and hg19-style inputs", {
  expect_identical(call_cyp2c9(build = "GRCh37", contigs = std_contigs("GRCh37", "NCBI")), cyp2c9_expected)
  expect_identical(
    call_cyp2c9(build = "GRCh37", chrom = "chr10", contigs = paste0("chr", c(1:22, "X", "Y", "M"))),
    cyp2c9_expected
  )
})

test_that("exact duplicate records do not break calling", {
  carriers <- list(S2 = c("CYP2C9*2", NA), S3 = c(NA, "CYP2C9*3"), REF = c(NA, NA))
  lines <- pv_lines("CYP2C9", carriers)
  pos <- grep("^chr10\t94942290\t", lines)
  lines <- append(lines, lines[pos], after = pos)
  f <- write_test_vcf(lines, contigs = std_contigs(), samples = names(carriers))
  expect_warning(calls <- calls_of(run_pipeline(f, "CYP2C9")), "duplicated record")
  expect_identical(calls, c(S2 = "*2|*1", S3 = "*1|*3", REF = "*1|*1"))
})
