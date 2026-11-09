# Each scenario reads CYP2C19 (single chromosome, no overlapping ranges) and
# checks that every reference position is recovered with canonical seqinfo.
expect_canonical <- function(p, gene, build) {
  ref <- availableGeneRanges(gene, build = build)
  expect_s4_class(p, "PGx")
  expect_equal(nrow(p), length(ref))
  expect_identical(
    GenomeInfoDb::seqlevelsInUse(p),
    GenomeInfoDb::seqlevelsInUse(ref)
  )
  expect_identical(GenomeInfoDb::seqinfo(p), GenomeInfoDb::seqinfo(ref))
  expect_true(all(IRanges::overlapsAny(SummarizedExperiment::rowRanges(p), ref, type = "equal")))
  # Downstream overlap operations must not fail on incompatible genomes
  expect_no_error(IRanges::subsetByOverlaps(SummarizedExperiment::rowRanges(p), ref, type = "equal"))
}

read_quiet <- function(f, gene = "CYP2C19", build = "GRCh38") {
  suppressMessages(readPGx(f, gene, build = build))
}

test_that("package reference ranges use UCSC (GRCh38) and NCBI (GRCh37) styles", {
  for (g in setdiff(availableGenes(), "CYP2D6")) {
    expect_identical(GenomeInfoDb::seqlevelsStyle(availableGeneRanges(g, "GRCh38")), "UCSC", info = g)
    expect_identical(GenomeInfoDb::seqlevelsStyle(availableGeneRanges(g, "GRCh37")), "NCBI", info = g)
  }
})

test_that("GRCh38 VCFs with UCSC names are read for any header", {
  analysis_set <- c(std_contigs("GRCh38", "UCSC"), "chr1_KI270706v1_random", "chrUn_KI270302v1", "chrEBV", "HLA-A*01:01:01:01")
  headers <- list(full = std_contigs("GRCh38", "UCSC"), target = "chr10", analysis = analysis_set, none = character())
  for (h in names(headers)) {
    f <- write_test_vcf(gene_lines("CYP2C19", chrom = "chr10", gt = "0|1"), contigs = headers[[h]])
    expect_no_warning(p <- read_quiet(f))
    expect_canonical(p, "CYP2C19", "GRCh38")
  }
})

test_that("GRCh38 VCFs with NCBI names are converted to UCSC", {
  headers <- list(
    full = std_contigs("GRCh38", "NCBI"),
    target = "10",
    none = character(),
    numeric_with_chrM = c(as.character(1:22), "X", "Y", "chrM"),
    unplaced = c(std_contigs("GRCh38", "NCBI"), "KI270706.1", "EBV")
  )
  for (h in names(headers)) {
    f <- write_test_vcf(gene_lines("CYP2C19", chrom = "10", gt = "0|1"), contigs = headers[[h]])
    expect_message(p <- readPGx(f, "CYP2C19"), "10 -> chr10|chr10 -> 10")
    expect_canonical(p, "CYP2C19", "GRCh38")
  }
})

test_that("GRCh37 VCFs with NCBI names are read for any header", {
  headers <- list(
    b37 = c(std_contigs("GRCh37", "NCBI"), "GL000192.1"),
    hs37d5 = c(std_contigs("GRCh37", "NCBI"), "NC_007605", "hs37d5"),
    none = character(),
    numeric_with_chrM = c(as.character(1:22), "X", "Y", "chrM")
  )
  for (h in names(headers)) {
    f <- write_test_vcf(gene_lines("CYP2C19", "GRCh37", chrom = "10", gt = "0|1"), contigs = headers[[h]])
    expect_no_warning(p <- read_quiet(f, build = "GRCh37"))
    expect_canonical(p, "CYP2C19", "GRCh37")
  }
})

test_that("hg19-style GRCh37 VCFs keep the GRCh37 genome", {
  headers <- list(
    full = paste0("chr", c(1:22, "X", "Y", "M")),
    no_chrM = paste0("chr", c(1:22, "X", "Y")),
    target = "chr10",
    randoms = c(paste0("chr", c(1:22, "X", "Y", "M")), "chr1_gl000191_random", "chrUn_gl000211"),
    none = character()
  )
  for (h in names(headers)) {
    f <- write_test_vcf(gene_lines("CYP2C19", "GRCh37", chrom = "chr10", gt = "0|1"), contigs = headers[[h]])
    expect_no_warning(p <- read_quiet(f, build = "GRCh37"))
    expect_canonical(p, "CYP2C19", "GRCh37")
    expect_identical(unique(GenomeInfoDb::genome(p)), "GRCh37")
  }
})

test_that("chromosome names are detected from the records, not a stale header", {
  for (b in c("GRCh38", "GRCh37")) {
    # Header says UCSC, records say NCBI
    f1 <- write_test_vcf(gene_lines("CYP2C19", b, chrom = "10", gt = "0|1"), contigs = paste0("chr", c(1:22, "X", "Y", "M")))
    expect_no_warning(p1 <- read_quiet(f1, build = b))
    expect_canonical(p1, "CYP2C19", b)
    # Header says NCBI, records say UCSC
    f2 <- write_test_vcf(gene_lines("CYP2C19", b, chrom = "chr10", gt = "0|1"), contigs = c(1:22, "X", "Y", "MT"))
    expect_no_warning(p2 <- read_quiet(f2, build = b))
    expect_canonical(p2, "CYP2C19", b)
  }
})

test_that("a header mixing UCSC and NCBI names resolves to the names in use", {
  contigs <- c(std_contigs("GRCh38", "UCSC"), std_contigs("GRCh38", "NCBI"))
  f <- write_test_vcf(gene_lines("CYP2C19", chrom = "chr10", gt = "0|1"), contigs = contigs)
  expect_canonical(read_quiet(f), "CYP2C19", "GRCh38")
})

test_that("VCFs without contig header lines work for every build and style", {
  combos <- expand.grid(build = c("GRCh38", "GRCh37"), chrom = c("chr10", "10"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(combos))) {
    b <- combos$build[i]
    f <- write_test_vcf(gene_lines("CYP2C19", b, chrom = combos$chrom[i], gt = "0|1"))
    expect_canonical(read_quiet(f, build = b), "CYP2C19", b)
  }
})

test_that("RefSeq, dbSNP and case-variant chromosome names are supported", {
  for (chrom in c("NC_000010.11", "ch10", "Chr10", "CHR10")) {
    f <- write_test_vcf(gene_lines("CYP2C19", chrom = chrom, gt = "0|1"), contigs = chrom)
    p <- read_quiet(f)
    expect_canonical(p, "CYP2C19", "GRCh38")
    expect_true(all(startsWith(rownames(p), "chr10:")), info = chrom)
  }
  f <- write_test_vcf(gene_lines("CYP2C19", "GRCh37", chrom = "NC_000010.10", gt = "0|1"))
  expect_canonical(read_quiet(f, build = "GRCh37"), "CYP2C19", "GRCh37")
})

test_that("output is identical whatever the input naming", {
  gts <- c("0|1", "1|1", "0|0")
  names_in <- c("chr10", "10", "NC_000010.11", "ch10")
  ps <- lapply(names_in, function(chrom) {
    read_quiet(write_test_vcf(gene_lines("CYP2C19", chrom = chrom, gt = gts), samples = c("A", "B", "C")))
  })
  for (i in 2:length(ps)) {
    expect_identical(SummarizedExperiment::rowRanges(ps[[i]]), SummarizedExperiment::rowRanges(ps[[1]]))
    expect_identical(rownames(ps[[i]]), rownames(ps[[1]]))
    expect_identical(VariantAnnotation::geno(ps[[i]])$GT, VariantAnnotation::geno(ps[[1]])$GT)
  }
})

test_that("a VCF from the wrong build is diagnosed from header contig lengths", {
  # GRCh37 chr10 length declared in a file read as GRCh38
  len37 <- GenomeInfoDb::seqlengths(availableGeneRanges("CYP2C19", "GRCh37"))[["10"]]
  f <- write_test_vcf(
    gene_lines("CYP2C19", "GRCh37", chrom = "chr10", gt = "0|1"),
    contigs = "chr10", contig_lengths = c(chr10 = len37)
  )
  expect_error(read_quiet(f, build = "GRCh38"), "GRCh37")

  # Matching lengths are fine
  len38 <- GenomeInfoDb::seqlengths(availableGeneRanges("CYP2C19", "GRCh38"))[["chr10"]]
  f2 <- write_test_vcf(
    gene_lines("CYP2C19", chrom = "chr10", gt = "0|1"),
    contigs = "chr10", contig_lengths = c(chr10 = len38)
  )
  expect_canonical(read_quiet(f2), "CYP2C19", "GRCh38")
})

test_that("an empty result still has canonical seqinfo", {
  ref <- availableGeneRanges("CYP2C19")
  # A record on the right chromosome, away from every PGx position
  f <- write_test_vcf(vrec("10", max(GenomicRanges::end(ref)) + 10000, "A", "G", "0|1"), contigs = "10")
  expect_warning(p <- read_quiet(f), "No overlapping ranges")
  expect_equal(nrow(p), 0)
  expect_identical(GenomeInfoDb::seqinfo(p), GenomeInfoDb::seqinfo(ref))
})
