# Helpers for building small, synthetic, bgzipped + tabix-indexed VCF files.
#
# These let readPGx() tests run without the large 1000 Genomes files in
# data-raw/ and make it easy to control chromosome naming, header contigs,
# record content and genotypes independently.

# Build one VCF body line. `gt` holds one genotype per sample; NULL produces a
# sites-only record.
vrec <- function(chr, pos, ref, alt, gt = NULL, id = ".", filter = "PASS", info = ".") {
  fixed <- c(chr, format(pos, scientific = FALSE, trim = TRUE), id, ref, alt, ".", filter, info)
  if (is.null(gt)) {
    return(paste(fixed, collapse = "\t"))
  }
  paste(c(fixed, "GT", gt), collapse = "\t")
}

# Standard contig names for a build in a given naming style
std_contigs <- function(build = "GRCh38", style = c("UCSC", "NCBI", "dbSNP")) {
  style <- match.arg(style)
  chrs <- c(as.character(1:22), "X", "Y")
  switch(style,
    UCSC = paste0("chr", c(chrs, "M")),
    NCBI = c(chrs, "MT"),
    dbSNP = paste0("ch", c(chrs, "MT"))
  )
}

# Write VCF lines to a bgzipped, indexed file and return its path.
#   contigs        character vector of ##contig IDs written to the header
#   contig_lengths optional named integer vector of lengths for header contigs
#   samples        sample names (must match the number of GT columns)
#   index          "tbi", "csi" (requires bcftools) or "none"
#   sites_only     write a header without FORMAT/sample columns
write_test_vcf <- function(lines, contigs = character(), samples = "S1",
                           contig_lengths = NULL, index = c("tbi", "csi", "none"),
                           sites_only = FALSE, extra_header = character(),
                           name = "test") {
  index <- match.arg(index)
  dir <- tempfile("ursaPGx-vcf-")
  dir.create(dir)

  contig_lines <- if (length(contigs)) {
    len <- if (!is.null(contig_lengths)) contig_lengths[contigs] else rep(NA, length(contigs))
    ifelse(
      is.na(len),
      paste0("##contig=<ID=", contigs, ">"),
      paste0("##contig=<ID=", contigs, ",length=", len, ">")
    )
  }
  cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  if (!sites_only) cols <- c(cols, "FORMAT", samples)

  header <- c(
    "##fileformat=VCFv4.2",
    contig_lines,
    '##FILTER=<ID=PASS,Description="All filters passed">',
    if (!sites_only) '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    extra_header,
    paste(cols, collapse = "\t")
  )

  txt <- file.path(dir, paste0(name, ".vcf"))
  writeLines(c(header, lines), txt)
  gz <- file.path(dir, paste0(name, ".vcf.gz"))
  Rsamtools::bgzip(txt, dest = gz, overwrite = TRUE)
  unlink(txt)

  if (index == "tbi") {
    Rsamtools::indexTabix(gz, format = "vcf")
  } else if (index == "csi") {
    status <- system2("bcftools", c("index", "-c", shQuote(gz)), stdout = FALSE, stderr = FALSE)
    if (status != 0) stop("bcftools index -c failed")
  }
  gz
}

has_bcftools <- function() nzchar(Sys.which("bcftools"))

# Records at every availableGeneRanges() position with REF = N * width and
# ALT = A. `chrom` replaces the package chromosome name. Positions are sorted.
#   offtarget adds records outside the gene that must be excluded
gene_lines <- function(gene, build = "GRCh38", chrom = NULL,
                       gt = c("0|1", "1|0", "0|0"), offtarget = FALSE) {
  ref <- availableGeneRanges(gene, build = build)
  chr <- if (is.null(chrom)) as.character(GenomicRanges::seqnames(ref)) else rep(chrom, length(ref))
  pos <- GenomicRanges::start(ref)
  w <- GenomicRanges::width(ref)
  if (offtarget) {
    # Far outside the gene: never returned by the tabix query
    extra_pos <- c(min(pos) - 5000L, max(GenomicRanges::end(ref)) + 5000L)
    extra_w <- c(1L, 1L)
    # Overlapping a PGx range without equalling it: returned by the tabix
    # query, so only the exact-match filter can remove them.
    #  - an SNV one base inside a multi-base reference range
    #  - a 2-bp record ending on a single-base reference position
    cand_pos <- c(pos[w > 1] + 1L, pos[w == 1] - 1L)
    cand_w <- c(rep(1L, sum(w > 1)), rep(2L, sum(w == 1)))
    canon <- as.character(GenomicRanges::seqnames(ref))[1]
    cand <- GenomicRanges::GRanges(canon, IRanges::IRanges(cand_pos, width = cand_w))
    ok <- !IRanges::overlapsAny(cand, GenomicRanges::granges(ref), type = "equal") &
      IRanges::overlapsAny(cand, GenomicRanges::granges(ref))
    pick <- utils::head(which(ok), 3)
    chr <- c(chr, rep(chr[1], 2 + length(pick)))
    pos <- c(pos, extra_pos, cand_pos[pick])
    w <- c(w, extra_w, cand_w[pick])
  }
  o <- order(pos)
  vapply(o, function(i) vrec(chr[i], pos[i], strrep("N", w[i]), "A", gt), character(1))
}

# All PharmVar-defined sites (range, REF, ALT) for a gene
pv_sites <- function(gene, build = "GRCh38") {
  haps <- grep(paste0("^", gene, "\\*"), availableHaplotypes(build), value = TRUE)
  defs <- do.call(rbind, lapply(haps, function(h) {
    v <- availableHaplotypeRanges(h, build = build)
    data.frame(
      hap = h,
      chr = as.character(GenomicRanges::seqnames(v)),
      start = GenomicRanges::start(v),
      end = GenomicRanges::end(v),
      ref = as.character(VariantAnnotation::ref(v)),
      alt = as.character(VariantAnnotation::alt(v)),
      stringsAsFactors = FALSE
    )
  }))
  defs
}

# VCF lines with real PharmVar REF/ALT alleles at every position of a gene.
#   carriers  named list: sample -> c(hap1, hap2); use NA for the reference
#             haplotype, e.g. list(S1 = c("CYP2C9*2", NA))
#   split     one record per ALT (TRUE) or a single multi-ALT record (FALSE)
#   reverse_split  reverse the record order at positions with multiple ALTs
#   chrom     replacement chromosome name
pv_lines <- function(gene, carriers, build = "GRCh38", split = TRUE,
                     reverse_split = FALSE, chrom = NULL) {
  defs <- pv_sites(gene, build)
  sites <- unique(defs[, c("chr", "start", "end", "ref", "alt")])
  sites <- sites[order(sites$start, sites$alt), ]
  has_allele <- function(hap, s) {
    !is.na(hap) && any(defs$hap == hap & defs$start == s$start & defs$end == s$end & defs$alt == s$alt)
  }
  gt_for <- function(s_list, alt_index) {
    vapply(carriers, function(h) {
      a <- vapply(1:2, function(k) {
        idx <- which(vapply(s_list, function(s) has_allele(h[k], s), logical(1)))
        if (length(idx)) as.character(alt_index[idx[1]]) else "0"
      }, character(1))
      paste(a, collapse = "|")
    }, character(1))
  }
  groups <- split(seq_len(nrow(sites)), paste(sites$start, sites$end))
  groups <- groups[order(vapply(groups, function(g) sites$start[g[1]], numeric(1)))]
  out <- character()
  for (g in groups) {
    s_list <- lapply(g, function(i) sites[i, ])
    chr <- if (is.null(chrom)) s_list[[1]]$chr else chrom
    if (split || length(g) == 1) {
      ord <- if (reverse_split) rev(seq_along(s_list)) else seq_along(s_list)
      for (k in ord) {
        s <- s_list[[k]]
        out <- c(out, vrec(chr, s$start, s$ref, s$alt, gt_for(list(s), 1L)))
      }
    } else {
      alts <- vapply(s_list, function(s) s$alt, character(1))
      out <- c(out, vrec(chr, s_list[[1]]$start, s_list[[1]]$ref, paste(alts, collapse = ","),
                         gt_for(s_list, seq_along(s_list))))
    }
  }
  out
}

# Run the full calling pipeline on a VCF and return the call DataFrame
run_pipeline <- function(f, gene, build = "GRCh38") {
  suppressMessages({
    p <- readPGx(f, gene, build = build)
    p <- determineCallableAlleles(p)
    p <- buildReferenceDataFrame(p)
    p <- convertGTtoNucleotides(p)
    callPhasedDiplotypes(p)
  })
}

# Calls as a named character vector (sample -> diplotype)
calls_of <- function(res) stats::setNames(as.character(res[[1]]), rownames(res))

# Path to the large 1000 Genomes test data; skip if unavailable
kgp_dir <- function() {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_not_installed("here")
  d <- here::here("data-raw", "sample_vcf")
  testthat::skip_if_not(dir.exists(d), "1000 Genomes sample VCFs not available")
  d
}

kgp_file <- function(chrom) {
  d <- kgp_dir()
  f <- file.path(d, paste0("1kGP_high_coverage_Illumina.", chrom, ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"))
  if (!file.exists(f) || !file.exists(paste0(f, ".tbi"))) {
    return(NA_character_)
  }
  f
}

# Genes covered by the 1000 Genomes sample VCFs and the chromosome file each
# lives on. Genes whose file is missing locally are dropped (with a skip note).
kgp_genes <- function() {
  df <- data.frame(
    chrom = c("chr19", "chr19", "chr19", "chr10", "chr10", "chr10",
              "chr7", "chr7", "chr19", "chr13", "chr12", "chr15", "chr8"),
    gene = c("CYP2A13", "CYP2A6", "CYP2B6", "CYP2C19", "CYP2C8", "CYP2C9",
             "CYP3A4", "CYP3A5", "CYP4F2", "NUDT15", "SLCO1B1", "CYP1A2", "NAT2"),
    stringsAsFactors = FALSE
  )
  df$file <- vapply(df$chrom, kgp_file, character(1))
  missing <- is.na(df$file)
  if (any(missing)) {
    message("1kGP files missing for: ", paste(df$gene[missing], collapse = ", "), " (skipped)")
  }
  df <- df[!missing, ]
  testthat::skip_if(nrow(df) == 0, "No 1000 Genomes sample VCFs available")
  df
}
