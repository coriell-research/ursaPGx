#' Locate the index for a bgzipped VCF file
#'
#' Returns the path to a tabix (.tbi) or CSI (.csi) index for a local file, or
#' NULL for remote files (where Rsamtools resolves the index itself).
#' @noRd
.findVcfIndex <- function(file, remote) {
  if (remote) {
    return(NULL)
  }
  candidates <- paste0(file, c(".tbi", ".csi"))
  idx <- candidates[file.exists(candidates)]
  if (length(idx) == 0) {
    stop(
      "No index found for '", file, "'. Expected '", candidates[1], "' or '",
      candidates[2], "'. Create one with `tabix -p vcf <file>` or ",
      "`bcftools index -t <file>`.",
      call. = FALSE
    )
  }
  idx[1]
}

#' Check that a local file is BGZF compressed
#'
#' Reading a plain-gzip (or uncompressed) file through a tabix index can crash
#' the R session, so verify the BGZF magic bytes before touching htslib.
#' @noRd
.isBgzf <- function(file) {
  con <- file(file, "rb")
  on.exit(close(con))
  b <- readBin(con, "raw", n = 18)
  length(b) >= 14 &&
    identical(b[1:4], as.raw(c(0x1f, 0x8b, 0x08, 0x04))) &&
    identical(b[13:14], as.raw(c(0x42, 0x43)))
}

#' Candidate chromosome names for a canonical package seqlevel
#'
#' Given a seqlevel used by the package reference ranges (e.g. "chr10" for
#' GRCh38 or "10" for GRCh37), return the names that the same chromosome may
#' carry in a user VCF: UCSC ("chr10"), NCBI/Ensembl ("10"), dbSNP ("ch10") and
#' the RefSeq accession without version ("NC_000010").
#' @noRd
.seqlevelCandidates <- function(level) {
  core <- sub("^chr", "", level, ignore.case = TRUE)
  cores <- if (toupper(core) %in% c("M", "MT")) c("M", "MT") else core
  refseq <- if (toupper(core) %in% c("M", "MT")) {
    "NC_012920"
  } else {
    num <- switch(toupper(core), X = 23L, Y = 24L, suppressWarnings(as.integer(core)))
    if (is.na(num)) character() else sprintf("NC_%06d", num)
  }
  list(names = unique(c(paste0("chr", cores), cores, paste0("ch", cores))), refseq = refseq)
}

#' Map the reference seqlevels onto the chromosome names used in the VCF
#'
#' The tabix index lists the sequences that actually contain records, so it is
#' used as the source of truth instead of the (optional, possibly stale) VCF
#' header. Each seqlevel in use by \code{ref} is matched case-insensitively to
#' one of its UCSC, NCBI, dbSNP or RefSeq spellings in the index.
#'
#' @return A named character vector mapping canonical reference seqlevels
#' (names) to the seqlevels used by the VCF file (values).
#' @noRd
.mapSeqlevelsToFile <- function(tab, ref, gene) {
  idx <- Rsamtools::seqnamesTabix(tab)
  idx_unversioned <- sub("\\.[0-9]+$", "", idx)
  in_use <- GenomeInfoDb::seqlevelsInUse(ref)

  map <- vapply(in_use, function(level) {
    cand <- .seqlevelCandidates(level)
    hits <- idx[tolower(idx) %in% tolower(cand$names) | idx_unversioned %in% cand$refseq]
    if (length(hits) == 0) {
      present <- if (length(idx) == 0) "none" else paste(utils::head(idx, 10), collapse = ", ")
      stop(
        "The VCF does not contain any records on chromosome '", level,
        "', which is required for ", gene, ". Sequences with records in the ",
        "VCF index: ", present, if (length(idx) > 10) ", ...", ".",
        call. = FALSE
      )
    }
    if (length(hits) > 1) {
      # Prefer the exact package spelling, otherwise take the first match
      exact <- hits[hits == level]
      if (length(exact) == 1) {
        return(exact)
      }
      message(
        "Multiple sequences in the VCF match chromosome '", level, "': ",
        paste(hits, collapse = ", "), ". Using '", hits[1], "'."
      )
    }
    hits[1]
  }, character(1))

  if (any(map != names(map))) {
    message(
      "Mapping reference chromosome names to VCF chromosome names: ",
      paste(names(map), "->", map, collapse = ", ")
    )
  }

  map
}

#' Check that VCF header contig lengths agree with the requested build
#'
#' A common user error is reading a GRCh37 VCF with the default GRCh38 build,
#' which otherwise produces an empty result with only a generic warning. When
#' the header declares contig lengths, compare them against the reference.
#' @noRd
.checkContigLengths <- function(file, ref, map, build) {
  hdr <- tryCatch(VariantAnnotation::scanVcfHeader(file), error = function(e) NULL)
  if (is.null(hdr)) {
    return(invisible(TRUE))
  }
  hdr_len <- tryCatch(
    GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(hdr)),
    error = function(e) NULL
  )
  if (is.null(hdr_len)) {
    return(invisible(TRUE))
  }
  ref_len <- GenomeInfoDb::seqlengths(ref)

  for (canonical in names(map)) {
    file_level <- map[[canonical]]
    hl <- if (file_level %in% names(hdr_len)) hdr_len[[file_level]] else NA
    rl <- ref_len[[canonical]]
    if (!is.na(hl) && !is.na(rl) && hl != rl) {
      other <- if (build == "GRCh38") "GRCh37" else "GRCh38"
      stop(
        "Contig '", file_level, "' has length ", hl, " in the VCF header but ",
        "length ", rl, " in the ", build, " reference. Is this VCF aligned to ",
        other, "? If so, use `build = \"", other, "\"`.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' Require samples and a GT field in the imported VCF
#' @noRd
.checkHasGenotypes <- function(vcf) {
  if (ncol(vcf) == 0) {
    stop("The VCF contains no samples. readPGx() requires sample genotypes (GT).", call. = FALSE)
  }
  if (is.null(VariantAnnotation::geno(vcf)$GT)) {
    stop("The VCF does not contain a GT FORMAT field. readPGx() requires sample genotypes.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Report unphased, missing or haploid genotypes at retained PGx positions
#' @noRd
.checkGenotypeQuality <- function(vcf) {
  gt <- VariantAnnotation::geno(vcf)$GT
  if (length(gt) == 0) {
    return(invisible(TRUE))
  }
  n_missing <- sum(grepl(".", gt, fixed = TRUE))
  n_unphased <- sum(grepl("/", gt, fixed = TRUE))
  n_haploid <- sum(!grepl("[|/]", gt) & !grepl(".", gt, fixed = TRUE))
  problems <- c(
    if (n_unphased > 0) paste(n_unphased, "unphased"),
    if (n_missing > 0) paste(n_missing, "missing"),
    if (n_haploid > 0) paste(n_haploid, "haploid")
  )
  if (length(problems) > 0) {
    warning(
      "Genotypes at PGx positions include ", paste(problems, collapse = ", "),
      " call(s) out of ", length(gt), ". Diplotype calling expects phased, ",
      "diploid, non-missing genotypes; affected samples may be called '*Amb'.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Drop copies of a record returned for more than one query block
#'
#' Tabix returns a record once for every query region it overlaps. Records
#' whose span (REF length or INFO/END) crosses several reduced reference blocks
#' are therefore returned several times. Copies coming from different query
#' blocks are removed; repeats within the same block are genuine duplicates in
#' the input and are kept.
#' @noRd
.dropSpanDuplicates <- function(vcf) {
  rr <- SummarizedExperiment::rowRanges(vcf)
  if (length(rr) > 1 && "paramRangeID" %in% names(S4Vectors::mcols(rr))) {
    pid <- as.character(rr$paramRangeID)
    gt <- VariantAnnotation::geno(vcf)$GT
    gt_key <- if (is.null(gt) || ncol(gt) == 0) "" else apply(gt, 1, paste, collapse = "\t")
    key <- paste(
      names(rr), as.character(GenomicRanges::seqnames(rr)), GenomicRanges::start(rr),
      GenomicRanges::end(rr), as.character(VariantAnnotation::ref(vcf)),
      .altString(vcf), gt_key, sep = "\r"
    )
    first_pid <- pid[match(key, key)]
    vcf <- vcf[!(duplicated(key) & pid != first_pid), ]
  }
  rr <- SummarizedExperiment::rowRanges(vcf)
  S4Vectors::mcols(rr)$paramRangeID <- NULL
  SummarizedExperiment::rowRanges(vcf) <- rr
  vcf
}

#' Comma-separated ALT alleles for each record
#' @noRd
.altString <- function(vcf) {
  vapply(VariantAnnotation::alt(vcf), function(a) paste(as.character(a), collapse = ","), character(1))
}

#' Rebuild VariantAnnotation-style rownames using canonical seqlevels
#'
#' readVcf() names records by ID, or by "CHROM:POS_REF/ALT" when ID is missing.
#' Generated names carry the file's chromosome spelling, so rewrite their
#' prefix and ensure names are unique.
#' @noRd
.canonicalRownames <- function(nms, rr, map) {
  if (length(nms) == 0) {
    return(nms)
  }
  canon <- as.character(GenomicRanges::seqnames(rr))
  pos <- GenomicRanges::start(rr)
  inverse <- structure(names(map), names = unname(map))
  for (file_level in names(inverse)) {
    prefix <- paste0(file_level, ":", pos, "_")
    generated <- startsWith(nms, prefix) & canon == inverse[[file_level]]
    nms[generated] <- paste0(
      inverse[[file_level]], ":", pos[generated], "_",
      substring(nms[generated], nchar(prefix[generated]) + 1)
    )
  }
  make.unique(nms, sep = "_")
}

#' Read in PGx Data from VCF File
#'
#' This function is a wrapper around \code{\link[VariantAnnotation]{readVcf}} that
#' reads only the PGx locations into a PGx object. The user must provide a
#' bgzipped and indexed VCF file. The path to the indexed VCF file gets passed to
#' \code{\link[Rsamtools]{TabixFile}} which is then used as the file argument to
#' \code{readVcf()}.
#'
#' @details The PGx locations used to subset the input VCF file are derived from
#' PharmVar VCF files available at \url{https://www.pharmvar.org/download}.
#' The "create-reference.R" file included in the data-raw directory of
#' this package describes the steps used to reproduce the construction of the
#' GRanges object used in this function. These GRanges objects are accessible
#' via the \code{availableGeneRanges()} and \code{availableHaplotypeRanges()}
#' functions.
#'
#' Chromosome names are detected from the sequences listed in the VCF index
#' (not the VCF header), so UCSC ("chr10"), NCBI/Ensembl ("10"), dbSNP
#' ("ch10") and RefSeq ("NC_000010.11") naming are all supported, matched
#' case-insensitively. If the VCF header declares contig lengths that disagree
#' with the requested build, an error is raised.
#'
#' Only records whose ranges exactly equal a PharmVar position are retained.
#' Records that overlap a PharmVar position without matching it exactly (for
#' example indels that are not left-aligned and normalized) are dropped and
#' reported with a message. Input VCFs should be normalized, e.g. with
#' \code{bcftools norm -m -any -f ref.fa}.
#'
#' The returned object always uses the package's canonical sequence
#' information: "UCSC" style seqlevels and the hg38 genome for GRCh38, and
#' "NCBI" style seqlevels and the GRCh37 genome for GRCh37. Rows are sorted by
#' position.
#'
#' @param file Path to a bgzipped, phased, and normalized VCF file. The file
#' must have a tabix (.tbi) or CSI (.csi) index in the same directory.
#' @param gene The PGx gene to subset from VCF file. See
#' \code{availableGenes()} for available genes.
#' @param build The genome build. One of "GRCh38" or "GRCh37".
#' @export
#' @return Object of class PGx
#' @examples
#' \dontrun{
#' # Specify the path to the indexed VCF file
#' vcf <- "/path/to/vcf.gz"
#'
#' # Create a PGx object by specifying the desired gene
#' p <- readPGx(vcf, "CYP2C9")
#' }
readPGx <- function(file, gene, build = c("GRCh38", "GRCh37")) {
  # Argument checks ------------------------------------------------------------
  stopifnot("Only a single file can be used as input" = length(file) == 1)
  if (!is.character(file) || is.na(file) || !nzchar(file)) {
    stop("`file` must be a single path to a bgzipped, indexed VCF file.", call. = FALSE)
  }

  builds <- c("GRCh38", "GRCh37")
  if (identical(build, builds)) {
    build <- builds[1]
  }
  if (!is.character(build) || length(build) != 1 || is.na(build) || !build %in% builds) {
    stop("Genome must be one of 'GRCh38' or 'GRCh37'", call. = FALSE)
  }

  if (length(gene) != 1 || is.na(gene)) {
    stop("Exactly one gene must be supplied.", call. = FALSE)
  }
  if (!gene %in% availableGenes(build)) {
    stop("Gene not in gene list. See `availableGenes()` for supported genes.", call. = FALSE)
  }
  if (gene == "CYP2D6") {
    stop("Please use `cyrius()` to perform calling for CYP2D6.", call. = FALSE)
  }

  ext <- tolower(tools::file_ext(file))
  if (ext == "vcf") {
    stop("File extension is '.vcf'. Are you sure this is an indexed VCF file?", call. = FALSE)
  }
  if (ext == "bcf") {
    stop(
      "BCF input is not supported. Convert to a bgzipped VCF with ",
      "`bcftools view -Oz -o out.vcf.gz in.bcf` and index it with `tabix -p vcf out.vcf.gz`.",
      call. = FALSE
    )
  }

  remote <- grepl("^(https?|ftp|s3)://", file, ignore.case = TRUE)
  if (!remote) {
    if (!file.exists(file)) {
      stop("VCF file not found: '", file, "'.", call. = FALSE)
    }
    if (!.isBgzf(file)) {
      stop(
        "'", file, "' is not BGZF compressed. Compress it with `bgzip` ",
        "(not `gzip`) and index it with `tabix -p vcf`.",
        call. = FALSE
      )
    }
  }
  index <- .findVcfIndex(file, remote)

  # Resolve chromosome names ---------------------------------------------------
  ref <- availableGeneRanges(gene, build = build)
  canonical_si <- GenomeInfoDb::seqinfo(ref)
  genome <- unique(GenomeInfoDb::genome(ref))[1]

  tab <- if (is.null(index)) Rsamtools::TabixFile(file) else Rsamtools::TabixFile(file, index = index)
  map <- .mapSeqlevelsToFile(tab, ref, gene)
  .checkContigLengths(file, ref, map, build)

  ref_in_use <- GenomeInfoDb::keepSeqlevels(ref, names(map))
  ref_file <- GenomeInfoDb::renameSeqlevels(ref_in_use, map)

  # Read the VCF ---------------------------------------------------------------
  # Query with reduced ranges: overlapping reference ranges would otherwise
  # cause tabix to return the same record once per overlapping range. Blocks
  # are named so records spanning several blocks can be de-duplicated.
  query <- GenomicRanges::reduce(ref_file)
  names(query) <- paste0("q", seq_along(query))
  # readVcf() warns when header contigs share no names with the query ranges
  # (e.g. a stale header); seqinfo is replaced below, so that warning is noise.
  vcf <- withCallingHandlers(
    VariantAnnotation::readVcf(tab, genome = genome, param = query),
    warning = function(w) {
      if (grepl("no sequence levels in common", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  .checkHasGenotypes(vcf)
  vcf <- .dropSpanDuplicates(vcf)

  # Convert to canonical seqinfo -----------------------------------------------
  # seqlevels<- (unlike keepSeqlevels) also adds levels absent from an empty result
  GenomeInfoDb::seqlevels(vcf, pruning.mode = "coarse") <- unname(map)
  vcf <- GenomeInfoDb::renameSeqlevels(vcf, structure(names(map), names = unname(map)))
  GenomeInfoDb::seqinfo(vcf, new2old = match(
    GenomeInfoDb::seqlevels(canonical_si), GenomeInfoDb::seqlevels(vcf)
  )) <- canonical_si

  # Keep only records exactly matching a reference position --------------------
  rr <- SummarizedExperiment::rowRanges(vcf)
  ov <- GenomicRanges::findOverlaps(query = rr, subject = ref, type = "equal")
  keep <- sort(unique(S4Vectors::queryHits(ov)))

  # Report likely representation mismatches: indels overlapping a PGx position,
  # or any record starting at a PGx position, that do not match it exactly.
  # (SNVs that merely fall inside a PharmVar indel span are not reported.)
  alt_chr <- .altString(vcf)
  symbolic <- grepl("^<", alt_chr)
  is_indel <- GenomicRanges::width(rr) != 1 |
    vapply(VariantAnnotation::alt(vcf), function(a) any(nchar(as.character(a)) != 1), logical(1))
  candidate <- IRanges::overlapsAny(rr, ref) & !symbolic &
    (is_indel | IRanges::overlapsAny(rr, ref, type = "start"))
  near_miss <- setdiff(which(candidate), keep)
  if (length(near_miss) > 0) {
    ex <- paste0(
      as.character(GenomicRanges::seqnames(rr[near_miss])), ":",
      GenomicRanges::start(rr[near_miss]), "_",
      as.character(VariantAnnotation::ref(vcf)[near_miss]), "/", alt_chr[near_miss]
    )
    message(
      length(near_miss), " record(s) overlap PGx positions but do not exactly ",
      "match a PharmVar definition and were dropped: ",
      paste(utils::head(ex, 5), collapse = ", "), if (length(ex) > 5) ", ...",
      ". If these are indels, ensure the VCF is normalized and left-aligned ",
      "(e.g. `bcftools norm -m -any -f ref.fa`)."
    )
  }

  vcf_filtered <- vcf[keep, ]
  vcf_filtered <- vcf_filtered[order(SummarizedExperiment::rowRanges(vcf_filtered)), ]
  .checkGenotypeQuality(vcf_filtered)

  # Warn about genuinely duplicated input records
  rr_f <- SummarizedExperiment::rowRanges(vcf_filtered)
  key <- paste(
    as.character(rr_f), as.character(VariantAnnotation::ref(vcf_filtered)),
    .altString(vcf_filtered)
  )
  if (anyDuplicated(key) > 0) {
    warning(
      "The VCF contains ", sum(duplicated(key)), " duplicated record(s) ",
      "(same position, REF and ALT) at PGx positions: ",
      paste(utils::head(unique(key[duplicated(key)]), 5), collapse = "; "),
      ". Duplicated records were kept; consider deduplicating the input ",
      "(e.g. `bcftools norm -d exact`).",
      call. = FALSE
    )
  }

  rownames(vcf_filtered) <- .canonicalRownames(rownames(vcf_filtered), rr_f, map)

  p <- PGx(
    vcf = vcf_filtered,
    pgxGene = gene,
    pgxBuild = build,
    pgxReferenceDataFrame = S4Vectors::DataFrame(),
    hasCallableAlleles = FALSE,
    hasReferenceDataFrame = FALSE
  )

  # Issue a warning if there are no overlapping ranges
  if (nrow(p) == 0) {
    warning("No overlapping ranges were detected between the sample VCF and the reference.")
  }

  return(p)
}
