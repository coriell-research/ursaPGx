# Create reference by extracting ranges from PharmVar VCF files
#
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(VariantAnnotation))

# message("Downloading PharmVar Database...")
# url <- "https://www.pharmvar.org/get-download-file?name=ALL&refSeq=ALL&fileType=zip&version=current"
# zipped <- tempfile()
# unzipped <- tempdir()
# download.file(url, destfile = zipped)
#
# message("Extracting PharmVarDB...")
# fpaths <- utils::unzip(zipped, overwrite = TRUE, exdir = unzipped)
fpaths <- list.files(
  "data-raw/pharmvar-5.2.13",
  pattern = "_[0-9]+\\.vcf",
  recursive = TRUE,
  full.names = TRUE
)

grch38_files <- grep("GRCh38", fpaths, value = TRUE)
grch37_files <- grep("GRCh37", fpaths, value = TRUE)
names(grch38_files) <- gsub("\\.vcf", "", basename(grch38_files))
names(grch37_files) <- gsub("\\.vcf", "", basename(grch37_files))

# Read each file into a VRanges object
grch38_haplotype_grl <- lapply(grch38_files, readVcfAsVRanges, genome = Seqinfo(genome = "hg38"))
grch37_haplotype_grl <- lapply(grch37_files, readVcfAsVRanges, genome = Seqinfo(genome = "hg19"))
grch38_haplotype_grl <- lapply(grch38_haplotype_grl, keepStandardChromosomes)
grch37_haplotype_grl <- lapply(grch37_haplotype_grl, keepStandardChromosomes)

# Combine VRanges per gene into a single VRanges object
grch38_haplotypes <- names(grch38_haplotype_grl)
grch37_haplotypes <- names(grch37_haplotype_grl)
names(grch38_haplotypes) <- gsub("_[0-9]+$", "", grch38_haplotypes)
names(grch37_haplotypes) <- gsub("_[0-9]+$", "", grch37_haplotypes)
grch38_genes <- unique(names(grch38_haplotypes))
grch37_genes <- unique(names(grch37_haplotypes))

getGrs <- function(x, hg38 = TRUE) {
  if (hg38) {
    l <- grch38_haplotype_grl[grch38_haplotypes[names(grch38_haplotypes) == x]]
  } else {
    l <- grch37_haplotype_grl[grch37_haplotypes[names(grch37_haplotypes) == x]]
  }
  names(l) <- NULL
  unique(do.call(c, l))
}

grch38_gene_grl <- lapply(grch38_genes, getGrs)
grch37_gene_grl <- lapply(grch37_genes, getGrs, hg38 = FALSE)
names(grch38_gene_grl) <- grch38_genes
names(grch37_gene_grl) <- grch37_genes

# Save package data
usethis::use_data(
  grch38_gene_grl, grch37_gene_grl, grch38_haplotype_grl, grch37_haplotype_grl,
  overwrite = TRUE, internal = TRUE
)
