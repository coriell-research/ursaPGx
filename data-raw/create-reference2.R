# Create Reference Ranges from PharmVar database
#
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(VariantAnnotation))


# File operations ---------------------------------------------------------
fpaths <- list.files(
    path = "data-raw/pharmvar-5.2.13", 
    pattern = "*.vcf",
    recursive = TRUE,
    full.names = TRUE
)

# Select GRCh38 and GRCh37 files separately
grch37_vcf <- fpaths[grepl("GRCh37", fpaths)]
grch38_vcf <- fpaths[grepl("GRCh38", fpaths)]

# Keep only major allele definitions
grch37_vcf <- grch37_vcf[!grepl("\\.[0-9]+\\.vcf$", grch37_vcf)]
grch38_vcf <- grch38_vcf[!grepl("\\.[0-9]+\\.vcf$", grch38_vcf)]

# Name by the haplotype
names(grch37_vcf) <- tools::file_path_sans_ext(basename(grch37_vcf))
names(grch38_vcf) <- tools::file_path_sans_ext(basename(grch38_vcf))

# Swap "_" for "*"
names(grch37_vcf) <- gsub("_", "*", names(grch37_vcf))
names(grch38_vcf) <- gsub("_", "*", names(grch38_vcf))

# Read in data as ranges --------------------------------------------------

# Read in separate lists for each genome build
grch37_haplotype_ranges <- parallel::mclapply(
    grch37_vcf,
    readVcfAsVRanges,
    genome = "hg19",
    mc.cores = 4
)
grch38_haplotype_ranges <- parallel::mclapply(
    grch38_vcf,
    readVcfAsVRanges,
    genome = "hg38",
    mc.cores = 4
)

# Keep only the standard chromosome information
grch37_haplotype_ranges <- parallel::mclapply(
    grch37_haplotype_ranges,
    keepStandardChromosomes,
    mc.cores = 4
)
grch38_haplotype_ranges <- parallel::mclapply(
    grch38_haplotype_ranges,
    keepStandardChromosomes,
    mc.cores = 4
)

# Coerce to GRangesLists
grch37_haplotype_grl <- GRangesList(grch37_haplotype_ranges)
grch38_haplotype_grl <- GRangesList(grch38_haplotype_ranges)

# Get per gene ranges ------------------------------------------------------

# Extract all unique gene names from the VCF file names
grch37_genes <- gsub("data-raw/pharmvar\\-[0-9]+\\.[0-9]+\\.[0-9]+/", "", grch37_vcf)
grch37_genes <- gsub("\\/GRCh37.*", "", grch37_genes)
grch37_uniq_genes <- unique(grch37_genes)

grch38_genes <- gsub("data-raw/pharmvar\\-[0-9]+\\.[0-9]+\\.[0-9]+/", "", grch38_vcf)
grch38_genes <- gsub("\\/GRCh38.*", "", grch38_genes)
grch38_uniq_genes <- unique(grch38_genes)

# Collapse haplotype ranges into single gene ranges
grch37_gene_gr <- vector("list", length(grch37_uniq_genes))
grch38_gene_gr <- vector("list", length(grch38_uniq_genes))
names(grch37_gene_gr) <- grch37_uniq_genes
names(grch38_gene_gr) <- grch38_uniq_genes

for (g in grch37_uniq_genes) {
    keep <- names(grch37_genes[grch37_genes == g])
    grch37_gene_gr[[g]] <- grch37_haplotype_ranges[keep]
}
for (g in grch38_uniq_genes) {
    keep <- names(grch38_genes[grch38_genes == g])
    grch38_gene_gr[[g]] <- grch38_haplotype_ranges[keep]
}

# Create a list of GRangesLists for each gene
grch37_gene_grl <- lapply(grch37_gene_gr, GRangesList)
grch38_gene_grl <- lapply(grch38_gene_gr, GRangesList)

# Save as package data ----------------------------------------------------

usethis::use_data(
    grch38_gene_grl, grch37_gene_grl, grch38_haplotype_grl, grch37_haplotype_grl,
    overwrite = TRUE, internal = TRUE
)
