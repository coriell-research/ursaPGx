# Create reference by extracting ranges from PharmVar haplotype.tsv files
#
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

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
    pattern = "*.tsv",
    recursive = TRUE,
    full.names = TRUE
)

grch38_files <- grep("GRCh38", fpaths, value = TRUE)
grch37_files <- grep("GRCh37", fpaths, value = TRUE)
names(grch38_files) <- gsub("\\.NC.*", "", basename(grch38_files))
names(grch37_files) <- gsub("\\.NC.*", "", basename(grch37_files))

grch38_dt <- rbindlist(lapply(grch38_files, fread))
grch37_dt <- rbindlist(lapply(grch37_files, fread))

# Filter for main alleles
grch38_dt <- grch38_dt[!`Haplotype Name` %like% "\\.[0-9]+$"]
grch37_dt <- grch37_dt[!`Haplotype Name` %like% "\\.[0-9]+$"]
grch38_dt <- grch38_dt[!`Haplotype Name` %like% "^rs[0-9]+"]
grch37_dt <- grch37_dt[!`Haplotype Name` %like% "^rs[0-9]+"]
grch38_dt <- grch38_dt[ReferenceSequence != "REFERENCE"]
grch37_dt <- grch37_dt[ReferenceSequence != "REFERENCE"]
grch38_dt <- grch38_dt[!`Haplotype Name` %like% "\\*1$"]
grch37_dt <- grch37_dt[!`Haplotype Name` %like% "\\*1$"]

# Map from RefSeq IDs to UCSC IDs
hg38_alias <- data.table::data.table(
    refseq = c("NC_000001.11", "NC_000010.11", "NC_000011.10", "NC_000012.12",
               "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10",
               "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000002.12",
               "NC_000020.11", "NC_000021.9", "NC_000022.11", "NC_000003.12",
               "NC_000004.12", "NC_000005.10", "NC_000006.12", "NC_000007.14",
               "NC_000008.11", "NC_000009.12", "NC_012920.1", "NC_000023.11",
               "NC_000024.10"),
    ucsc = c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
             "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
             "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
             "chrM", "chrX", "chrY"),
    key = "refseq"
)

hg19_alias <- data.table::data.table(
    refseq = c("NC_000001.10", "NC_000010.10", "NC_000011.9", "NC_000012.11",
               "NC_000013.10", "NC_000014.8", "NC_000015.9", "NC_000016.9",
               "NC_000017.10", "NC_000018.9", "NC_000019.9", "NC_000002.11",
               "NC_000020.10", "NC_000021.8", "NC_000022.10", "NC_000003.11",
               "NC_000004.11", "NC_000005.9", "NC_000006.11", "NC_000007.13",
               "NC_000008.10", "NC_000009.11", "NC_000023.10", "NC_000024.9"),
    ucsc = c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
             "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
             "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
             "chrX", "chrY"),
    key = "refseq"
)
grch38_dt <- hg38_alias[grch38_dt, on = c("refseq"="ReferenceSequence")]
grch37_dt <- hg19_alias[grch37_dt, on = c("refseq"="ReferenceSequence")]

# Relabel insertions and deletions to match VCF standard
grch38_dt[Type == "insertion", ReferenceAllele := "<INS>"]
grch37_dt[Type == "insertion", ReferenceAllele := "<INS>"]
grch38_dt[Type == "deletion", VariantAllele := "<DEL>"]
grch37_dt[Type == "deletion", VariantAllele := "<DEL>"]

# Convert data.tables to GRanges for package
grch38_gr <- GenomicRanges::makeGRangesFromDataFrame(
    grch38_dt,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = "ucsc",
    start.field = "Variant Start",
    end.field = "Variant Stop",
    starts.in.df.are.0based = FALSE,
    seqinfo = GenomeInfoDb::Seqinfo(genome = "hg38")
)
grch37_gr <- GenomicRanges::makeGRangesFromDataFrame(
    grch37_dt,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = "ucsc",
    start.field = "Variant Start",
    end.field = "Variant Stop",
    starts.in.df.are.0based = FALSE,
    seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
)
grch38_gr <- GenomeInfoDb::keepStandardChromosomes(grch38_gr)
grch37_gr <- GenomeInfoDb::keepStandardChromosomes(grch37_gr)

# Split by haplotye and gene
grch38_gene_grl <- split(grch38_gr, f = grch38_gr$Gene)
grch37_gene_grl <- split(grch37_gr, f = grch37_gr$Gene)
grch38_haplotype_grl <- split(grch38_gr, f = grch38_gr$`Haplotype Name`)
grch37_haplotype_grl <- split(grch37_gr, f = grch37_gr$`Haplotype Name`)

# Save package data
usethis::use_data(
  grch38_gene_grl, grch37_gene_grl, grch38_haplotype_grl, grch37_haplotype_grl,
  overwrite = TRUE, internal = TRUE
)
