#  Create Haplotype reference GRanges from PharmVar db
#
# The script pulls whatever is the latest version of the Pharmvar db and creates
# reference files that get used internally by the package.
# ------------------------------------------------------------------------------
library(data.table)
library(GenomicRanges)
library(GenomeInfoDb)


# message("Downloading PharmVar Database...")
# url <- "https://www.pharmvar.org/get-download-file?name=ALL&refSeq=ALL&fileType=zip&version=current"
# zipped <- tempfile()
# unzipped <- tempdir()
# download.file(url, destfile = zipped)
#
# message("Extracting PharmVarDB...")
# fpaths <- utils::unzip(zipped, overwrite = TRUE, exdir = unzipped)
fpaths <- list.files(
  "data-raw/pharmvar-5.2.13/",
  pattern = "haplotypes.tsv",
  recursive = TRUE,
  full.names = TRUE
  )

message("Reading haplotype definitions...")
haplotypes <- fpaths[grepl("haplotypes.tsv", fpaths)]
grch38_files <- haplotypes[grepl("GRCh38", haplotypes)]
grch37_files <- haplotypes[grepl("GRCh37", haplotypes)]
names(grch38_files) <- gsub("\\..*haplotypes\\.tsv", "", basename(grch38_files))
names(grch37_files) <- gsub("\\..*haplotypes\\.tsv", "", basename(grch37_files))

# Bind into single table first to clean, then split into list again
grch38_dt <- rbindlist(lapply(grch38_files, fread))
grch37_dt <- rbindlist(lapply(grch37_files, fread))
names(grch38_dt) <- gsub(" ", "", names(grch38_dt))
names(grch37_dt) <- gsub(" ", "", names(grch37_dt))
grch38_dt <- grch38_dt[ReferenceSequence != "REFERENCE"]
grch37_dt <- grch37_dt[ReferenceSequence != "REFERENCE"]
setkey(grch38_dt, "ReferenceSequence")
setkey(grch37_dt, "ReferenceSequence")

# Recode the insertions and deletions to match VCF spec
grch38_dt[ReferenceAllele == "-", ReferenceAllele := "<INS>"]
grch38_dt[VariantAllele == "-", VariantAllele := "<DEL>"]
grch37_dt[ReferenceAllele == "-", ReferenceAllele := "<INS>"]
grch37_dt[VariantAllele == "-", VariantAllele := "<DEL>"]

# Define the chromosome alias map
# https://hgdownload.cse.ucsc.edu/goldenPath/hg{19|38}/database/chromAlias.txt.gz
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

# Join on UCSC names
grch38_def <- hg38_alias[grch38_dt]
grch37_def <- hg19_alias[grch37_dt]

# Create new ID column
grch38_def[, id := paste(ucsc, VariantStart, VariantStop, ReferenceAllele, VariantAllele, sep = ".")]
grch37_def[, id := paste(ucsc, VariantStart, VariantStop, ReferenceAllele, VariantAllele, sep = ".")]

# Split
grch38_haplotypes <- split(grch38_def, by = "HaplotypeName")
grch37_haplotypes <- split(grch37_def, by = "HaplotypeName")
grch38_genes <- split(grch38_def, by = "Gene")
grch37_genes <- split(grch37_def, by = "Gene")

# Create GRanges objects for each haplotype
grch38_haplotype_grl <- lapply(
  grch38_haplotypes,
  makeGRangesFromDataFrame,
  seqnames.field = "ucsc",
  start.field = "VariantStart",
  end.field = "VariantStop",
  ignore.strand = TRUE,
  seqinfo = Seqinfo(genome = "hg38"),
  keep.extra.columns = TRUE
  )
grch38_haplotype_grl <- lapply(grch38_haplotype_grl, keepStandardChromosomes)

grch37_haplotype_grl <- lapply(
  grch37_haplotypes,
  makeGRangesFromDataFrame,
  seqnames.field = "ucsc",
  start.field = "VariantStart",
  end.field = "VariantStop",
  ignore.strand = TRUE,
  seqinfo = Seqinfo(genome = "hg19"),
  keep.extra.columns = TRUE
)
grch37_haplotype_grl <- lapply(grch37_haplotype_grl, keepStandardChromosomes)

# Create GRanges for each gene
grch38_gene_grl <- lapply(
  grch38_genes,
  makeGRangesFromDataFrame,
  seqnames.field = "ucsc",
  start.field = "VariantStart",
  end.field = "VariantStop",
  ignore.strand = TRUE,
  seqinfo = Seqinfo(genome = "hg38"),
  keep.extra.columns = TRUE
)
grch38_gene_grl <- lapply(grch38_gene_grl, keepStandardChromosomes)

grch37_gene_grl <- lapply(
  grch37_genes,
  makeGRangesFromDataFrame,
  seqnames.field = "ucsc",
  start.field = "VariantStart",
  end.field = "VariantStop",
  ignore.strand = TRUE,
  seqinfo = Seqinfo(genome = "hg19"),
  keep.extra.columns = TRUE
)
grch37_gene_grl <- lapply(grch37_gene_grl, keepStandardChromosomes)

# Save data for package use
usethis::use_data(
  grch38_def, grch37_def,
  grch38_gene_grl, grch37_gene_grl,
  grch38_haplotype_grl, grch37_haplotype_grl,
  overwrite = TRUE, internal = TRUE
)
