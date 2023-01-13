---
editor_options: 
  markdown: 
    wrap: 79
---

## ursaPGx <img src="man/figures/logo.png" align="right" height="139"/>

The goal of this package is to use phased VCF data to assign star alleles to
samples using existing frameworks from the
[Bioconductor](https://www.bioconductor.org/) ecosystem and
[PharmVar](https://www.pharmvar.org) database.

## Data Sources

Reference alleles and haplotype definitions are extracted from the most recent
version of [PharmVar](https://www.pharmvar.org/download). See the
`create-reference.R` function in the data-raw directory for the exact script.

## Installation

This package is still in active development but can be installed with:

    # Install the required Bioconductor packages
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    # Install ursaPGx from github
    devtools::install_github("coriell-research/ursaPGx")

## Example Usage

The simplest use case (and currently the only implemented usage) is for calling
phased diplotypes. Below, we read in the VCF file as a `PGx` object and call
the wrapper method `callPhasedDiplotypes()` on it. The result is a `DataFrame`
with a row for every sample in the VCF and a column for every callable defined
allele.

``` r
library(ursaPGx)


# Use 1kGP VCF file as example (.tbi index must be located in same directory)
vcf <- "1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# Read in PGx data for CYP2C8
p <- readPGx(vcf, "CYP2C8")

# Run the phased diplotype caller pipeline for CYP2C8
df <- callPhasedDiplotypes(p)

# Optionally add the results back to the colData slot of the PGx object
colData(p) <- df

# View the results
head(colData(p))
```

The `callPhasedDiplotypes()` method is a wrapper around several functions
useful for generating diplotype calls for all defined alleles. If you wish to
examine all individual steps of the pipeline that can be done like so:

``` r
# 1. Determine which of the defined alleles are fully represented in the sample VCF
p <- getCallableAlleles(p)

# Use the getter method on the PGx object to retreive the full list
callableAlleles(p)

# For each callable allele convert the genotype matrix to nucleotides
# Only CYP2C8*3 is shown below
p3 <- pgxGenotypeCodesToNucleotides(p, "CYP2C8_3")

# 3. Generate CYP2C8*11 allele calls for all samples
df3 <- callPhasedDiplotype(p3)
```

## Accessing the definitions

Allele definitions are downloaded from PharmVar and converted into `VRanges`
objects. The user can access all of the allele definitions with the following
functions:

``` r
# List all available PGx genes with definitions
availableGenes()

# List all available star alleles with definitions
availableHaplotypes()

# Access the VRanges for a given gene
availableGeneRanges("CYP2C8")

# Access VRanges for the given the star allele
availableHaplotypeRanges("CYP2C8_3")
```

## Working with the `PGx` object

`PGx` objects inherit from `VariantAnnotation::VCF` object so anything you can
do with a `VCF` object you can do with a `PGx` object. See the
[VariantAnnotation
documentation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
for more details.

``` r
# Contains all sample-level data
colData(PGx)

# Contains information about the genomic positions
rowRanges(PGx)

# Converts the genotype calls into a matrix
gt <- geno(PGx)$GT

# Subsets the PGx object by the CYP2C8*3 defined positions
cyp2c8_3 <- subsetByOverlaps(
    p, 
    availableHaplotypeRanges("CYP2C8_3"), 
    type = "equal"
    )
```
