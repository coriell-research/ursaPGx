## ursaPGx <img src='man/figures/logo.png' align="right" height="139" />

The goal of this package is to use phased VCF data to assign star alleles 
to samples using existing frameworks from the 
[Bioconductor](https://www.bioconductor.org/) ecosystem and 
[PharmVar](https://www.pharmvar.org) database.

## Data Sources

Reference alleles and haplotype definitions are extracted from the most recent
version of [PharmVar](https://www.pharmvar.org/download). See the 
`create-reference.R` function in the data-raw directory for the exact script.

## Installation

This package is still in active development but can be installed with:

`devtools::install_github("coriell-research/ursaPGx")`

## Example Usage

```r
library(ursaPGx)


# Use 1kGP VCF file as example (.tbi index must be located in same directory)
vcf <- "1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# Read in PGx data for CYP2C19
p <- readPGx(vcf, "CYP2C19")

# Display Class information
p
> Class: PGx
> PGx gene: CYP2C19 
> PGx build: GRCh38 
> Number of samples: 3202
```
