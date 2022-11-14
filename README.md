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

# To call phased star alleles use the callPhasedDiplotypes() function
df <- callPhasedDiplotypes(p)
>     sample  H1 H2
> 1: HG00096  *2 *1
> 2: HG00097  *1 *1
> 3: HG00099 *17 *1
> 4: HG00100  *1 *1
> 5: HG00101  *1 *1
> 6: HG00102 *17 *1
```

By default `callPhasedDiplotypes()` returns exact matches between the observed
genotypes and allele definitions. The next best match can be returned by setting
the argument `exact = FALSE`

```
df <- callPhasedDiplotypes(p, exact = FALSE)
```

Sub-allele level calls can be returned by setting `summarize = FALSE`

```
df <- callPhasedDiplotypes(p, summarize = FALSE)
>     sample            H1            H2
> 1: HG00096 *2.002,*2.012 *1.002,*1.009
> 2: HG00097 *1.002,*1.009 *1.002,*1.009
> 3: HG00099       *17.001 *1.002,*1.009
> 4: HG00100 *1.002,*1.009 *1.002,*1.009
> 5: HG00101 *1.002,*1.009 *1.002,*1.009
> 6: HG00102       *17.001 *1.002,*1.009
```
