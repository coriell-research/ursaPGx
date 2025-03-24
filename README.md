## ursaPGx <img src="man/figures/logo.png" align="right" height="139"/>

The goal of this package is to use phased VCF data to assign star alleles to
samples using existing frameworks from the
[Bioconductor](https://www.bioconductor.org/) ecosystem and
[PharmVar](https://www.pharmvar.org) database. This package was purpose-built for
annotating the 1000 Genomes Project 30X phased VCF call sets from the [NYGC](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) but we have designed the functions to be generally applicable to 
other datasets.

### Publication

The details of the ursaPGx star-calling algorithm along with benchmark results are 
presented in our [manuscript](https://www.frontiersin.org/articles/10.3389/fbinf.2024.1351620/full)

## Data Sources

Reference alleles and haplotype definitions are extracted from [PharmVar](https://www.pharmvar.org/download). See the
`create-reference.R` function in the data-raw directory for the exact script.

The current version of the reference haplotypes from PharmVar is: **Version 6.2.3**

## Installation

This package is still in active development and can be installed with:

```r
if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pak("coriell-research/ursaPGx")
```

## Quick Start

The `callDiplotypes()` function is a wrapper for calling the main pipeline 
steps and returning a DataFrame of the allele calls. Generating allele calls
for all samples in a VCF file for CYP2C8, for example, can be done with:

```r
# Specify the path to the VCF file
vcf <- "1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# Call phased diplotypes for CYP2C8
result <- callDiplotypes(vcf, gene = "CYP2C8", phased = TRUE)

result
>DataFrame with 3202 rows and 1 column
>             CYP2C8
>        <character>
>HG00096       *4|*1
>HG00097       *1|*1
>HG00099       *1|*1
>HG00100       *1|*1
>HG00101       *1|*1
>...             ...
>NA21137       *1|*1
>NA21141       *1|*1
>NA21142       *1|*3
>NA21143       *1|*1
>NA21144       *1|*1
```

## Full pipeline

Each of steps wrapped in the `callDiplotypes()` function above can be run 
individually so that the results of each step can be checked. The full caller 
pipeline for CYP2C19, for example:

```r
# Specify the path the the VCF object
vcf <- "1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# Read in the VCF data as a PGx object for CYP2C19
CYP2C19 <- readPGx(vcf, gene = "CYP2C19")

# Determine what alleles can be called from the data
CYP2C19 <- determineCallableAlleles(CYP2C19)

# To return a vector of the callable alleles
pgxCallableAlleles(CYP2C19)

# Create a reference of all positions from the callable alleles 
CYP2C19 <- buildReferenceDataFrame(CYP2C19)

# To return the reference DataFrame
pgxReferenceDataFrame(CYP2C19)

# Convert the genotype code to nucleotides
CYP2C19 <- convertGTtoNucleotides(CYP2C19)

# To return the genotype matrix for all samples
pgxGenotypeMatrix(CYP2C19)

# Create diplotype calls for every sample
result <- callPhasedDiplotypes(CYP2C19)
head(result)

>        CYP2C19
>HG00096   *2|*1
>HG00097   *1|*1
>HG00099  *17|*1
>HG00100   *1|*1
>HG00101   *1|*1
>HG00102  *17|*1

# Examine exactly which positions matched all samples and all haplotypes
details <- detailPhasedCalls(CYP2C19)

# Returns a position by star allele matrix for all variants on the first 
# haplotype of sample HG00096
details[["HG00096"]][["H1"]]

# Checking this result yields the same call as above
which(apply(details[["HG00096"]][["H1"]], 2, all))  # CYP2C19*2
```

## CYP2D6 

CYP2D6 allele calling is performed using an interface to [Ilumina Cyrius CYP2D6 star allele caller](https://github.com/Illumina/Cyrius). Since CYP2D6 calling needs copy 
number information, BAM/CRAM files are used as input to the function instead of 
VCF. Please refer to the function documentation (`?cyrius()`) for more 
information about calling CYP2D6.

Cyrius is a Python program and needs certain Python dependencies to run 
successfully. In order to run `cyrius()`, first install the necessary 
dependencies using the `install_cyrius()` function and then activate the environment:

```r
# Install the Cyrius dependencies into a conda env called 'r-ursaPGx'
install_cyrius()

# Activate the 'r-ursaPGx' environment
reticulate::use_condaenv("r-ursaPGx")
```

Now you're ready to use `cyrius()` to call CYP2D6 alleles from R:

```r
# Create a vector of BAM/CRAM file paths
cram <- c("HG00276.final.cram", "HG00436.final.cram", "HG00589.final.cram")

# Optionally name the input files
names(cram) <- c("HG00276", "HG00436", "HG00589")

# Specify the path to the reference fasta file used in CRAM creation
fa <- "GRCh38_full_analysis_set_plus_decoy_hla.fa"

# Call CYP2D6 for each of the samples using Cyrius
result <- cyrius(cram, reference = fa)

result

>DataFrame with 3 rows and 3 columns
>             Sample    Genotype      Filter
>        <character> <character> <character>
>HG00276     HG00276       *4/*5        PASS
>HG00436     HG00436    *2x2/*71        PASS
>HG00589     HG00589      *1/*21        PASS
```
