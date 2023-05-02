test_that("callableAlleles can be extracted for every sample", {
    chroms <- c("chr10", "chr19", "chr22", "chr1", "chr11", "chr7", "chr13", "chr12")
    d <- here::here("data-raw", "sample_vcf")
    f <- paste0("1kGP_high_coverage_Illumina.", chroms, ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
    files <- file.path(d, f)
    names(files) <- chroms

    # Mapping for chromosome to genes
    df <- data.frame(
        chrom = c("chr19", "chr19", "chr19", "chr10", "chr10", "chr10", 
                  "chr7", "chr7", "chr19", "chr13", "chr12"
        ),
        gene = c("CYP2A13", "CYP2A6", "CYP2B6", "CYP2C19", "CYP2C8", "CYP2C9", 
                 "CYP3A4", "CYP3A5", "CYP4F2", "NUDT15", "SLCO1B1"
        )
    )

    # Read in all gene data and test if each gene has callable alleles
    for (i in 1:nrow(df)) {
        chrom <- df[i, "chrom", drop = TRUE]
        gene <- df[i, "gene", drop = TRUE]
        p <- readPGx(files[chrom], gene = gene)
        print(pgxGene(p))


        p <- determineCallableAlleles(p)
        expect_type(pgxCallableAlleles(p), "character")
        expect_equal(length(pgxCallableAlleles(p)) > 0, TRUE, label = "There exists callable alleles")
    }
})
