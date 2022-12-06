#' Check PGx object for overlapping ranges
#' 
#' Return a boolean indicating whether the allele definition ranges are 
#' completely represented in the sample VCF. i.e. all positions in the 
#' reference are also present in the sample VCF. This subset represents the 
#' callable' alleles
#' @param pgx PGx object
#' @param ref VRanges object containing haplotype definitions
#' @return boolean value indicating whether all defined positions are present in
#' the sample PGx object
checkOverlaps <- function(pgx, ref) {
    s <- SummarizedExperiment::rowRanges(pgx)
    ov <- IRanges::subsetByOverlaps(ref, s, type = "equal")
    if (length(ov) == length(ref)) {
        return(TRUE)
    }
    return(FALSE)
}

#' Return all haplotype ranges for the given PGx object
#' 
#' @param pgx PGx object
#' @return list of VRanges
getHaplotypeRanges <- function(pgx) {
    haplotypes <- grep(pgxGene(pgx), pgxHaplotypes(), value = TRUE)
    ranges <- ursaPGx:::grch38_haplotype_grl[haplotypes]
    return(ranges)
}

#' Subset the PGx object to contain only the positions in the reference allele
#' 
#' @param pgx PGx object
#' @param ref VRanges object of callable ranges
#' @return list of PGx objects that contain only the range data from each 
#' callable allele
subsetCallablePGx <- function(pgx, ref) {
    pgx_sub <- IRanges::subsetByOverlaps(pgx, ref, type = "equal")
    return(pgx_sub)
}

#' Call phased diplotypes
#' 
#' @param pgx PGx object
#' @param verbose Display verbose output? Default TRUE.
#' @return DataFrame of phased allele calls for each sample in PGx
callPhasedDiplotypes <- function(pgx) {
    message("Checking for callable alleles...")
    haplotype_ranges <- getHaplotypeRanges(pgx)
    is_ov <- vapply(
        X = haplotype_ranges, 
        FUN = checkOverlaps, 
        pgx = pgx,
        FUN.VALUE = logical(1)
    )
    
    message("Keeping callable ranges...")
    haplotype_ranges <- haplotype_ranges[is_ov]
    
    message("Creating PGx objects from callable ranges...")
    callable <- lapply(haplotype_ranges, subsetCallablePGx, pgx = pgx)
}