test_that("Cyrius runs", {
    # Ensure environment is loaded -- this can be modified later
    reticulate::use_condaenv("r-ursaPGx")
    
    cram <- list.files(
        path = here::here("data-raw", "sample_cram"), 
        pattern = "*.cram$",
        full.names = TRUE
        )
    names(cram) <- gsub(".final", "", tools::file_path_sans_ext(basename(cram)))
    fa <- here::here("data-raw", "fasta", "GRCh38_full_analysis_set_plus_decoy_hla.fa")
    result <- cyrius(cram, fa) 
})



