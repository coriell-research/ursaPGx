.onLoad <- function() {
    message("Checking if miniconda is installed...")
    if (!file.exists(file.path(reticulate::miniconda_path(), "bin/conda"))) {
        reticulate::install_miniconda()
    }
    
    message("Checking if 'r-reticulate' environment exists...")
    envs <- reticulate::conda_list()
    if (!"r-reticulate" %in% envs) {
        reticulate::conda_create("r-reticulate")
        
        # Install necessary packages
        reticulate::conda_install(
            "r-reticulate", 
            packages = c("pysam", "statsmodels", "scipy"),
            channel = c("conda-forge", "anaconda", "bioconda"))
    }
}
