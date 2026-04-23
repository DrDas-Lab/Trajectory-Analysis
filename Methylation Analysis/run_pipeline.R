############################################################
## AML MULTIOMICS TRAJECTORY PIPELINE
## Author: Shradha Garnaik
## Description:
## Bulk RNA + DNA-methylation integration in AML
############################################################

rm(list = ls())

# -----------------------------
# Load package manager
# -----------------------------
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

if (file.exists("renv.lock")) {
  message("Checking package synchronization...")
  renv::restore(prompt = FALSE)
}

required_dirs <- c("input", "scripts")
output_dir <- "output"

for (d in required_dirs) {
  if (!dir.exists(d)) {
    stop(sprintf("CRITICAL ERROR -> Required directory '%s' is missing. Pipeline aborted.", d))
  }
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

source("scripts/setup.R")

# -----------------------------
# Run analysis modules
# -----------------------------
modules <- c(
  "scripts/01_preprocess_counts.R",
  "scripts/02_CIBERSORT_pipeline.R",
  "scripts/03_Processing_RNA_Methylation_datasets.R",
  "scripts/04_ModuleScore_trajectory.R",
  "scripts/05_RNA_heatmap.R",
  "scripts/06_Methylation_heatmap.R",
  "scripts/07_Correlation_plots.R"
)

for (script in modules) {
  if (file.exists(script)) {
    
    message("Running: ", script)
    
    start_time <- Sys.time()
    source(script)
    end_time <- Sys.time()
    
    duration <- difftime(end_time, start_time, units = 'secs')
    message(sprintf("* Finished %s in %.2f seconds\n", basename(script), 
                    as.numeric(duration)))
  } else {
    warning("Missing script file: ", script)
  }
}

cat("Pipeline completed successfully\n")
