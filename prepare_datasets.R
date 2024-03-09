#prepare_datasets.R

# Extract clinical and expression data from SE experiment objs dataset.
# Input: Extract all datasets from "https://github.com/bhklab/PredictioR/tree/main/data/ICB" excluding (Ravi (two cohorts), Zhao, and Thibaudin) datasets
# Output: Clinical data and expression data 

# Load required packages
library(SummarizedExperiment)
library(dplyr)
library(gh)
library(httr)

# Set Result directory 
results_path <- "~/BHK lab/ICB/PredictioR/results/"

# Fetch directory bhk GitHub source (bhklab/PredictioR/tree/main/data/ICB)
contents <- gh::gh("https://api.github.com/repos/bhklab/PredictioR/contents/data/ICB")

# Exclude three datasets: "Ravi", "Zhao", "Thibaudin"
exclusions <- c("Ravi", "Zhao", "Thibaudin")
rda_files <- vapply(contents, function(x) x$name, character(1))
rda_files <- rda_files[grepl("\\.rda$", rda_files) & !grepl(paste(exclusions, collapse = "|"), rda_files)]

# Downloads and loads an RDA file from a given GitHub URL, then deletes the temporary file.
download_and_load_rda <- function(file_name) {
  local_file_path <- tempfile()
  download.file(sprintf("https://raw.githubusercontent.com/bhklab/PredictioR/main/data/ICB/%s", file_name), local_file_path, mode = "wb")
  load(local_file_path)
  unlink(local_file_path)
}

sapply(rda_files, download_and_load_rda)

# Process each dataset
for (file_name in rda_files) {
  load(file_name)
  
  # Store the base name for each dataset
  dat_name <- sub("^ICB_", "", gsub("\\.rda$", "", file_name))
  
  # Extract expression, clinical data
  expr <- assays(dat_icb)[["gene_expression"]]
  clin <- data.frame(colData(dat_icb))
  clin <- clin[, c("patientid", "sex", "age", "cancer_type", "treatment", "response", "survival_time_pfs", "event_occurred_pfs", "survival_time_os", "event_occurred_os")]
  
  # Check if the clin rownames are identical and in the same order as colnames(expr) using all.equal
  if (!isTRUE(all.equal(rownames(clin), colnames(expr)))) {
    stop("The rownames of clin are not identical or in the same order as colnames(expr).")
  }
  
  # Create directory for each dataset 
  dataset_dir_path <- file.path(results_path, dat_name) 
  if (!dir.exists(dataset_dir_path)) {
    dir.create(dataset_dir_path, recursive = TRUE) 
  }
  
  # Save expression and clinical datasets as CSV files with row names included
  write.csv(expr, file = file.path(dataset_dir_path, paste0(dat_name, "_expr.csv")), row.names = TRUE)
  write.csv(clin, file = file.path(dataset_dir_path, paste0(dat_name, "_clin.csv")), row.names = TRUE)
}
