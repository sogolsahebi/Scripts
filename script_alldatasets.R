# Libraries
library(MultiAssayExperiment)
library(Hmisc)
library(survival)
library(stringr)

# Set result dir
results_dir <- "~/BHK lab/ICB/PredictioR/results/"

# Define the genes of interest
genes <- c("EGFR", "ERBB2", "CXCL9")

# Function to split study name 
get_study_name <- function(study_icb) {
  parts <- unlist(strsplit(study_icb, "__"))
  paste(parts[1], parts[2], parts[3], sep="__")
}

# Function to perform survival analysis
perform_survival_analysis <- function(expr, clin, study_name, genes) {
  # 1. Association with OS
  res.os <- geneSurvCont(dat.icb = expr,
                         clin = clin,
                         time.censor = 36,
                         n.cutoff = 15,
                         missing.perc = 0.50, 
                         const.int = 0.001,
                         study = study_name,
                         feature = genes,
                         surv.outcome = "OS")
  res.os$FDR <- p.adjust(res.os$Pval, method="BH")
  
  # 2. Association with PFS
  res.pfs <- geneSurvCont(dat.icb = expr,
                          clin = clin,
                          time.censor = 36,
                          n.cutoff = 15,
                          missing.perc = 0.50, 
                          const.int = 0.001,
                          study = study_name,
                          feature = genes,
                          surv.outcome = "PFS")
  res.pfs$FDR <- p.adjust(res.pfs$Pval, method="BH")
  
  # 3. Association with response
  res.logreg <- geneLogReg(dat.icb = expr,
                           clin = clin,
                           n.cutoff = 10,
                           missing.perc = 0.001, 
                           study = study_name,
                           feature = genes,
                           n0.cutoff = 3, 
                           n1.cutoff = 3)
  res.logreg$FDR <- p.adjust(res.logreg$Pval, method="BH")
  
  # Return results
  return(list(OS = res.os, PFS = res.pfs, Response = res.logreg))
}

# Go through each study directory, process clin and expr, and perform analyses
study_results <- lapply(list.dirs(results_dir, full.names = TRUE, recursive = FALSE), function(dir) {
  data.files <- list.files(dir)
  expr.path <- grep("_expr", data.files, value = TRUE)
  clin.path <- grep("_clin", data.files, value = TRUE)
  
  if (length(expr.path) > 0 && length(clin.path) > 0) {
    expr <- read.csv(file.path(dir, expr.path))
    rownames(expr) <- expr$X
    expr <- expr[, -1]
    
    clin <- read.csv(file.path(dir, clin.path))
    
    study_icb <- substr(expr.path, 1, nchar(expr.path) - 10)
    study_name <- get_study_name(study_icb)
    
    # Perform survival analysis
    results <- perform_survival_analysis(expr, clin, study_name, genes)
    
    # Return study name and results
    list(study_name = study_name, results = results)
  }
})

# Print results for each study
for (study_result in study_results) {
  cat("Study:", study_result$study_name, "\n")
  cat("OS Results:", "\n")
  print(study_result$results$OS)
  cat("PFS Results:", "\n")
  print(study_result$results$PFS)
  cat("Response Results:", "\n")
  print(study_result$results$Response)
  cat("\n")
}
