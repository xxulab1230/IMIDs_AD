# Tissue-Specific eQTL Analysis Pipeline
# Integrates GTEx data with Mendelian Randomization
# Author: Assistant
# Date: 2023-10-15

# Load Required Packages --------------------------------------------------
library(TwoSampleMR)    # For MR analysis
library(MRInstruments)   # For GWAS data tools
library(ggplot2)          # For visualization
library(dplyr)           # For data manipulation

# Set Global Parameters ---------------------------------------------------
# Customize these parameters for different analyses
INPUT_GENES_FILE <- "/path/to/FUMA_gene_list.txt"  # FUMA gene list path
TISSUES_OF_INTEREST <- c("Whole Blood", "Brain Cortex")  # Target tissues
OUTCOME_ID_AD <- "finn-b-AD_AM"       # AD outcome ID
OUTCOME_ID_IMID <- "ebi-a-GCSTXXXXX" # IMID outcome ID
LD_CLUMP_R2 <- 0.001                 # LD clumping threshold
MAF_THRESHOLD <- 0.3                  # MAF filter
PROXY_RSQ <- 0.8                      # LD proxy threshold

# Data Preparation --------------------------------------------------------
# Load FUMA-annotated significant genes
gene_data <- read.table(INPUT_GENES_FILE, header = TRUE) %>% 
  distinct(symbol, .keep_all = TRUE)

# Load GTEx eQTL dataset
data("gtex_eqtl")

# Tissue-Specific eQTL Processing -----------------------------------------
process_eqtl_data <- function(gene_list, tissues) {
  # Subset GTEx data for target genes and tissues
  eqtl_subset <- gtex_eqtl %>%
    filter(gene_name %in% gene_list$symbol,
           tissue %in% tissues)
  
  # Format and clump eQTL data
  formatted_data <- format_gtex_eqtl(eqtl_subset) %>%
    clump_data(clump_r2 = LD_CLUMP_R2)
  
  return(formatted_data)
}

# Main Analysis Function -------------------------------------------------
run_mr_analysis <- function(exposure_dat, outcome_id) {
  # Extract outcome data with LD proxies
  outcome_dat <- extract_outcome_data(
    snps = exposure_dat$SNP,
    outcomes = outcome_id,
    proxies = 1,
    rsq = PROXY_RSQ,
    maf_threshold = MAF_THRESHOLD
  )
  
  # Harmonize exposure-outcome effects
  harmonized_dat <- harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat
  ) %>% 
    tryCatch(power_prune(.), error = function(e) .)  # Safe pruning
  
  # Perform MR analysis
  mr_results <- mr(harmonized_dat, 
                   method_list = c("mr_ivw", "mr_egger_regression",
                                   "mr_weighted_median", "mr_weighted_mode"))
  
  # Convert to odds ratios for binary outcomes
  if(all(harmonized_dat$outcome.type == "binary")) {
    mr_results <- generate_odds_ratios(mr_results)
  }
  
  return(list(results = mr_results, 
              data = harmonized_dat))
}

# Quality Control Checks -------------------------------------------------
perform_qc_checks <- function(harmonized_data) {
  # 1. Heterogeneity test
  heterogeneity <- mr_heterogeneity(harmonized_data)
  
  # 2. Pleiotropy test
  pleiotropy <- mr_pleiotropy_test(harmonized_data)
  
  # 3. Leave-one-out analysis
  loo_analysis <- mr_leaveoneout(harmonized_data)
  
  return(list(heterogeneity = heterogeneity,
              pleiotropy = pleiotropy,
              loo = loo_analysis))
}

# Visualization Functions ------------------------------------------------
plot_mr_results <- function(mr_results, title_suffix = "") {
  # Create comprehensive visualization panel
  p <- mr_forest_plot(mr_results) +
    theme_minimal() +
    labs(title = paste("MR Effect Estimates", title_suffix),
         x = "Odds Ratio (95% CI)",
         y = "Genetic Instrument") +
    scale_x_continuous(trans = "log10")
  
  return(p)
}

# Main Workflow ----------------------------------------------------------
# Step 1: Process eQTL data for AD analysis
ad_eqtl_data <- process_eqtl_data(gene_data, TISSUES_OF_INTEREST)

# Step 2: Run MR analysis for AD outcome
ad_analysis <- run_mr_analysis(ad_eqtl_data, OUTCOME_ID_AD)

# Step 3: Process eQTL data for IMID analysis
imid_eqtl_data <- process_eqtl_data(gene_data, TISSUES_OF_INTEREST)

# Step 4: Run MR analysis for IMID outcome
imid_analysis <- run_mr_analysis(imid_eqtl_data, OUTCOME_ID_IMID)

# Step 5: Perform quality control checks
ad_qc <- perform_qc_checks(ad_analysis$data)
imid_qc <- perform_qc_checks(imid_analysis$data)

# Step 6: Generate visualizations
ad_plot <- plot_mr_results(ad_analysis$results, "for AD Outcome")
imid_plot <- plot_mr_results(imid_analysis$results, "for IMID Outcome")

# Step 7: Save outputs
write.csv(ad_analysis$results, "AD_MR_Results.csv", row.names = FALSE)
write.csv(imid_analysis$results, "IMID_MR_Results.csv", row.names = FALSE)

# Save QC reports
saveRDS(ad_qc, "AD_QC_Results.rds")
saveRDS(imid_qc, "IMID_QC_Results.rds")

# Save visualizations
ggsave("AD_MR_Plot.pdf", ad_plot, width = 10, height = 8)
ggsave("IMID_MR_Plot.pdf", imid_plot, width = 10, height = 8)

# Optional: Generate interactive HTML report
# Requires additional setup with rmarkdown
# render("MR_Report.Rmd", output_file = "Analysis_Report.html")