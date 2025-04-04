# Load Required Packages --------------------------------------------------
library(TwoSampleMR)    # Core MR analysis
library(ggplot2)         # Advanced visualization
library(ggpubr)          # Publication-ready plots
library(plotly)          # Interactive plots
library(MRPRESSO)        # Pleiotropy test
library(gridExtra)       # Multi-panel plots
library(RColorBrewer)    # Color schemes
library(kableExtra)      # Results tables
library(metap)           # P-value combination

# Custom Functions --------------------------------------------------------
#' Enhanced MR Analysis with Diagnostic Suite
advanced_mr_analysis <- function(exposure_id, outcome_id, clump_r2 = 0.001, clump_kb = 10000) {
  
  # Step 1: Data Extraction & QC
  exposure_dat <- extract_instruments(exposure_id, 
                                      p1 = 5e-8, 
                                      clump = TRUE,
                                      r2 = clump_r2,
                                      kb = clump_kb)
  
  # Quality Check 1: Instrument Sufficiency
  if(nrow(exposure_dat) < 3) {
    warning(paste("Insufficient instruments for", exposure_id))
    return(NULL)
  }
  
  # Step 2: Outcome Data Extraction
  outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcome_id)
  
  # Step 3: Harmonization with strict QC
  dat <- harmonise_data(exposure_dat, outcome_dat, 
                        action = 2,   # Remove ambiguous/palindromic
                        tolerance = 0.08) # Strand tolerance
  
  # Quality Check 2: Post-harmonization
  if(nrow(dat) < 3) {
    warning(paste("Insufficient SNPs after harmonization:", exposure_id, "->", outcome_id))
    return(NULL)
  }
  
  # Step 4: Main MR Analysis
  methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", 
               "mr_raps", "mr_sign")
  res <- mr(dat, method_list = methods) %>%
    mutate(Direction = paste(exposure_id, "→", outcome_id))
  
  # Step 5: Sensitivity Analyses
  sensitivity <- list(
    heterogeneity = mr_heterogeneity(dat),
    pleiotropy = mr_pleiotropy_test(dat),
    leaveoneout = mr_leaveoneout(dat),
    presso = mr_presso(BetaOutcome = "beta.outcome", 
                       BetaExposure = "beta.exposure",
                       SdOutcome = "se.outcome", 
                       SdExposure = "se.exposure",
                       data = dat,
                       NbDistribution = 1000,
                       SignifThreshold = 0.05)
  )
  
  # Step 6: Instrument Strength
  f_stats <- dat %>%
    mutate(F = (beta.exposure^2)/(se.exposure^2)) %>%
    summarise(mean_F = mean(F),
              min_F = min(F),
              prop_weak = mean(F < 10))
  
  return(list(data = dat,
              results = res,
              sensitivity = sensitivity,
              f_stats = f_stats))
}

# Visualization Toolkit ---------------------------------------------------
#' Create Diagnostic Plots Panel
mr_diagnostic_plots <- function(mr_results, plot_title = "") {
  # Color scheme setup
  qual_col <- brewer.pal(8, "Set2")
  
  # 1. Scatter Plot with Regression Lines
  p_scatter <- ggplot(mr_results$data, aes(x = beta.exposure, y = beta.outcome)) +
    geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome,
                      ymax = beta.outcome + 1.96*se.outcome),
                  color = qual_col[1], alpha = 0.5) +
    geom_errorbarh(aes(xmin = beta.exposure - 1.96*se.exposure,
                       xmax = beta.exposure + 1.96*se.exposure),
                   color = qual_col[2], alpha = 0.5) +
    geom_point(aes(text = paste("SNP:", SNP, "\nGene:", gene)), 
               color = qual_col[3], size = 2.5) +
    geom_abline(intercept = 0, slope = mr_results$results$b[1], 
                color = qual_col[4], linetype = "dashed") +
    labs(title = paste("MR Scatter Plot:", plot_title),
         x = "Genetic Effect on Exposure", 
         y = "Genetic Effect on Outcome") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  # 2. Forest Plot with Multiple Methods
  p_forest <- ggplot(mr_results$results, 
                     aes(x = reorder(method, -b), y = b, 
                         ymin = b - 1.96*se, ymax = b + 1.96*se)) +
    geom_pointrange(fatten = 2, color = qual_col[5]) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    coord_flip() +
    labs(title = "Effect Size Comparison", x = "MR Method", y = "Beta (95% CI)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank())
  
  # 3. Leave-One-Out Analysis Plot
  p_loo <- mr_forest_plot(mr_results$sensitivity$leaveoneout)[[1]] +
    scale_color_manual(values = qual_col[6:7]) +
    labs(title = "Leave-One-Out Sensitivity") +
    theme(legend.position = "none")
  
  # 4. Funnel Plot with Asymmetry Test
  p_funnel <- mr_funnel_plot(mr_results$sensitivity$pleiotropy)[[1]] +
    geom_point(color = qual_col[8]) +
    labs(title = "Funnel Plot of SNP Effects") +
    theme_bw(base_size = 10)
  
  # Combine into panel
  ggarrange(p_scatter, 
            ggarrange(p_forest, p_loo, p_funnel, ncol = 3), 
            nrow = 2,
            heights = c(2, 1)) %>%
    annotate_figure(top = text_grob(plot_title, 
                                    face = "bold", 
                                    size = 14))
}

# Analysis Pipeline -------------------------------------------------------
# GWAS ID Configuration (Update with actual IDs)
disease_ids <- list(
  RA = "ebi-a-GCST004245",
  IBD = "ebi-a-GCST004131",
  PSO = "ebi-a-GCST005956",
  SLE = "ebi-a-GCST003156",
  AD = "ebi-a-GCST005647"
)

# Perform Bidirectional Analysis
full_results <- lapply(names(disease_ids)[1:4], function(disease) {
  list(
    forward = advanced_mr_analysis(disease_ids[[disease]], disease_ids$AD),
    reverse = advanced_mr_analysis(disease_ids$AD, disease_ids[[disease]])
  )
}) %>% setNames(names(disease_ids)[1:4])

# Generate Interactive Report ---------------------------------------------
generate_interactive_report <- function(results) {
  # Create plotly objects
  plotly_plots <- lapply(names(results), function(disease) {
    p <- mr_diagnostic_plots(results[[disease]]$forward, 
                             paste(disease, "→ AD")) 
    ggplotly(p) %>% 
      layout(hoverlabel = list(bgcolor = "white"),
             annotations = list(text = disease,
                                xref = "paper", yref = "paper",
                                x = 0.5, y = 1.1, 
                                showarrow = FALSE))
  })
  
  # Create results table
  res_table <- do.call(rbind, lapply(results, function(x) {
    x$forward$results %>% 
      select(Direction, method, b, se, pval)
  })) %>%
    mutate(`FDR Adjusted` = p.adjust(pval, method = "fdr"),
           across(where(is.numeric), ~round(., 4)))
  
  # Build HTML report
  htmltools::browsable(
    tagList(
      h1("Interactive MR Analysis Report"),
      h2("Diagnostic Plots"),
      lapply(plotly_plots, htmltools::tags$div),
      h2("Statistical Results"),
      kable(res_table, "html") %>%
        kable_styling("striped", full_width = FALSE) %>%
        add_header_above(c(" " = 3, "Primary Results" = 2, "Adjustment" = 1))
    )
  )
}

# Execute Analysis & Generate Outputs --------------------------------------
# 1. Generate Static Reports
pdf("MR_Analysis_Report.pdf", width = 14, height = 10)
lapply(names(full_results), function(disease) {
  mr_diagnostic_plots(full_results[[disease]]$forward, 
                      paste(disease, "→ AD"))
  mr_diagnostic_plots(full_results[[disease]]$reverse, 
                      paste("AD →", disease))
})
dev.off()

# 2. Create Interactive HTML Report
generate_interactive_report(full_results)

# 3. Save Comprehensive Results
saveRDS(full_results, "Full_MR_Results.rds")

# 4. Sensitivity Analysis Summary
sensitivity_report <- do.call(rbind, lapply(full_results, function(res) {
  data.frame(
    Exposure = res$forward$results$exposure[1],
    Outcome = res$forward$results$outcome[1],
    Egger_Pleiotropy = res$forward$sensitivity$pleiotropy$pval,
    Cochran_Q = res$forward$sensitivity$heterogeneity$Q_pval,
    PRESSO_Global = res$forward$sensitivity$presso$`MR-PRESSO results`$`Global Test`$Pvalue
  )
}))
write.csv(sensitivity_report, "Sensitivity_Analysis_Summary.csv")

# 5. Combine p-values using Fisher's method
combined_pvalues <- do.call(rbind, lapply(full_results, function(res) {
  data.frame(
    Direction = paste(res$forward$results$exposure[1], "→", res$forward$results$outcome[1]),
    Fisher_p = sumlog(res$forward$results$pval)$p
  )
}))