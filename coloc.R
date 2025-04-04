# Load Required Packages --------------------------------------------------
library(coloc)          # Co-localization analysis
library(ieugwasr)       # GWAS data access
library(dplyr)          # Data manipulation
library(ggplot2)        # Visualization
library(locuscomparr)   # Locus comparison plots
library(GenomicRanges)  # Genomic range operations

# 1. Data Preparation -----------------------------------------------------
#' Fetch and harmonize GWAS data for two traits
#' @param trait1_id GWAS ID for first trait
#' @param trait2_id GWAS ID for second trait
#' @param chr Chromosome number
#' @param start Start position (bp)
#' @param end End position (bp)
#' @param r2 LD threshold for clumping
prepare_coloc_data <- function(trait1_id, trait2_id, chr, start, end, r2=0.01) {
  
  # 1.1 Extract genomic region for both traits
  trait1_dat <- ieugwasr::associations(trait1_id, chr, start, end) %>%
    clump(r2 = r2) %>%
    mutate(MAF = ifelse(eaf > 0.5, 1-eaf, eaf)) %>%
    dplyr::select(chr = chr, pos = position, snp = rsid,
                  p = p, beta = beta, se = se, eaf, MAF)
  
  trait2_dat <- ieugwasr::associations(trait2_id, chr, start, end) %>%
    clump(r2 = r2) %>%
    mutate(MAF = ifelse(eaf > 0.5, 1-eaf, eaf)) %>%
    dplyr::select(chr = chr, pos = position, snp = rsid, 
                  p = p, beta = beta, se = se, eaf, MAF)
  
  # 1.2 Harmonize datasets
  merged_dat <- inner_join(trait1_dat, trait2_dat, 
                           by = c("chr", "pos", "snp"),
                           suffix = c("_t1", "_t2"))
  
  if(nrow(merged_dat) < 50) {
    warning("Insufficient overlapping SNPs (<50) for reliable coloc analysis")
    return(NULL)
  }
  
  return(list(trait1 = trait1_dat,
              trait2 = trait2_dat,
              merged = merged_dat))
}

# 2. Colocalization Analysis ----------------------------------------------
perform_coloc_analysis <- function(dataset, 
                                   p1 = 1e-4, p2 = 1e-4, 
                                   p12 = 1e-5) {
  
  # 2.1 Prepare coloc format
  coloc_data <- list(
    pvalues = dataset$merged$p_t1,
    beta = dataset$merged$beta_t1,
    varbeta = (dataset$merged$se_t1)^2,
    type = "quant",
    snp = dataset$merged$snp,
    position = dataset$merged$pos,
    MAF = dataset$merged$MAF_t1
  )
  
  # 2.2 Run coloc.abf
  result <- coloc.abf(
    dataset1 = coloc_data,
    pvalues = dataset$merged$p_t2,
    beta = dataset$merged$beta_t2,
    varbeta = (dataset$merged$se_t2)^2,
    type = "quant",
    p1 = p1, 
    p2 = p2, 
    p12 = p12
  )
  
  # 2.3 Calculate Posterior Probabilities
  pp <- result$summary %>%
    t() %>%
    as.data.frame() %>%
    mutate(Probability = c("PP0", "PP1", "PP2", "PP3", "PP4")) %>%
    arrange(desc(.))
  
  return(list(full_results = result,
              posterior_probs = pp))
}

# 3. Visualization Functions ---------------------------------------------
#' Generate regional association plot
plot_regional_association <- function(dataset, title = "") {
  ggplot(dataset$merged, aes(x = pos)) +
    geom_point(aes(y = -log10(p_t1), color = "Trait 1"), alpha = 0.7) +
    geom_point(aes(y = -log10(p_t2), color = "Trait 2"), alpha = 0.7) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(x = "Genomic Position (bp)", 
         y = "-log10(p-value)",
         title = title) +
    theme_minimal() +
    theme(legend.position = "top")
}

#' Plot posterior probabilities
plot_posterior <- function(coloc_result) {
  ggplot(coloc_result$posterior_probs, aes(x = Probability, y = V1)) +
    geom_col(fill = "#009E73", alpha = 0.8) +
    geom_text(aes(label = round(V1, 3)), vjust = -0.5) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(y = "Posterior Probability", 
         title = "Coloc Posterior Probabilities") +
    theme_classic()
}

# 4. Full Analysis Pipeline -----------------------------------------------
# Example parameters
analysis_params <- list(
  trait1_id = "ebi-a-GCST005647",  # Alzheimer's Disease
  trait2_id = "ebi-a-GCST004245",  # Rheumatoid Arthritis
  chr = 6,
  start = 28000000,
  end = 34000000,  # MHC region
  p1 = 1e-4,
  p2 = 1e-4,
  p12 = 1e-5
)

# Run analysis
dataset <- prepare_coloc_data(analysis_params$trait1_id,
                              analysis_params$trait2_id,
                              analysis_params$chr,
                              analysis_params$start,
                              analysis_params$end)

coloc_result <- perform_coloc_analysis(dataset,
                                       p1 = analysis_params$p1,
                                       p2 = analysis_params$p2,
                                       p12 = analysis_params$p12)

# 5. Generate Report ------------------------------------------------------
# 5.1 Regional Plot
p_regional <- plot_regional_association(dataset, "MHC Region Association")

# 5.2 Posterior Probability Plot
p_posterior <- plot_posterior(coloc_result)

# 5.3 Locus Comparison Plot
p_locus <- locuscomparr::locuscompare(
  in_fn1 = dataset$merged %>% 
    mutate(pval = p_t1) %>%
    dplyr::select(rsid = snp, pval),
  in_fn2 = dataset$merged %>% 
    mutate(pval = p_t2) %>%
    dplyr::select(rsid = snp, pval),
  title1 = "Alzheimer's Disease",
  title2 = "Rheumatoid Arthritis"
)

# 6. Save Outputs ---------------------------------------------------------
# Save plots
ggsave("regional_association.pdf", p_regional, width = 10, height = 6)
ggsave("posterior_probabilities.pdf", p_posterior, width = 6, height = 6)
ggsave("locus_comparison.png", p_locus, width = 10, height = 8)

# Save numerical results
write.csv(coloc_result$posterior_probs, "coloc_posterior_probabilities.csv")
saveRDS(coloc_result, "full_coloc_results.rds")

# 7. Sensitivity Analysis -------------------------------------------------
#' Run coloc with different prior probabilities
sensitivity_analysis <- function(dataset) {
  priors_grid <- expand.grid(
    p1 = c(1e-4, 1e-5),
    p2 = c(1e-4, 1e-5),
    p12 = c(1e-5, 1e-6)
  )
  
  sensitivity_results <- apply(priors_grid, 1, function(params) {
    res <- perform_coloc_analysis(dataset,
                                  p1 = params["p1"],
                                  p2 = params["p2"],
                                  p12 = params["p12"])
    c(params, PP4 = res$posterior_probs$V1[5])
  })
  
  return(bind_rows(sensitivity_results))
}

# Run sensitivity analysis
sens_results <- sensitivity_analysis(dataset)
write.csv(sens_results, "coloc_sensitivity_analysis.csv")

### 结果解释指南 -----------------------------------------------------------
# PP4 > 0.8: 强共定位证据
# PP4 0.5-0.8: 中等共定位证据 
# PP4 < 0.5: 无显著共定位

### 高级功能说明 -----------------------------------------------------------
# 1. 多重检验校正
if(length(coloc_result$full_results$results) > 1) {
  coloc_result$full_results$results <- coloc_result$full_results$results %>%
    mutate(FDR = p.adjust(PP.H4.abf, method = "fdr"))
}

# 2. 条件分析
conditional_analysis <- function(dataset, lead_snp) {
  condition_snps <- dataset$merged %>%
    filter(snp != lead_snp) %>%
    pull(snp)
  
  dataset$merged %>%
    filter(snp %in% condition_snps) %>%
    perform_coloc_analysis()
}

# 3. 跨种族验证
cross_population_coloc <- function(eur_dataset, asn_dataset) {
  common_snps <- intersect(eur_dataset$merged$snp, asn_dataset$merged$snp)
  
  combined_data <- list(
    merged = bind_rows(
      eur_dataset$merged %>% filter(snp %in% common_snps) %>% mutate(pop = "EUR"),
      asn_dataset$merged %>% filter(snp %in% common_snps) %>% mutate(pop = "ASN")
    )
  )
  
  perform_coloc_analysis(combined_data)
}