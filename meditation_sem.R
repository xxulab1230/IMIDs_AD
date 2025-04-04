# Load required libraries
# Install packages if not already installed
if (!require("mediation")) install.packages("mediation")
if (!require("piecewiseSEM")) install.packages("piecewiseSEM")
if (!require("scales")) install.packages("scales")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tibble")) install.packages("tibble")
if (!require("stringr")) install.packages("stringr")
if (!require("networkD3")) install.packages("networkD3")
if (!require("webshot")) install.packages("webshot")

# Load libraries
library(mediation)      # For mediation analysis
library(piecewiseSEM)   # For structural equation modeling
library(scales)         # For rescaling variables
library(ggplot2)        # For visualization
library(dplyr)          # For data manipulation
library(tibble)         # For creating tibbles
library(stringr)        # For string operations
library(networkD3)      # For network visualizations
library(webshot)        # For saving plots as PDF

# --- Step 1: Data Preprocessing ---
# Assume 'msigma' is your dataset (replace with your actual data loading code)
# Example: msigma <- read.csv("path/to/your/data.csv")

# Rescale continuous variables to a 0-1 range for consistency
msigma$MRPRS_scaled <- rescale(msigma$MRPRS, to = c(0, 1))
msigma$ABETA_scaled <- rescale(msigma$ABETA, to = c(0, 1))

# Recode binary outcome variable 'status' (e.g., 1 = case, 0 = control)
msigma$status <- ifelse(msigma$status == 1, 0, 1)

# Remove missing values to ensure clean data for modeling
msigma_clean <- na.omit(msigma)

# --- Step 2: Basic Mediation Analysis ---
# Define mediation path models
# Path a: Independent variable (MRPRS) → Mediator (ABETA)
model_m <- lm(ABETA ~ MRPRS + AGE + GENDER, data = msigma_clean)

# Path b: Mediator (ABETA) → Outcome (status), controlling for MRPRS
model_y <- glm(status ~ ABETA + MRPRS + AGE + GENDER, 
               family = binomial(link = "logit"), 
               data = msigma_clean)

# Conduct mediation analysis with bootstrapping for robust estimates
set.seed(123)  # Set seed for reproducibility
med_result <- mediate(model_m, model_y, 
                      treat = "MRPRS",        # Treatment/independent variable
                      mediator = "ABETA",     # Mediator variable
                      boot = TRUE,            # Use bootstrapping
                      sims = 5000)            # Number of bootstrap simulations

# Output results
summary(med_result)  # Display mediation effects (ACME, ADE, etc.)
plot(med_result)     # Plot mediation results

# --- Step 3: Complex Mediation with Piecewise SEM ---
# Define a more complex model with multiple predictors
model_mediator <- lm(ABETA ~ STK4 + CDH3 + PABPC1L + TRIM27 + MRPRS + DAB2 + RNH1 + 
                       YWHAB + GOSR1 + S1PR5 + ABT1 + HIST1H1C + ZKSCAN3 + SMPD3 + 
                       HIST1H1T + AP1M2 + ZKSCAN4, data = msigma_clean)

model_outcome <- lm(status ~ ABETA + STK4 + CDH3 + PABPC1L + TRIM27 + MRPRS + DAB2 + RNH1 + 
                      YWHAB + GOSR1 + S1PR5 + ABT1 + HIST1H1C + ZKSCAN3 + SMPD3 + 
                      HIST1H1T + AP1M2 + ZKSCAN4, data = msigma_clean)

# Combine into a piecewise SEM model
mediation_result <- psem(model_mediator, model_outcome, data = msigma_clean)

# Summarize and visualize the SEM model
summary(mediation_result)  # Display path coefficients and model fit
plot(mediation_result)     # Plot the SEM structure

# --- Step 4: Extract and Process Path Coefficients ---
# Get coefficient table from SEM results
coef_table <- summary(mediation_result)$coefficients

# Define expected column names for the coefficient table
correct_colnames <- c("Response", "Predictor", "Estimate", "Std.Error", 
                      "DF", "Crit.Value", "P.Value", "Std.Estimate")

# Subset to relevant columns
coef_table <- coef_table[, colnames(coef_table) %in% correct_colnames]

# --- Step 5: Identify Significant Mediators ---
predictor <- "ABETA"  # Define the mediator of interest

# Filter for significant predictors of the mediator (p < 0.05)
sig_mediators <- coef_table %>%
  filter(Response == predictor & P.Value < 0.05) %>%
  pull(Predictor)

# --- Step 6: Optimize Models Based on Significant Mediators ---
# Create formulas using only significant mediators
optimized_mediator_formula <- as.formula(paste0(predictor, " ~ ", paste(sig_mediators, collapse = " + ")))
optimized_outcome_formula <- as.formula(paste("status ~ ", predictor, " + ", paste(sig_mediators, collapse = " + ")))

# Fit optimized SEM model
final_model <- psem(
  lm(optimized_mediator_formula, data = msigma_clean),
  glm(optimized_outcome_formula, family = binomial, data = msigma_clean)
)

# Summarize and visualize the optimized model
summary(final_model)
plot(final_model)

# --- Step 7: Advanced Mediation Analysis with Multiple Variables ---
# Example with a different treatment variable (STK4) and mediator (ABETA)
med_result_advanced <- mediate(
  model.m = lm(ABETA ~ COLGALT2 + STK4 + SMPD3 + RNH1 + PABPC1L + AGE + GENDER, 
               data = msigma_clean),
  model.y = glm(status ~ ABETA + COLGALT2 + STK4 + SMPD3 + RNH1 + PABPC1L + AGE + GENDER, 
                family = binomial, data = msigma_clean),
  treat = "STK4",
  mediator = "ABETA",
  boot = TRUE,
  sims = 1000
)

# Summarize results
summary(med_result_advanced)

# --- Step 8: Extract Effects for Visualization ---
extract_effects <- function(med_obj) {
  # Define a mapping table for mediation effects
  effect_map <- tribble(
    ~term, ~coef_field, ~ci_field, ~p_field,
    "ACME (control)", "d0", "d0.ci", "d0.p",
    "ACME (treated)", "d1", "d1.ci", "d1.p",
    "ADE (control)", "z0", "z0.ci", "z0.p",
    "ADE (treated)", "z1", "z1.ci", "z1.p",
    "Total Effect", "tau.coef", "tau.ci", "tau.p",
    "Prop. Mediated (control)", "n0", "n0.ci", NA,
    "Prop. Mediated (treated)", "n1", "n1.ci", NA,
    "ACME (average)", "d.avg", "d.avg.ci", "d.avg.p",
    "ADE (average)", "z.avg", "z.avg.ci", "z.avg.p",
    "Prop. Mediated (average)", "n.avg", "n.avg.ci", NA
  )
  
  # Safe extraction function to handle missing fields
  safe_extract <- function(obj, field) {
    if (field %in% names(obj)) return(obj[[field]])
    else if (exists(field, where = obj)) return(get(field, obj))
    else return(NA)
  }
  
  # Extract effects into a data frame
  eff_data <- effect_map %>%
    rowwise() %>%
    mutate(
      estimate = safe_extract(med_obj, coef_field),
      ci = list(safe_extract(med_obj, ci_field)),
      ci_low = if (is.null(ci) || all(is.na(ci))) NA else ci[1],
      ci_high = if (is.null(ci) || all(is.na(ci))) NA else ci[2],
      p_value = safe_extract(med_obj, p_field),
      effect_type = case_when(
        str_detect(term, "ACME") ~ "Mediation",
        str_detect(term, "ADE") ~ "Direct",
        str_detect(term, "Prop") ~ "Proportion",
        TRUE ~ "Total"
      )
    ) %>%
    dplyr::select(-ci) %>%
    ungroup()
  
  return(eff_data)
}

# Extract effects and add significance labels
forest_data <- extract_effects(med_result_advanced) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    estimate_label = if_else(
      effect_type == "Proportion",
      sprintf("%.2f (%.2f–%.2f)", estimate, ci_low, ci_high),
      sprintf("%.1e (%.1e–%.1e)", estimate, ci_low, ci_high)
    )
  )

# --- Step 9: Create Forest Plot with Dual Scaling ---
ggplot(forest_data, aes(x = estimate, y = reorder(term, estimate))) +
  # Error bars for confidence intervals
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high, color = effect_type),
    height = 0.2, linewidth = 1, show.legend = FALSE
  ) +
  # Points for effect estimates
  geom_point(aes(fill = effect_type), size = 3, shape = 21, stroke = 0.5) +
  # Significance markers
  geom_text(aes(label = significance, x = ci_high * 1.2), 
            color = "red", size = 5, hjust = 0) +
  # Percentage labels for proportion mediated
  geom_text(
    aes(label = if_else(effect_type == "Proportion", 
                        sprintf("%.0f%%", estimate * 100), "")),
    x = 1.05, color = "gray30", size = 3.5
  ) +
  # Dual scaling with log scale and percentage
  scale_x_continuous(
    name = "Effect Size (Log Scale)", trans = "log10",
    labels = trans_format("log10", math_format(10^.x)),
    sec.axis = sec_axis(~., name = "Proportion Mediated (%)", 
                        breaks = c(0.01, 0.1, 1), 
                        labels = percent_format(accuracy = 1))
  ) +
  # Custom colors for effect types
  scale_color_manual(values = c("Direct" = "#E69F00", "Mediation" = "#56B4E9", 
                                "Proportion" = "#009E73", "Total" = "#CC79A7")) +
  scale_fill_manual(values = c("Direct" = "#E69F00", "Mediation" = "#56B4E9", 
                               "Proportion" = "#009E73", "Total" = "#CC79A7"), 
                    guide = guide_legend(title = NULL)) +
  # Plot titles and theme
  labs(title = "Mediation Analysis with Dual Scaling", y = "") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "cm"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

# --- Step 10: Save Results ---
# Save the SEM network plot as HTML and convert to PDF
saveNetwork(plot(final_model), paste0("ADNI_", predictor, "_piecewiseSEM_ggplot2.html"))
webshot(paste0("ADNI_", predictor, "_piecewiseSEM_ggplot2.html"), 
        paste0("ADNI_", predictor, "_piecewiseSEM_ggplot2.pdf"), 
        delay = 2, zoom = 2)

# Optionally save the processed dataset
# save(msigma_clean, file = "processed_data.Rdata")