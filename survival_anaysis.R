# Load required libraries for survival analysis and model explanation
library(mlr3verse)        # Core package for machine learning in R
library(mlr3extralearners) # Additional learners for mlr3
library(mlr3proba)        # Probabilistic prediction extensions
library(mlr3measures)      # Performance measures for mlr3
library(pseudo)            # Utility functions (if needed)
library(survivalmodels)    # Survival-specific models
library(survival)          # Core survival analysis functions
library(survex)            # Survival model explanation package
library(cowplot)           # Plotting utilities

# Assuming 'msigma' (data) and 'meta' (metadata) are pre-loaded in the environment

# Function to preprocess ADNI dataset
preprocess_adni_data <- function(data, meta) {
  # Match rows with metadata to set time based on AGE
  data$time <- as.numeric(meta[match(rownames(data), meta$RID), "AGE"])
  # Normalize time by subtracting the floor of the minimum time divided by 10
  data$time <- data$time - floor(min(data$time, na.rm = TRUE) / 10) * 10
  # Convert status to integer, ensuring character coercion first
  data$status <- as.integer(as.character(data$status))
  # Ensure column names are syntactically valid
  colnames(data) <- make.names(colnames(data), unique = TRUE)
  # Remove rows with NA in critical columns
  data <- na.omit(data[, c("time", "status", setdiff(colnames(data), c("time", "status")))])
  return(data)
}

# Preprocess the data
survival_data <- preprocess_adni_data(msigma, meta)

# Set seed for reproducibility
set.seed(111)

# Create survival task
survival_task <- as_task_surv(
  survival_data,
  time = "time",
  event = "status",
  type = "right"  # Specify right-censored data
)

# Split data into training and testing sets (70% train, 30% test)
data_split <- partition(survival_task, ratio = 0.7)

# Define a list of survival models with specific configurations
survival_models <- list(
  "CoxPH" = lrn("surv.coxph", id = "CoxPH"),  # Cox Proportional Hazards
  "KaplanMeier" = lrn("surv.kaplan", id = "KaplanMeier"),  # Kaplan-Meier estimator
  "Ranger" = lrn("surv.ranger", id = "Ranger", num.trees = 500),  # Random Forest
  "DeepSurv" = lrn("surv.deepsurv", id = "DeepSurv",
                   epochs = 50, batch_size = 32, early_stopping = TRUE,
                   patience = 5),  # DeepSurv with early stopping
  "DeepHit" = lrn("surv.deephit", id = "DeepHit",
                  epochs = 50, batch_size = 32, early_stopping = TRUE,
                  patience = 5),  # DeepHit with early stopping
  "GAMBoost" = lrn("surv.gamboost", id = "GAMBoost"),  # Gradient Boosting GAM
  "XGBoostCox" = lrn("surv.xgboost.cox", id = "XGBoostCox", nrounds = 100),  # XGBoost Cox
  "LogHazard" = lrn("surv.loghaz", id = "LogHazard"),  # Logistic Hazard
  "NelsonAalen" = lrn("surv.nelson", id = "NelsonAalen"),  # Nelson-Aalen estimator
  "CVGLMNet" = lrn("surv.cv_glmnet", id = "CVGLMNet")  # Cross-validated GLMNet
)

# Train models with progress indication
trained_models <- list()
for (model_name in names(survival_models)) {
  cat(sprintf("Training model: %s\n", model_name))
  trained_models[[model_name]] <- survival_models[[model_name]]$train(
    survival_task,
    row_ids = data_split$train
  )
}

# Prepare test data for explainers
test_features <- survival_task$data(
  rows = data_split$test,
  cols = survival_task$feature_names
)
test_surv <- Surv(
  survival_task$data(rows = data_split$test, cols = "time"),
  survival_task$data(rows = data_split$test, cols = "status")
)

# Create explainers for each trained model
model_explainers <- lapply(names(trained_models), function(model_name) {
  cat(sprintf("Creating explainer for: %s\n", model_name))
  explain(
    model = trained_models[[model_name]],
    data = test_features,
    y = test_surv,
    label = model_name,
    verbose = FALSE
  )
})
names(model_explainers) <- names(trained_models)

# Evaluate model performance with error handling
performance_results <- lapply(names(model_explainers), function(explainer_name) {
  cat(sprintf("Evaluating performance for: %s\n", explainer_name))
  tryCatch(
    model_performance(model_explainers[[explainer_name]]),
    error = function(e) {
      message(sprintf("Error evaluating %s: %s", explainer_name, e$message))
      NULL
    }
  )
})
names(performance_results) <- names(model_explainers)

# Filter out NULL performance results
performance_results <- Filter(Negate(is.null), performance_results)

# Define weights for composite score
metric_weights <- c(
  "Integrated C/D AUC" = 0.7,    # Higher is better
  "Integrated Brier score" = -0.3 # Lower is better
)

# Function to extract key metrics safely
extract_key_metrics <- function(perf_obj) {
  result <- perf_obj$result
  list(
    "Integrated C/D AUC" = if (is.null(result$`Integrated C/D AUC`) || is.na(result$`Integrated C/D AUC`)) NA else result$`Integrated C/D AUC`,
    "Integrated Brier score" = if (is.null(result$`Integrated Brier score`) || is.na(result$`Integrated Brier score`)) NA else result$`Integrated Brier score`
  )
}

# Function to compute composite score with NA handling
compute_composite_score <- function(metrics, weights) {
  score <- 0
  used_weight <- 0
  for (metric in names(weights)) {
    val <- metrics[[metric]]
    if (!is.na(val) && is.numeric(val)) {
      score <- score + val * weights[metric]
      used_weight <- used_weight + abs(weights[metric])
    }
  }
  # Scale score if some metrics are missing
  total_weight <- sum(abs(weights))
  if (used_weight > 0) score <- score / used_weight * total_weight else NA
  return(score)
}

# Compile performance metrics into a dataframe
metrics_df <- do.call(rbind, lapply(names(performance_results), function(model_name) {
  metrics <- extract_key_metrics(performance_results[[model_name]])
  data.frame(
    Model = model_name,
    Integrated_CD_AUC = metrics$`Integrated C/D AUC`,
    Integrated_Brier_score = metrics$`Integrated Brier score`,
    Composite_Score = compute_composite_score(metrics, metric_weights),
    stringsAsFactors = FALSE
  )
}))

# Sort by composite score in descending order
metrics_df <- metrics_df[order(metrics_df$Composite_Score, decreasing = TRUE), ]

# Display and save results
print(metrics_df)
setwd("/DATA/data/MRAD/module/model")
dataset_name <- "ADNI-eQTLADgene-MRPRS"
write.csv(metrics_df, file = paste0(dataset_name, "_mlr3model_performance_metrics.csv"), row.names = FALSE)
save(survival_data, performance_results, model_explainers,
     file = paste0(dataset_name, "_mlr3model_performance_results.Rdata"))

# Select top 5 models
top_5_models <- head(metrics_df$Model, 5)

# Function to compute weighted feature importance
compute_weighted_importance <- function(top_models, explainers, imp_weights) {
  # Calculate feature importance for each model
  importance_list <- lapply(top_models, function(model_name) {
    model_parts(explainers[[model_name]])
  })
  names(importance_list) <- top_models
  
  # Extract unique feature names
  feature_names <- unique(unlist(lapply(importance_list, function(x) x$result$variable)))
  imp_df <- data.frame(row.names = feature_names)
  
  # Assign weighted importance
  for (i in seq_along(top_models)) {
    model_name <- top_models[i]
    imp_result <- importance_list[[model_name]]$result
    if ("variable" %in% colnames(imp_result) && "dropout_loss" %in% colnames(imp_result)) {
      imp_values <- setNames(imp_result$dropout_loss, imp_result$variable)
      imp_df[[model_name]] <- imp_values[feature_names] * imp_weights[i]
    } else {
      stop(sprintf("Unexpected importance structure for %s", model_name))
    }
  }
  
  # Compute weighted sum and sort
  imp_df$Weighted_Sum <- rowSums(imp_df, na.rm = TRUE)
  imp_df <- imp_df[order(imp_df$Weighted_Sum, decreasing = TRUE), ]
  return(imp_df)
}

# Calculate and save weighted feature importance
importance_weights <- c(1.0, 0.8, 0.6, 0.4, 0.2)
weighted_importance <- compute_weighted_importance(top_5_models, model_explainers, importance_weights)
write.csv(weighted_importance, file = paste0("mlr3weighted_variable_importance_", dataset_name, ".csv"), row.names = TRUE)

# Optional: Parallelize model training (uncomment to use)
# library(future)
# plan(multisession)  # Use multiple CPU cores
# trained_models <- future_lapply(survival_models, function(model) {
#   model$train(survival_task, row_ids = data_split$train)
# }, future.seed = TRUE)