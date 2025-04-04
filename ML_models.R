# Load Required Packages --------------------------------------------------
library(tidymodels)    # Core modeling framework
library(workflowsets)  # Managing multiple workflows
library(DALEXtra)      # Model interpretation
library(finetune)      # Advanced tuning methods
library(doParallel)    # Parallel processing
library(vip)           # Variable importance
library(ggpubr)        # Publication-ready plots

# Initialize Parallel Processing ------------------------------------------
# Set up parallel processing to accelerate tuning
cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(cores)
registerDoParallel(cl)

# Data Preparation --------------------------------------------------------
# Set random seed for reproducibility
set.seed(2023)

# Clean column names
colnames(msigma) <- make.names(colnames(msigma))

# Create formula
formula_string <- as.formula(paste("status ~", 
                                   paste(setdiff(colnames(msigma), "status"), 
                                         collapse = " + ")))

# Data Splitting ----------------------------------------------------------
# Stratified split to maintain class distribution
split_data <- initial_split(msigma, prop = 0.7, strata = status)
train_data <- training(split_data)
test_data <- testing(split_data)

# Create 10-fold cross-validation object with 5 repeats
cv_folds <- vfold_cv(train_data, v = 10, repeats = 5, strata = status)

# Preprocessing Recipe ----------------------------------------------------
preprocess_recipe <- recipe(formula_string, data = train_data) %>%
  # Handle numeric predictors
  step_corr(all_numeric(), threshold = 0.8) %>%     # Remove correlated features
  step_YeoJohnson(all_numeric()) %>%                # Power transformation
  step_center(all_numeric()) %>%                    # Center data
  step_scale(all_numeric()) %>%                     # Scale data
  step_zv(all_predictors()) %>%                     # Remove zero-variance
  step_impute_knn(all_numeric()) %>%                # KNN imputation
  step_dummy(all_nominal_predictors())              # Handle categorical

# Model Definitions -------------------------------------------------------
# XGBoost with extended parameter grid
xgb_spec <- boost_tree(
  trees = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  min_n = tune(),
  loss_reduction = tune()
) %>%
  set_engine("xgboost", eval_metric = "logloss") %>%
  set_mode("classification")

# Random Forest with feature importance
rf_spec <- rand_forest(
  mtry = tune(),
  trees = tune(),
  min_n = tune()
) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("classification")

# Neural Network with regularization
nnet_spec <- mlp(
  hidden_units = tune(),
  penalty = tune(),
  epochs = tune()
) %>%
  set_engine("nnet") %>%
  set_mode("classification")

# SVM with radial basis kernel
svm_spec <- svm_rbf(
  cost = tune(),
  rbf_sigma = tune(),
  margin = tune()
) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# Elastic Net Regularized Regression
glmnet_spec <- logistic_reg(
  penalty = tune(),
  mixture = tune()
) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# Model Workflow Set ------------------------------------------------------
model_set <- workflow_set(
  preproc = list(recipe = preprocess_recipe),
  models = list(
    xgboost = xgb_spec,
    random_forest = rf_spec,
    svm = svm_spec,
    glmnet = glmnet_spec,
    nnet = nnet_spec
  )
)

# Tuning Configuration ----------------------------------------------------
tune_ctrl <- control_grid(
  save_pred = TRUE,
  save_workflow = TRUE,
  parallel_over = "resamples",
  verbose = TRUE
)

# Parameter Grids ---------------------------------------------------------
# Custom parameter ranges for each model
param_grids <- list(
  xgboost = parameters(
    trees(range = c(500, 1000)),
    tree_depth(range = c(3, 8)),
    learn_rate(range = c(-2, -1)),  # 0.01 to 0.1
    min_n(range = c(5, 15)),
    loss_reduction(range = c(-5, -1))
  ),
  random_forest = parameters(
    mtry(range = c(1, floor(sqrt(ncol(msigma))))),
    trees(range = c(500, 1000)),
    min_n(range = c(5, 20))
  ),
  svm = parameters(
    cost(range = c(-5, 3)),
    rbf_sigma(range = c(-5, -1)),
    margin(range = c(0.1, 0.5))
  )
)

# Hyperparameter Tuning ---------------------------------------------------
tuning_results <- workflow_map(
  model_set,
  fn = "tune_grid",
  resamples = cv_folds,
  grid = 25,  # Number of parameter combinations
  metrics = metric_set(roc_auc, accuracy, sens, spec),
  verbose = TRUE,
  control = tune_ctrl
)

# Model Evaluation --------------------------------------------------------
# Collect performance metrics
performance_stats <- rank_results(tuning_results, 
                                  rank_metric = "roc_auc", 
                                  select_best = TRUE)

# Visualize Model Comparison
performance_plot <- autoplot(tuning_results) +
  theme_minimal() +
  ggtitle("Model Performance Comparison")

# Select Best Model -------------------------------------------------------
best_model <- tuning_results %>%
  extract_workflow_set_result(id = "recipe_xgboost") %>%
  select_best(metric = "roc_auc")

# Finalize Workflow -------------------------------------------------------
final_wf <- tuning_results %>%
  extract_workflow(id = "recipe_xgboost") %>%
  finalize_workflow(best_model)

# Final Model Training ----------------------------------------------------
final_fit <- last_fit(final_wf, split_data)

# Model Interpretation ----------------------------------------------------
# Variable Importance
vip_plot <- final_fit %>%
  extract_fit_parsnip() %>%
  vip(num_features = 15, geom = "point") +
  theme_minimal()

# DALEX Explanations
explainer <- explain_tidymodels(
  final_fit,
  data = select(train_data, -status),
  y = train_data$status,
  label = "XGBoost"
)

# Partial Dependence Plots
pdp <- model_profile(explainer, variables = "TOP3_features") %>%
  plot(geom = "profiles") +
  theme_minimal()

# Final Evaluation --------------------------------------------------------
# Confusion Matrix
conf_mat <- final_fit %>%
  collect_predictions() %>%
  conf_mat(truth = status, estimate = .pred_class)

# ROC Curve
roc_curve <- final_fit %>%
  collect_predictions() %>%
  roc_curve(truth = status, .pred_Class1) %>%
  autoplot() +
  theme_minimal()

# Save Outputs ------------------------------------------------------------
# Save trained model
saveRDS(final_fit, "final_model.rds")

# Generate Report
report_plots <- ggarrange(performance_plot, vip_plot, roc_curve, 
                          ncol = 2, nrow = 2)
ggsave("model_report.pdf", report_plots, width = 16, height = 12)

# Clean Up ----------------------------------------------------------------
stopCluster(cl)
registerDoSEQ()