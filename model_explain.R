# Install and Load Necessary Packages
if (!require("DALEXtra")) install.packages("DALEXtra")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("patchwork")) install.packages("patchwork")

library(DALEXtra)
library(ggplot2)
library(patchwork)

# Create the Explainer Object
explainer <- explain_tidymodels(
  model = final_fit,
  data = train_data,
  y = train_data$status,
  label = best_model
)

# Model Performance
model_perf <- model_performance(explainer)
roc_plot <- plot(model_perf, geom = "roc")
lift_plot <- plot(model_perf, geom = "lift")
performance_plots <- roc_plot + lift_plot + plot_layout(ncol = 2)
print(performance_plots)

# Variable Importance
var_imp <- model_parts(explainer)
var_imp_plot <- plot(var_imp)
print(var_imp_plot)

# Partial Dependence Plots
top_features <- var_imp$variable[1:3]
pdp_plots <- lapply(top_features, function(feature) {
  pdp <- model_profile(explainer, variables = feature)
  plot(pdp)
})
pdp_combined <- wrap_plots(pdp_plots, ncol = 1)
print(pdp_combined)

# Accumulated Local Effects (ALE) Plots
ale_plots <- lapply(top_features, function(feature) {
  ale <- model_profile(explainer, variables = feature, type = "accumulated")
  plot(ale)
})
ale_combined <- wrap_plots(ale_plots, ncol = 1)
print(ale_combined)

# Breakdown Plots for Individual Predictions
new_observation <- test_data[1, ]
bd <- predict_parts(explainer, new_observation = new_observation, type = "break_down")
bd_plot <- plot(bd)
print(bd_plot)

# SHAP Values for Individual Predictions
shap <- predict_parts(explainer, new_observation = new_observation, type = "shap")
shap_plot <- plot(shap)
print(shap_plot)

# Ceteris Paribus Profiles
cp_plots <- lapply(top_features, function(feature) {
  cp <- predict_profile(explainer, new_observation = new_observation, variables = feature)
  plot(cp)
})
cp_combined <- wrap_plots(cp_plots, ncol = 1)
print(cp_combined)

# Model Diagnostics
model_diag <- model_diagnostics(explainer)
residuals_plot <- plot(model_diag, variable = "y", type = "residuals")
pred_obs_plot <- plot(model_diag, variable = "y", type = "prediction")
diag_plots <- residuals_plot + pred_obs_plot + plot_layout(ncol = 2)
print(diag_plots)

# Save and Export Explanations
saveRDS(explainer, file = "model_explainer.rds")
ggsave("variable_importance.png", var_imp_plot, width = 10, height = 8)
ggsave("pdp_plots.png", pdp_combined, width = 10, height = 12)
ggsave("shap_plot.png", shap_plot, width = 10, height = 8)