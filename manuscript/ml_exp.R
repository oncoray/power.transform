experiment_dir <- r"(C:\Users\alexz\Documents\GitHub\power.transform\manuscript\ml_experiment)"

configurations <- list(
  "config_no_transformation_glm" = list(
    "transformation_method" = "none",
    "learner" = "glm"
  ),
  "config_new_test_glm" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "glm"
  ),
  "config_old_test_glm" = list(
    "transformation_method" = "yeo_johnson_non_shift",
    "learner" = "glm"
  ),
  "config_no_transformation_rf" = list(
    "transformation_method" = "none",
    "learner" = "random_forest_ranger"
  ),
  "config_new_test_rf" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "random_forest_ranger"
  ),
  "config_old_test_rf" = list(
    "transformation_method" = "yeo_johnson_non_shift",
    "learner" = "random_forest_ranger"
  )
)

results <- familiar.experiment::run_experiment(
  experiment_dir = experiment_dir,
  experiment_configs = configurations,
  fixed_familiar_parameters = list(
    "fs_method" = "mim"
  ),
  n_repetitions = 5L,
  n_nodes = 18L
)

results[, "experiment_parameters" := sub(pattern = "_glm", replacement = "", x = experiment_parameters)]
results[, "experiment_parameters" := sub(pattern = "_rf", replacement = "", x = experiment_parameters)]

results <- data.table::dcast(
  data = results,
  dataset + iteration_id + learner + dataset_split + metric + outcome_type ~ experiment_parameters,
  value.var = "value")

results <- results[dataset_split == "test"]

results[, list("improved" = familiar.experiment::assess_improvement(x = config_new_test, y = config_no_transformation, metric = metric)), by = c("dataset", "iteration_id", "learner")][, list(n = .N), by = c("improved", "learner")][order(improved, learner)]

wilc_data <- results[, list("improved" = familiar.experiment::assess_improvement(x = config_new_test, y = config_no_transformation, metric = metric)), by = c("dataset", "iteration_id", "learner")]
wilcox.test(x=as.numeric(wilc_data$improved) - 2.0, alternative = "greater")
ggplot2::ggplot(data = wilc_data, mapping = ggplot2::aes(x = improved)) + ggplot2::geom_bar()
