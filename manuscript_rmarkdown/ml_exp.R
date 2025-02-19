configurations <- list(
  "config_no_transformation_glm" = list(
    "transformation_method" = "none",
    "learner" = "glm"
  ),
  "config_no_transformation_glm_no_normalisation" = list(
    "transformation_method" = "none",
    "learner" = "glm",
    "normalisation_method" = "none"
  ),
  "config_invariant_robust_glm" = list(
    "transformation_method" = "yeo_johnson_robust",
    "learner" = "glm"
  ),
  "config_invariant_robust_glm_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_robust",
    "learner" = "glm",
    "normalisation_method" = "none"
  ),
  "config_invariant_robust_gof_glm" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "glm"
  ),
  "config_invariant_robust_gof_glm_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "glm",
    "normalisation_method" = "none"
  ),
  "config_conventional_glm" = list(
    "transformation_method" = "yeo_johnson_conventional",
    "learner" = "glm"
  ),
  "config_conventional_glm_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_conventional",
    "learner" = "glm",
    "normalisation_method" = "none"
  ),
  "config_no_transformation_rf" = list(
    "transformation_method" = "none",
    "learner" = "random_forest_ranger"
  ),
  "config_no_transformation_rf_no_normalisation" = list(
    "transformation_method" = "none",
    "learner" = "random_forest_ranger",
    "normalisation_method" = "none"
  ),
  "config_invariant_robust_rf" = list(
    "transformation_method" = "yeo_johnson_robust",
    "learner" = "random_forest_ranger"
  ),
  "config_invariant_robust_rf_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_robust",
    "learner" = "random_forest_ranger",
    "normalisation_method" = "none"
  ),
  "config_invariant_robust_gof_rf" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "random_forest_ranger"
  ),
  "config_invariant_robust_gof_rf_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "random_forest_ranger",
    "normalisation_method" = "none"
  ),
  "config_conventional_rf" = list(
    "transformation_method" = "yeo_johnson_conventional",
    "learner" = "random_forest_ranger"
  ),
  "config_conventional_rf_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_conventional",
    "learner" = "random_forest_ranger",
    "normalisation_method" = "none"
  ),
  "config_no_transformation_xgboost" = list(
    "transformation_method" = "none",
    "learner" = "xgboost_lm"
  ),
  "config_no_transformation_xgboost_no_normalisation" = list(
    "transformation_method" = "none",
    "learner" = "xgboost_lm",
    "normalisation_method" = "none"
  ),
  "config_invariant_robust_xgboost" = list(
    "transformation_method" = "yeo_johnson_robust",
    "learner" = "xgboost_lm"
  ),
  "config_invariant_robust_xgboost_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_robust",
    "learner" = "xgboost_lm",
    "normalisation_method" = "none"
  ),
  "config_invariant_robust_gof_xgboost" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "xgboost_lm"
  ),
  "config_invariant_robust_gof_xgboost_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "xgboost_lm",
    "normalisation_method" = "none"
  ),
  "config_conventional_xgboost" = list(
    "transformation_method" = "yeo_johnson_conventional",
    "learner" = "xgboost_lm"
  ),
  "config_conventional_xgboost_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_conventional",
    "learner" = "xgboost_lm",
    "normalisation_method" = "none"
  )
)

# Experiment without offset ----------------------------------------------------

experiment_dir <- r"(C:\Users\alexz\Documents\GitHub\power.transform\manuscript_rmarkdown\ml_experiment)"

results <- familiar.experiment::run_experiment(
  experiment_dir = experiment_dir,
  experiment_configs = configurations,
  fixed_familiar_parameters = list(
    "fs_method" = "mim"
  ),
  n_repetitions = 5L,
  n_nodes = 18L
)
