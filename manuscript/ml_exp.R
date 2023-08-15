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
  "config_new_test_glm" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "glm"
  ),
  "config_new_test_glm_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "glm",
    "normalisation_method" = "none"
  ),
  "config_old_test_glm" = list(
    "transformation_method" = "yeo_johnson_non_shift",
    "learner" = "glm"
  ),
  "config_old_test_glm_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_non_shift",
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
  "config_new_test_rf" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "random_forest_ranger"
  ),
  "config_new_test_rf_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "random_forest_ranger",
    "normalisation_method" = "none"
  ),
  "config_old_test_rf" = list(
    "transformation_method" = "yeo_johnson_non_shift",
    "learner" = "random_forest_ranger"
  ),
  "config_old_test_rf_no_normalisation" = list(
    "transformation_method" = "yeo_johnson_non_shift",
    "learner" = "random_forest_ranger",
    "normalisation_method" = "none"
  )
)

# Experiment without offset ----------------------------------------------------

# experiment_dir <- r"(C:\Users\alexz\Documents\GitHub\power.transform\manuscript\ml_experiment)"
#
# results <- familiar.experiment::run_experiment(
#   experiment_dir = experiment_dir,
#   experiment_configs = configurations,
#   fixed_familiar_parameters = list(
#     "fs_method" = "mim"
#   ),
#   n_repetitions = 5L,
#   n_nodes = 18L
# )

# Experiment with offset -------------------------------------------------------

experiment_dir <- r"(C:\Users\alexz\Documents\GitHub\power.transform\manuscript\ml_experiment_offset)"

results <- familiar.experiment::run_experiment(
  experiment_dir = experiment_dir,
  experiment_configs = configurations,
  fixed_familiar_parameters = list(
    "fs_method" = "mim"
  ),
  n_repetitions = 5L,
  numeric_offset = 1E4,
  n_nodes = 18L
)
