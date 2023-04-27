experiment_dir <- r"(C:\Users\alexz\Documents\GitHub\power.transform\manuscript\ml_experiment)"

configurations <- list(
  "config_no_transformation_glm" = list(
    "transformation_method" = "none",
    "learner" = "glm"
  ),
  "config_emp_test_glm" = list(
    "transformation_method" = "yeo_johnson_robust",
    "transformation_gof_test_p_value" = 0.01,
    "learner" = "glm"
  )
)

familiar.experiment::run_experiment(
  experiment_dir = experiment_dir,
  experiment_configs = configurations,
  fixed_familiar_parameters = list(
    "fs_method" = "mim"
  ),
  n_repetitions = 1L
)
