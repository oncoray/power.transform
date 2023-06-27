// Created 2023-06-26, Alex Zwanenburg
data {
	int n;
	int n_dataset;
	int n_transformer;
	int n_learner;
	int n_iteration
	int id_transformer[n];
	int id_dataset[n];
	int id_learner[n];
	int id_iteration[n];
	real y[n];
}
parameters {
	// Level 0 - transformation (top level)
	vector<lower=-1, upper=1>[n_transformer] beta_transformer;

	// Level 1 - interaction of learners and dataset
	matrix<lower=-1, upper=1>[n_dataset, n_learner] beta_learner;

	// Level 2 - interaction of dataset and iteration
	matrix<lower=-1, upper=1>[n_dataset, n_iteration] beta_iteration;
}
model {
  vector[n] y_measured;

  // The idea is to capture the global influence of power transformations.
  // What is the effect of power transformation, compensating for differences in
  // learner, and variability in iteration. We therefore capture as much local
  // interaction as possible: iteration within each dataset and learner within
  // each dataset. y are ranks normalised to [-1, 1] to remove the effect of
  // choice of metric and any direct dependence on the dataset.

  // Level 2 - iteration and dataset
  for (jj in 1:dataset) {
    for (ii in 1:n_iteration) {
      beta_iteration[jj, ii] ~ normal(0.0, 0.5);
    }
  }

	// Level 1 - learners and dataset
	for (jj in 1:n_dataset) {
	  for (ii in 1:n_learner) {
		  beta_learner[jj, ii] ~ normal(0.0, 0.5);
	  }
	}

  // Level 0 - transformation
  for (ii in 1:n_transformer) {
    beta_transformer[ii] ~ normal(0.0, 0.5);
  }

  for (ii in 1:n) {
    y_measured[ii] = beta_iteration[id_dataset[ii], id_iteration[ii]] + beta_learner[id_dataset[ii], id_learner[ii]] + beta_transformer[id_transformer[ii]];
  }

	y ~ normal(y_measured, 0.2);
}
//generated quantities {
	//real coef_factor[n_factor];
	//real coef_ib[n_ib];
	//
	//coef_factor = b_factor;
	//coef_ib     = b_ib;
//}
