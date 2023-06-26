// Created 2023-06-26, Alex Zwanenburg
data {
	int n;
	int n_dataset;
	int n_transformer;
	int n_learner;
	int id_transformer[n];
	int id_dataset[n];
	int id_learner[n];
	real y[n];
}
parameters {
	// Level 0 - transformation (top level)
	vector<lower=-1, upper=1>[n_transformer] beta_transformer;

	// Level 1 - interaction of learners and dataset
	matrix<lower=-1, upper=1>[n_dataset, n_learner] beta_learner;
}
model {
  vector[n] y_measured;

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
    y_measured[ii] = beta_learner[id_dataset[ii], id_learner[ii]] + beta_transformer[id_transformer[ii]];
  }
	y ~ normal(y_measured, 0.5);
}
//generated quantities {
	//real coef_factor[n_factor];
	//real coef_ib[n_ib];
	//
	//coef_factor = b_factor;
	//coef_ib     = b_ib;
//}
