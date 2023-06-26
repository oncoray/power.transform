// Created 2023-06-26, Alex Zwanenburg
// https://stackoverflow.com/questions/29183577/how-to-represent-a-categorical-predictor-rstan
data {
	int<lower=0> n;
	int<lower=0> n_data = 285;
	int<lower=0> n_transformer = 2;
	int<lower=0> n_learner = 1;
	int id_transformer[n];
	int id_learner[n];
	real y[n]
}
parameters {
	// Level 0 - top level
	real<lower=0> beta_transformer[n_transformer];
	real<lower=0> beta_learner[n_learner * n_data];
	real<lower=0> sigma_y;
}
model {
	real y_measured[n]
	sigma_y ~ normal(0.0, 0.5);

	// Learners
	for (jj in 1:(n_data)) {
	  for (ii in 1:(n_learner)) {
		  beta_learner[ii] ~ normal(0.0 , 0.5);
	  }
	}

	// Type of transformation.

	y_measured = y_data + b_transformer

	for(ii in 1:n_meas){
		y_row[ii] = b_factor[factor_id[ii]] + b_ib[ib_id[ii]] + b_ib_fbs[ib_id[ii]] * discr_id[ii];
	}

	y ~ normal(y_row, sigma_y);
}
//generated quantities {
	//real coef_factor[n_factor];
	//real coef_ib[n_ib];
	//
	//coef_factor = b_factor;
	//coef_ib     = b_ib;
//}
