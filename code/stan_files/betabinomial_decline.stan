functions {

  // The function returns the log probability p(X | alpha) of
  // a Dirichlet-Multinomial model. We integrate out the
  // Multinomial probabilities which are not of interest in our inference.
  //
  // Parameters
  // ----------
  // vector x: Observed abundances. Sum to n
  // vector alphas: Dirichlet parameters of the form I*p_i, where I is the
  //                dispersion parameter and p_i are the mean parameters
  //
  // Returns
  // -------
  // real logpdf : The logpdf of the data
  real dirichlet_multinomial(vector x, vector alphas) {
    real n;
    real alpha0;
    real norm_factor;
    real prod_factor;
    real logpdf;

    n = sum(x);
    alpha0 = sum(alphas);

    norm_factor = lgamma(n + 1) + lgamma(alpha0) - lgamma(n + alpha0);
    prod_factor = sum(lgamma(x + alphas) - (lgamma(x + 1.0) + lgamma(alphas)));
    logpdf = norm_factor + prod_factor;

    return logpdf;

  }

} data {

  int N; // Number of data points
  int P; // Number of predictor variables including intercept
  int start_abund[N]; 
  int end_abund[N];
  matrix[N, P] X;

} parameters {

  real beta0; // Intercept
  vector[P] beta;
  real<lower=0> phi; // Dispersal effects

} transformed parameters{

  vector[N] meanp;
  vector[N] alpha1;
  vector[N] alpha2;

  meanp = inv_logit(beta0 + X * beta);
  alpha1 = meanp * phi;
  alpha2 = (1 - meanp) * phi;

} model{

  vector[2] tdecline;
  vector[2] talphas;

  // Weakly regularizing priors
  beta ~ normal(0, 5);
  phi ~ cauchy(0, 2);

  for(i in 1:N){

    tdecline[1] = start_abund[i] - end_abund[i];
    tdecline[2] = end_abund[i];
    talphas[1] = alpha1[i];
    talphas[2] = alpha2[i];

    target += dirichlet_multinomial(tdecline, talphas);
  }
} generated quantities {

  vector[2] tdecline;
  vector[2] talphas;
  real log_lik[N];

  for(i in 1:N){

    tdecline[1] = start_abund[i] - end_abund[i];
    tdecline[2] = end_abund[i];
    talphas[1] = alpha1[i];
    talphas[2] = alpha2[i];

    log_lik[i] = dirichlet_multinomial(tdecline, talphas);
  }


}