functions{

    row_vector create_obs_vector(int index){
        // From index, make a 0 and 1 vector 
        // -1, that is missing data and assign as all 1s

        row_vector[4] new_rv = [0, 0, 0, 0];

        if(index != - 1){
            new_rv[index + 1] = 1;
        } else{
            new_rv = [1, 1, 1, 1];
        }

        return(new_rv);

    }


    real undetected_prob(real n, real a, real b){
        // Probability of detecting no Bd given some range of prevalence
        //
        // n : number of frogs sampled
        // a : lower probability bound
        // b : upper probability bound 

        return(-(1 / (n + 1))*((1 - b)^(n + 1) - (1 - a)^(n + 1)));
    }

    real summed_beta_cdfs_uniform(int lower, int upper, real n, real lb, real ub){
        // Observation error with uniform distribution on true prevalence
        //
        // lower : Lower integer on binomial sum
        // upper : Upper integer on binomial sum
        // n : Sample size
        // lb : Lower bound on beta
        // ub : upper bound on beta


        real sumvals;

        sumvals = 0;

        if(lower <= upper){
            for(k in lower:upper){
                sumvals += (1 / (n + 1)) * (beta_cdf(ub, k + 1, n - k + 1) - beta_cdf(lb, k + 1, n - k + 1));
            }
        }   
        return(fabs(sumvals));
    }

    real summed_beta_cdfs_sample_dist(int lower, int upper, real n, real lb, real ub, real y){
        // Observation error using the sampling distribution of true prevalence
        //
        // lower : Lower integer on binomial sum
        // upper : Upper integer on binomial sum
        // n : Sample size
        // lb : Lower bound on beta
        // ub : upper bound on beta
        // y : number of infected hosts observed

        real sumvals;
        real mult_factor;

        sumvals = 0;

        if(lower <= upper){
            for(k in lower:upper){

                mult_factor = exp((lgamma(n + 1) + lgamma(n + 2) + lgamma(k + y + 1) + lgamma(2*n - k - y + 1)) - 
                                  (lgamma(k + 1) + lgamma(n - k + 1) + lgamma(y + 1) + lgamma(n - y + 1) + lgamma(2*n + 2)));
                sumvals += mult_factor * (beta_cdf(ub, k + y + 1, 2*n - k - y + 1) - beta_cdf(lb, k + y + 1, 2*n - k - y + 1));
            }
        }   
        return(fabs(sumvals));
    }

    real single_missing_imputation(real a, real beta, real mu, real sigma){
        // Missing data imputation when one piece of data is missing
        return(Phi((a + beta*mu) / (sqrt(1 + (beta*sigma)^2))));

    }

    real double_missing_imputation(real a, real beta1, real beta2, real mu1, real mu2, real sigma1, real sigma2){
        // Missing data imputation when two pieces of dtaa are missing
        real A;
        real gamma;

        gamma = sqrt(1 + (beta1*sigma1)^2);
        A = (a + beta1*mu1) / gamma;

        return(Phi((A + (beta2 / gamma)*mu2) / (sqrt(1 + ((beta2 / gamma)*sigma2)^2))));

    }


}
data {

    int L; // Number of lakes
    int Y; // Number of years lakes were sampled
    real cutoff;
    int cutoff_inv_int;

    int obs_mat[L, Y]; // Lakes by years observation matrix, -1: Missing, 0: Bd-free, 1: Low prevalence, 2: High prevalence, 3: Extirpated
    int sample_size[L, Y]; // Number of frogs sampled
    int num_inf[L, Y]; // Number of frogs infected
    int missing_density[L, Y];
    real density[L, Y]; // Frog density
    real intensity[L, Y]; // Bd intensity given infection
    real temperature[L, Y]; // Maximum summer temperature data
    real tadpole[L, Y]; // Tadpole: 0/0.5/1

} transformed data {

    matrix[L, Y] obs_mat_matrix;
    matrix[Y, 4] long_obs_dat[L]; // Long form data
    matrix[4, 4] Pmat[L, Y];
    matrix[L, Y] sample_size_real;
    matrix[L, Y] num_inf_real;

    for(l in 1:L){
        for(y in 1:Y){

            long_obs_dat[l][y] = create_obs_vector(obs_mat[l, y]);
            obs_mat_matrix[l, y] = obs_mat[l, y];
            sample_size_real[l, y] = sample_size[l, y];
            num_inf_real[l, y] = num_inf[l, y];

            // Initialize Pmat with zeros
            for(i in 1:4){
                for(j in 1:4){
                    Pmat[l][y][i, j] = 0;
                }
            }


            if(obs_mat[l, y] == -1){

                Pmat[l][y][1, 1] = 1;
                Pmat[l][y][2, 2] = 1;
                Pmat[l][y][3, 3] = 1;
                Pmat[l][y][4, 4] = 1;

            } else {

                // 1: No Bd, 2: Bd low prev, 3: Bd medium prev, 4: Extinction
                // Each row given the prob of observing 1, 2, 3, 4 given 1, 2, 3, or 4. 
                Pmat[l][y][1, 1] = ((1)^long_obs_dat[l][y][1] * 
                                    (0)^long_obs_dat[l][y][2] * 
                                    (0)^long_obs_dat[l][y][3] * 
                                    (0)^long_obs_dat[l][y][4]);


                Pmat[l][y][2, 2] = (undetected_prob(sample_size_real[l, y], 0.0, cutoff)^long_obs_dat[l][y][1] *
                                    summed_beta_cdfs_sample_dist(1, sample_size[l, y] / cutoff_inv_int, sample_size_real[l, y], 0.0, cutoff, num_inf_real[l, y])^long_obs_dat[l][y][2] * 
                                    summed_beta_cdfs_sample_dist((sample_size[l, y] / cutoff_inv_int) + 1, sample_size[l, y], sample_size_real[l, y], 0.0, cutoff, num_inf_real[l, y])^long_obs_dat[l][y][3] *
                                    (0)^long_obs_dat[l][y][4]);


                Pmat[l][y][3, 3] = (undetected_prob(sample_size_real[l, y], cutoff, 1.0)^long_obs_dat[l][y][1] *
                              (summed_beta_cdfs_sample_dist(1, sample_size[l, y] / cutoff_inv_int, sample_size_real[l, y], cutoff, 1.0, num_inf_real[l, y]))^long_obs_dat[l][y][2] * 
                              (summed_beta_cdfs_sample_dist((sample_size[l, y] / cutoff_inv_int) + 1, sample_size[l, y], sample_size_real[l, y], cutoff, 1.0, num_inf_real[l, y]))^long_obs_dat[l][y][3] * 
                              (0)^long_obs_dat[l][y][4]);


                Pmat[l][y][4, 4] = ((0)^long_obs_dat[l][y][1] * 
                                    (0)^long_obs_dat[l][y][2] * 
                                    (0)^long_obs_dat[l][y][3] * 
                                    (1)^long_obs_dat[l][y][4]);
            }

        }
    }

} parameters {

    // Arrival parameter (a)
    real logit_a;

    // Given no failed invasion at low prevalence, transition among low and high prevalence states (omega)
    real beta_omega0;
    real beta_omega_density;
    real beta_omega_temp;
    real beta_omega_tadpole;

    // Extinction from high prevalence (ext)
    real logit_ext;

    // failed invasion parameter phi
    real beta_phi0;
    real beta_phi_density;
    real beta_phi_temp;
    real beta_phi_tadpole;

    real<lower=0> sigma_density;

} transformed parameters {

    matrix[4, 4] Tmat[L, Y - 1]; // Y - 1 one because this takes previous covariates and projects forward
    matrix[L, Y - 1] omega;
    matrix[L, Y - 1] phi;

    real<lower=0, upper=1> ext;
    real<lower=0, upper=1> a;

    ext = inv_logit(logit_ext);
    a = inv_logit(logit_a);

    for(l in 1:L){
        for(y in 1:(Y - 1)){

            if(missing_density[l, y] == 0){

                // Intensity and density are both present
                omega[l, y] = Phi(beta_omega0 + 
                                        beta_omega_density*density[l, y] + 
                                        beta_omega_temp*temperature[l, y] +
                                        beta_omega_tadpole*tadpole[l, y]);
                phi[l, y] = Phi(beta_phi0 + 
                                      beta_phi_density*density[l, y] + 
                                      beta_phi_temp*temperature[l, y] +
                                      beta_phi_tadpole*tadpole[l, y]);

            } else{
                // Just VES survey is missing.

                // We marginalize over tadpole, using our imputed tadpole probability 
                omega[l, y] = tadpole[l, y]*single_missing_imputation(beta_omega0 + beta_omega_temp*temperature[l, y] + beta_omega_tadpole,
                                                                      beta_omega_density,
                                                                      density[l, y],
                                                                      sigma_density) +
                             (1 - tadpole[l, y])*single_missing_imputation(beta_omega0 + beta_omega_temp*temperature[l, y],
                                                                           beta_omega_density,
                                                                           density[l, y],
                                                                           sigma_density);

                phi[l, y] = tadpole[l, y]*single_missing_imputation(beta_phi0 + beta_phi_temp*temperature[l, y] + beta_phi_tadpole,
                                                      beta_phi_density,
                                                      density[l, y],
                                                      sigma_density) + 
                            (1 - tadpole[l, y])*single_missing_imputation(beta_phi0 + beta_phi_temp*temperature[l, y],
                                                                          beta_phi_density,
                                                                          density[l, y],
                                                                          sigma_density);

            } 

            Tmat[l][y][1] = [1 - a, a*(omega[l, y]), a*(1 - omega[l, y]), 0];
            Tmat[l][y][2] = [phi[l, y], (omega[l, y])*(1 - phi[l, y]), (1 - omega[l, y])*(1 - phi[l, y]), 0];
            Tmat[l][y][3] = [0, 0, (1 - ext), ext];
            Tmat[l][y][4] = [0, 0, 0, 1];

        }
    }

} model {

    // matrix[4, 4] Tmat;
    matrix[Y, 4] stepwise_probs[L];
    row_vector[4] init_occupancy = [1, 0, 0, 0];

    beta_omega0 ~ normal(0, 5);
    beta_omega_density ~ normal(0, 3);
    beta_omega_temp ~ normal(0, 3);
    beta_omega_tadpole ~ normal(0, 3);

    beta_phi0 ~ normal(0, 5);
    beta_phi_density ~ normal(0, 3);
    beta_phi_temp ~ normal(0, 3);
    beta_phi_tadpole ~ normal(0, 3);

    logit_a ~ normal(0, 5);
    logit_ext ~ normal(0, 5);

    sigma_density ~ cauchy(0, 2);

    for(l in 1:L){
        for(y in 1:Y){

            if(y == 1){
                stepwise_probs[l][y] = init_occupancy;
            }
            else{

                stepwise_probs[l][y] = stepwise_probs[l][y - 1] * Tmat[l][y - 1] * Pmat[l][y];
            }
        }
    }

    // log-likelihood of Bd occupancy vector at time Y
    for(l in 1:L){
        target += log(sum(stepwise_probs[l][Y]));
    }
} generated quantities {

    matrix[Y, 4] stepwise_probs[L];
    real log_lik[L];
    row_vector[4] init_occupancy = [1, 0, 0, 0];

    for(l in 1:L){
        for(y in 1:Y){

            if(y == 1){
                stepwise_probs[l][y] = init_occupancy;
            }
            else{

                stepwise_probs[l][y] = stepwise_probs[l][y - 1] * Tmat[l][y - 1] * Pmat[l][y];
            }
        }
    }

    // log-likelihood of Bd occupancy vector at time Y
    for(l in 1:L){
        log_lik[l] = log(sum(stepwise_probs[l][Y]));
    }
}
