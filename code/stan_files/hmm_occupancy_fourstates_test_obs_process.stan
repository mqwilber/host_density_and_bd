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

    matrix create_obs_mat_pa(int L, int Y, matrix obs_mat){
        // Create a simplified matrix that distinguished between 
        // obs Bd present and not

        matrix[L, Y] obs_mat_pa;

        for(l in 1:L){
            for(y in 1:Y){

                if(obs_mat[l, y] == -1){

                    obs_mat_pa[l, y] = -1;

                } else if(obs_mat[l, y] == 0){

                    obs_mat_pa[l, y] = 0;

                } else {

                    obs_mat_pa[l, y] = 1;
                }
            }
        }

        return(obs_mat_pa);
    }

    matrix create_obs_mat_extinct(int L, int Y, matrix obs_mat){
        // Create a simplified matrix that distinguished between 
        // obs Bd present and not

        matrix[L, Y] obs_mat_extinct;

        for(l in 1:L){
            for(y in 1:Y){

                if(obs_mat[l, y] == -1){

                    obs_mat_extinct[l, y] = -1;

                } else if(obs_mat[l, y] == 3){

                    obs_mat_extinct[l, y] = 1;

                } else {

                    obs_mat_extinct[l, y] = 0;
                }
            }
        }

        return(obs_mat_extinct);
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


}
data {

    int L; // Number of lakes
    int Y; // Number of years lakes were sampled
    real cutoff;

    int obs_mat[L, Y]; // Lakes by years observation matrix, -1: Missing, 0: Bd-free, 1: Low prevalence, 2: High prevalence, 3: Extirpated
    int sample_size[L, Y]; // Number of frogs sampled
    int num_inf[L, Y]; // Number of frogs infected

} transformed data {

    matrix[L, Y] obs_mat_matrix;
    matrix[L, Y] obs_mat_pa;
    matrix[L, Y] obs_mat_extinct;
    matrix[Y, 4] long_obs_dat[L]; // Long form data
    matrix[4, 4] Pmat[L, Y];
    matrix[L, Y] sample_size_real;
    matrix[L, Y] num_inf_real;


    obs_mat_pa = create_obs_mat_pa(L, Y, obs_mat_matrix);
    obs_mat_extinct = create_obs_mat_extinct(L, Y, obs_mat_matrix);

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

                // Account for prevalence likelihood
                // if(sample_size[l, y] > 0)
                //     target += binomial_logit_lpmf(num_inf[l, y] | sample_size[l, y], logit_prev[l, ty]);

                // 1: No Bd, 2: Bd low prev, 3: Bd medium prev, 4: Extinction
                // Each row given the prob of observing 1, 2, 3, 4 given 1, 2, 3, or 4. 
                Pmat[l][y][1, 1] = ((1)^long_obs_dat[l][y][1] * 
                                    (0)^long_obs_dat[l][y][2] * 
                                    (0)^long_obs_dat[l][y][3] * 
                                    (0)^long_obs_dat[l][y][4]);

                // Inference with uniform distribution
                // Pmat[l][y][2, 2] = (undetected_prob(sample_size_real[l, y], 0.0, cutoff)^long_obs_dat[l][y][1] *
                //                     summed_beta_cdfs_uniform(1, sample_size[l, y] / 2, sample_size_real[l, y], 0.0, 0.5)^long_obs_dat[l][y][2] * 
                //                     summed_beta_cdfs_uniform((sample_size[l, y] / 2) + 1, sample_size[l, y], sample_size_real[l, y], 0.0, 0.5)^long_obs_dat[l][y][3] *
                //                     (0)^long_obs_dat[l][y][4]);

                Pmat[l][y][2, 2] = (undetected_prob(sample_size_real[l, y], 0.0, cutoff)^long_obs_dat[l][y][1] *
                                    summed_beta_cdfs_sample_dist(1, sample_size[l, y] / 2, sample_size_real[l, y], 0.0, 0.5, num_inf_real[l, y])^long_obs_dat[l][y][2] * 
                                    summed_beta_cdfs_sample_dist((sample_size[l, y] / 2) + 1, sample_size[l, y], sample_size_real[l, y], 0.0, 0.5, num_inf_real[l, y])^long_obs_dat[l][y][3] *
                                    (0)^long_obs_dat[l][y][4]);


                // Inference with uniform distribution
                // Pmat[l][y][3, 3] = (undetected_prob(sample_size_real[l, y], cutoff, 1.0)^long_obs_dat[l][y][1] *
                //               (summed_beta_cdfs_uniform(1, sample_size[l, y] / 2, sample_size_real[l, y], 0.5, 1.0))^long_obs_dat[l][y][2] * 
                //               (summed_beta_cdfs_uniform((sample_size[l, y] / 2) + 1, sample_size[l, y], sample_size_real[l, y], 0.5, 1.0))^long_obs_dat[l][y][3] * 
                //               (0)^long_obs_dat[l][y][4]);

                Pmat[l][y][3, 3] = (undetected_prob(sample_size_real[l, y], cutoff, 1.0)^long_obs_dat[l][y][1] *
                              (summed_beta_cdfs_sample_dist(1, sample_size[l, y] / 2, sample_size_real[l, y], 0.5, 1.0, num_inf_real[l, y]))^long_obs_dat[l][y][2] * 
                              (summed_beta_cdfs_sample_dist((sample_size[l, y] / 2) + 1, sample_size[l, y], sample_size_real[l, y], 0.5, 1.0, num_inf_real[l, y]))^long_obs_dat[l][y][3] * 
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

    // Given arrival probability, transition to low prevalence or high prevalence (zeta)
    real logit_zeta;

    // Given no failed invasion at low prevalence, transition among low and high prevalence states (omega)
    real logit_omega;

    // Extinction from high prevalence (ext)
    real logit_ext;

    // failed invasion parameter phi
    real logit_phi;


} transformed parameters {

    matrix[4, 4] Tmat;

    real<lower=0, upper=1> a;
    real<lower=0, upper=1> zeta;
    real<lower=0, upper=1> omega;
    real<lower=0, upper=1> ext;
    real<lower=0, upper=1> phi;

    a = inv_logit(logit_a);
    zeta = inv_logit(logit_zeta);
    omega = inv_logit(logit_omega);
    ext = inv_logit(logit_ext);
    phi = inv_logit(logit_phi);

    // Building the transition matrix by row
    Tmat[1] = [1 - a, a*(zeta), a*(1 - zeta), 0];
    Tmat[2] = [phi, (omega)*(1 - phi), (1 - omega)*(1 - phi), 0];
    Tmat[3] = [0, 0, 1 - ext, ext];
    Tmat[4] = [0, 0, 0, 1];

} model {

    matrix[Y, 4] stepwise_probs[L];
    row_vector[4] init_occupancy = [1, 0, 0, 0];

    // Weakly regularizing priors
    logit_a ~ normal(0, 5);
    logit_zeta ~ normal(0, 5);
    logit_omega ~ normal(0, 5);
    logit_ext ~ normal(0, 5);
    logit_phi ~ normal(0, 5);

    for(l in 1:L){
        for(y in 1:Y){


            if(y == 1){
                stepwise_probs[l][y] = init_occupancy;
            }
            else{

                stepwise_probs[l][y] = stepwise_probs[l][y - 1] * Tmat * Pmat[l][y];
            }
        }
    }

    // log-likelihood of Bd occupancy vector at time Y
    for(l in 1:L){
        target += log(sum(stepwise_probs[l][Y]));
    }
}
