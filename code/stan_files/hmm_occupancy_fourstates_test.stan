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

}
data {

    int L; // Number of lakes
    int Y; // Number of years lakes were sampled

    int obs_mat[L, Y]; // Lakes by years observation matrix

} transformed data {

    // Create a long form of the observed data
    matrix[Y, 4] long_obs_dat[L];

    for(l in 1:L){
        for(y in 1:Y){
            long_obs_dat[l][y] = create_obs_vector(obs_mat[l, y]);
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

    // failed invasion parmaeter phi
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
                stepwise_probs[l][y] = stepwise_probs[l][y - 1] * Tmat .* long_obs_dat[l][y];
            }

        }
    }

    // log-likelihood of Bd occupancy vector at time Y
    for(l in 1:L){
        target += log(sum(stepwise_probs[l][Y]));
    }
}
