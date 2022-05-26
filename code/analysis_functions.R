library(data.table)
library(rstan)
library(ggplot2)
library(lubridate)
library(loo)
library(patchwork)


impute_vals = function(trow, missing_vals){
  # Impute missing values in an array by using prior and subsequent values
  #
  # Parameters
  # ----------
  #
  # trow: array with missing values
  # missing_vals: What values are considered missing
  #
  # Returns
  # ------
  # : Imputed array


  impute_row = trow
  ind = !(trow %in% missing_vals)

  if((sum(ind) < length(trow)) & (sum(ind) > 0)){

    tcounts = rle(ind)
    vals = tcounts$values
    lens = tcounts$lengths

    # First is missing, fill with next
    if(vals[1] == FALSE){
      start_false = 1:lens[1]
      next_val = trow[lens[1] + 1]
      impute_row[start_false] = next_val
    }

    # 2. Last is missing, fill with previous
    if(vals[length(vals)] == FALSE){
      end_ind = length(trow) - lens[length(lens)]
      end_false = (end_ind + 1):length(trow)
      impute_row[end_false] = trow[end_ind]
    }

    # 3. Fill in between values

    if(length(vals) >= 3){

      truth_mat = cbind(vals[1:(length(vals) - 2)],
                        vals[2:(length(vals) - 1)],
                        vals[3:(length(vals))])
      which_true = truth_mat[, 1] == TRUE

      if(any(which_true)){

        csum = cumsum(lens)
        ind_mat = cbind(csum[1:(length(csum) - 2)],
                        csum[2:(length(csum) - 1)],
                        csum[3:(length(csum))])
        rep_mat = ind_mat[which_true, ,drop=F]

        for(i in 1:nrow(rep_mat)){
          start_val = trow[rep_mat[i, 1]]
          end_val = trow[rep_mat[i, 2] + 1]
          fill_range = (rep_mat[i, 1] + 1):rep_mat[i, 2]
          impute_row[fill_range] = (end_val + start_val) / 2
        }

      }

    }
  } else {
    impute_row = trow
  }

  return(impute_row)
}

inv_logit = function(x){
  return(1 / (1 + exp(-x)))
}

max_na = function(x){

  if(all(is.na(x))) {
    return(NA)
  }else{
    return(max(x, na.rm=TRUE))
  }
}

flat_scale_unflat = function(mat, missing=-1){
  # Flatten, scale (center, divide by standard deviation) and return matrix
  #
  # Parameters
  # ----------
  # mat : matrix to scale
  # missing: The value to consider missing
  #
  # Returns
  # -------
  # : matrix with values scaled

  dmat = dim(mat)
  vmat = as.vector(mat)
  tmean = mean(vmat[vmat != missing])
  tsd = sd(vmat[vmat != missing])

  # Not that we also "scale" the missing but that is OK
  scaled_flat_mat = (vmat - tmean) / tsd
  mat_scaled = matrix(scaled_flat_mat, nrow=dmat[1], ncol=dmat[2])
  return(mat_scaled)

}

flat_scale_unflat_mean_sd = function(mat, missing=-1){
  # Mean and standard deviation of matrix prior to scaling
  #
  # Parameters
  # ----------
  # mat : matrix to scale
  # missing: The value to consider missing
  #
  # Returns
  # -------
  # : c(matrix mean, matrix standard deviation)

  dmat = dim(mat)
  vmat = as.vector(mat)
  tmean = mean(vmat[vmat != missing])
  tsd = sd(vmat[vmat != missing])
  return(c(tmean, tsd))

}

convert_to_states = function(num_inf, num_sampled, abundance, rho){
  # Convert to observed data to observed states in hidden markov model
  #
  # Missing values are -1
  # Observed states are defined as follows
  # 0 : Bd absent
  # 1: Bd present and prevalence less than rho
  # 2: Bd present and prevalence greater than rho
  # 3: Host population extirpated (no hosts observed)
  #
  # Parameters
  # ----------
  # num_inf : matrix (sites by years) of number of frogs observed as infected
  # num_sampled : matrix (sites by years) of number of frogs swabbed
  # abundance : matrix (sites by years) of number of frogs detected on VES surveys
  # rho : cutoff prevalence between high and low prevalence states
  #
  # Returns
  # -------
  # : matrix
  #   Same dimensions as num_inf, but contains observed states and missing values

  state_matrix = array(NA, dim=dim(num_inf))

  state_matrix[(num_inf == -1) | ((num_inf == 0) & (num_sampled == 0) & (abundance > 0))] = -1
  state_matrix[num_inf == 0 & num_sampled > 0] = 0
  state_matrix[((num_inf / num_sampled) <= rho) & (num_inf != -1) & (num_inf != 0)] = 1
  state_matrix[((num_inf / num_sampled) > rho) & (num_inf != -1) & (num_inf != 0)] = 2
  state_matrix[(num_inf == 0) & (num_sampled == 0) & (abundance == 0)] = 3

  # If threes are not at the end then this is a missing sample
  state_matrix = t(apply(state_matrix, 1, adjust_threes))

  return(state_matrix)


}

adjust_threes = function(x){
  # Account for the fact that extinction is irreversible so observed
  # extinction followed by amphibian presence is really a missing observation
  #
  # Parameters
  # ----------
  # x: vector of observed states (e.g., from a site where each entry is a year)
  #
  # Returns
  # -------
  # : vector
  #   Same as length as x, but threes adjusted to missing as necessary

  non_missing = which(x != -1)
  threes = which(x == 3)

  if(length(threes) == 1){

    if(threes < max(non_missing)){

      x[threes] = -1

    }

  } else if (length(threes) > 1){

    x[threes[1:(length(threes) - 1)]] = -1

    if(max(threes) < max(non_missing)){

      x[max(threes)] = -1

    }
  }
  return(x)
  
}


naive_probs = function(sim_res){
  # Given a matrix with observed states, calculate the naive estimates of
  # transition probabilities that DO NOT account for observation error.
  # 
  # A useful function for checking how much power the model might have to infer 
  # particular transition probabilities
  #
  # Parameters
  # ----------
  # sim_res : matrix of states
  #
  # Returns
  # -------
  # : list
  #   "estimates": Parameter estimates a, phi, omega, e, and zeta
  #   "sample_size": The sample size used to compute each of the transition probs
  #                  naively


  # Loop through each row
  n = nrow(sim_res)
  cols = ncol(sim_res)
  all_stacks = list()
  for(i in 1:n){

    traj = sim_res[i, ]
    first = traj[1:(cols - 1)]
    second = traj[2:(cols)]
    comb = cbind(first, second)
    all_stacks[[i]] = comb
  }

  all_df = as.data.table(do.call(rbind, all_stacks))

  # If any are equal to -1, drop them
  all_df = all_df[first != -1 & second != -1]

  # Arrival
  arrival_df = all_df[first == 0]
  a = 1 - (nrow(arrival_df[second == 0]) / nrow(arrival_df))
  arrival_df2 = arrival_df[second != 0]
  zeta = nrow(arrival_df2[second == 1]) / nrow(arrival_df2)

  # Loss df
  loss_df = all_df[first == 1]
  phi = nrow(loss_df[second == 0]) / nrow(loss_df)
  loss_df2 = loss_df[second != 0]
  omega = (nrow(loss_df2[second == 1]) + nrow(arrival_df2[second == 1])) / (nrow(loss_df2) + nrow(arrival_df2))

  # Extinction df
  ext_df = all_df[first == 2]
  ext = nrow(ext_df[second == 3]) / nrow(ext_df)

  sample_sizes = c(a=nrow(arrival_df), 
                phi=nrow(loss_df),
                omega=nrow(loss_df2) + nrow(arrival_df2),
                e=nrow(ext_df),
                zeta=nrow(arrival_df2))
  params = c(a=a, phi=phi, omega=omega, e=ext, zeta=zeta)
  return(list(estimates=params, sample_size=sample_sizes))

}

naive_probs_with_logistic_regression = function(state_mat, density_mat, density_param, 
                                                temp_mat=NULL, tad_mat=NULL, snow_mat=NULL, 
                                                lake_ids=NULL){
  # Given a matrix with observed states, calculate the naive estimates of
  # transition probabilities that DO NOT account for observation error.
  # Account for the effect of density on the naive transition parameter estimates.
  #
  # Parameters
  # ----------
  # state_mat : matrix of observed states
  # density_mat : matrix of observed densities (same dimension as state_mat)
  # density_param : Either "phi" or "omega". Specifies which 
  #                 transition probability has a density effect.
  # temp_mat : matrix of observed temperatures. Optional
  # snow_mat : Matrix of observed winter severity
  # tad_mat : Matrix of observed tadpole presence/absence. Optional
  # lake_ids : vector of lake_ids, same length as nrow(state_mat)
  #
  # Returns
  # -------
  # list
  #   "estimates": The naive estimates of transition probabilities
  #   "sample_size": The sample size used to compute each transition prob.
  #   "data": data used to estimate density effects on phi or omega


  # Loop through each row
  n = nrow(state_mat)
  cols = ncol(state_mat)
  all_stacks = list()

  for(i in 1:n){

    traj = state_mat[i, ]
    traj_density = density_mat[i, ]
    first = traj[1:(cols - 1)]
    second = traj[2:(cols)]
    density = traj_density[1:(cols - 1)]

    if(!is.null(temp_mat)) {

      traj_temp = temp_mat[i, ]
      temp = traj_temp[1:(cols - 1)]

    } else{

      temp = rep(NA, length(first))

    }

    if(!is.null(tad_mat)) {

      traj_tad = tad_mat[i, ]
      tad = traj_tad[1:(cols - 1)]

    } else{

      tad = rep(NA, length(first))

    }

    if(!is.null(snow_mat)) {

      traj_snow = snow_mat[i, ]
      snow = traj_snow[1:(cols - 1)]

    } else{

      snow = rep(NA, length(first))

    }

    if(!is.null(lake_ids)){
      lid = rep(lake_ids[i], length(first))
    } else{
      lid = rep(i, length(first))
    }

    comb = cbind(first, second, density, temp, tad, snow, lid)


    all_stacks[[i]] = comb

  }

  all_df = as.data.table(do.call(rbind, all_stacks))

  # If any are equal to -1, drop them
  all_df = all_df[first != -1 & second != -1]

  # Arrival
  arrival_df = all_df[first == 0]
  a = 1 - (nrow(arrival_df[second == 0]) / nrow(arrival_df))
  arrival_df2 = arrival_df[second != 0]
  names(a) = "a"
  # zeta = nrow(arrival_df2[second == 1]) / nrow(arrival_df2)

  # Loss df
  loss_df = all_df[first == 1]

  # Fit phi
  y = as.integer(loss_df$second == 0)
  loss_dt = data.frame(y=y, density=loss_df$density, temp=loss_df$temp, tad=loss_df$tad, snow=loss_df$snow, lid=loss_df$lid)

  if(density_param == 'phi'){
    fit_phi = glm(y ~ density, data=loss_dt, family='binomial')
    phi = coef(fit_phi)
    names(phi) = c("beta0", "beta1")
    return_dt = loss_dt
  } else{
    fit_phi = glm(y ~ 1, data=loss_dt, family='binomial')
    phi = 1 / (1 + exp(-coef(fit_phi)))
    names(phi) = "phi"
  }
  
  # Fit omega
  loss_df2 = loss_df[second != 0]
  omega_dt = as.data.table(cbind(c(as.integer(loss_df2$second == 1), as.integer(arrival_df2$second == 1)),
                                 c(loss_df2$density, arrival_df2$density),
                                 c(loss_df2$temp, arrival_df2$temp),
                                 c(loss_df2$tad, arrival_df2$tad),
                                 c(loss_df2$snow, arrival_df2$snow),
                                 c(loss_df2$lid, arrival_df2$lid)))
  # omega_dt = as.data.table(cbind(c(as.integer(as.integer(arrival_df2$second == 1))),
  #                                c(arrival_df2$density)))
  colnames(omega_dt) = c('y', 'density', "temp", "tad", "snow", "lid")

  if(density_param == 'omega'){

    fit_omega = glm(y ~ density, data=omega_dt, family='binomial')
    omega = coef(fit_omega)
    names(omega) = c("beta0", "beta1")
    return_dt = omega_dt

  } else{

    fit_omega = glm(y ~ 1, data=omega_dt, family='binomial')
    omega = 1 / (1 + exp(-coef(fit_omega)))
    names(omega) = c("omega")

  }

  # Extinction df
  ext_df = all_df[first == 2]
  ext = nrow(ext_df[second == 3]) / nrow(ext_df)
  names(ext) = "e"

  sample_sizes = c(a=nrow(arrival_df), 
                phi=nrow(loss_df),
                omega=nrow(omega_dt),
                e=nrow(ext_df))
  params = c(a, phi, omega, ext)
  return(list(estimates=params, sample_size=sample_sizes, data=return_dt))

}

phi_omega_naive_coefs = function(state_mat, abund_mat, temp_mat, tad_mat, snow_mat){
  # Calculate estimates of 
  #
  # Parameters
  # ----------
  # state_mat : matrix of observed states
  # density_mat : matrix of observed densities (same dimension as state_mat)
  # density_param : Either "phi" or "omega". Specifies which 
  #                 transition probability has a density effect.
  # temp_mat : matrix of observed temperatures.
  # tad_mat : Matrix of observed tadpole presence/absence.
  # snow_mat : Matrix of observed winter severity

  #
  # Returns
  # -------
  # : data.table with estimated effects of density, temperature, and tadpole presence
  #   one omega and phi for the NAIVE case

  # fd = fitted_data[[2]]
  # state_mat = fd$obs_mat
  # abund_mat = fd$density
  # temp_mat = fd$temperature
  # tad_mat = fd$tadpole
  # snow_mat = fd$snow
  omega_data = naive_probs_with_logistic_regression(state_mat, abund_mat, "omega", 
                                                    temp_mat=temp_mat, tad_mat=tad_mat, snow=snow_mat)$data
  phi_data = naive_probs_with_logistic_regression(state_mat, abund_mat, "phi", 
                                                  temp_mat=temp_mat, tad_mat=tad_mat, snow=snow_mat)$data
  omega_dt = as.data.table(omega_data)
  phi_dt = as.data.table(phi_data)

  dts = list(omega_dt, phi_dt)
  vars = c('omega', "phi")
  fitted_coefs = list()

  for(i in 1:2){

    dt = dts[[i]]
    var = vars[i]

    # Compute logistic regressions
    map_names = c('density'=paste("beta", var, "density", sep="_"),
                      'temp'=paste("beta", var, "temp", sep="_"),
                      'tad'=paste("beta", var, 'tadpole', sep="_"),
                      'snow'=paste("beta", var, 'snow', sep="_"),
                      'temp:snow'=paste("beta", var, 'tempsnow', sep="_"))
    fit = glm(y ~ density + temp + tad + snow + temp:snow, data=dt, family=binomial(link="probit"))

    ests = as.data.frame(confint(fit))
    colnames(ests) = c("lower", "upper")
    ests$variable = plyr::revalue(rownames(ests), map_names)
    ests$med = coef(fit)
    ests = ests[ests$variable != "(Intercept)", ]

    nm_density = paste("beta", var, "density", sep="_")
    nm_temp = paste("beta", var, "temp", sep="_")
    nm_tadpole = paste("beta", var, "tadpole", sep="_")
    nm_snow = paste("beta", var, "snow", sep="_")
    nm_tempsnow = paste("beta", var, "tempsnow", sep="_")
    pretty_map = c(paste0(var, ", density"),
                   paste0(var, ", temperature"),
                   paste0(var, ", tadpole"),
                   paste0(var, ", snow depth"),
                   paste0(var, ", temp. x snow"))
    names(pretty_map) = c(nm_density, nm_temp, nm_tadpole, nm_snow, nm_tempsnow)
    ests$pretty_variable = plyr::revalue(ests$variable, pretty_map)

    fitted_coefs[[i]] = ests

  }

  fitted_coefs_dt = as.data.table(do.call(rbind, fitted_coefs))
  return(fitted_coefs_dt)

}

format_for_plotting_naive = function(fitted_data, cutoffs){
  # Given formatted data used to fit the model, extract the
  # the phi and omega coefficients for naive estimates of parameters
  # for each level of the cutoff (rho)
  #
  # Parameters
  # ----------
  # fitted_data : list
  #   Each entry in list is a list that contains obs_mat, density, temperature, and tadpole
  #   matrices that correspond to the ith value in cutoffs
  # cutoffs : vector
  #   Cutoff (rho) values that delineate low and high prevalence states.
  #
  # Returns
  # -------
  # : data.table
  #   Estimated coefficients as returned from phi_omega_naive_coefs for each cutoff value


  # Loop through each cutoff and return naive coefficient estimates
  all_coefs = list()
  for(i in 1:length(cutoffs)){

    dat = fitted_data[[i]]
    coefs = phi_omega_naive_coefs(dat$obs_mat, dat$density, dat$temperature, dat$tadpole, dat$snow)
    coefs$cutoff = cutoffs[i]
    all_coefs[[i]] = coefs

  }

  all_coefs_dt = as.data.table(do.call(rbind, all_coefs))
  return(all_coefs_dt)

}

rarify_dat = function(dat, subsample){
  # Rarify data to partially uncorrelate sample size and density
  # 
  # Parameters
  # ----------
  # dat : list
  #   Should contain num_samps (number of individual sampled) 
  #   and num_inf (nubmer of individuals detected as infected)
  # subsample : int 
  #   Maximum number of frogs collected at a location
  #
  # Returns
  # -------
  # : dat
  #  Same as dat, but with 

  ind = dat$num_samps > subsample # Identify sites with greater than subsample samples
  num_inf_sub = dat[ind]$num_inf
  num_samps_sub = dat[ind]$num_samps

  new_inf = array(NA, dim=length(num_inf_sub))
  new_samps = array(NA, dim=length(num_inf_sub))

  for(i in 1:length(num_inf_sub)){
    pa_vect = rep(c(1, 0), c(num_inf_sub[i], num_samps_sub[i] - num_inf_sub[i]))
    samp = sample(pa_vect, subsample, replace=FALSE)
    new_inf[i] = sum(samp == 1)
    new_samps[i] = length(samp)
  }

  dat$num_samps[ind] = new_samps
  dat$num_inf[ind] = new_inf
  return(dat)

}

format_for_plotting = function(param_estimates, cutoffs){
  # Format Bayesian HMM model parameter estimates for plotting
  #
  # Parameters
  # ----------
  # param_estimates : list
  #   List of fitted stan models for each cutoff
  # cutoffs : vector
  #   Cutoffs (rho) used for each model to distinguish high and low prevalence states
  #
  # Returns
  # -------
  # : data.frame
  #   Median, lower (2.5%), and upper (97.5%) credible intervals for parameters across 
  #   model fits

  params = c("beta_omega_density", "beta_phi_density", 
             "beta_omega_temp", "beta_phi_temp",
             "beta_omega_snow", "beta_phi_snow", 
             "beta_omega_tempsnow", "beta_phi_tempsnow",
             "beta_omega_tadpole", "beta_phi_tadpole")

  cis = list()
  for(j in 1:length(cutoffs)){

    ests = data.frame(sapply(params, function(x) quantile(rstan::extract(param_estimates[[j]], pars=params)[[x]], 
                                                             c(0.025, 0.5, 0.975))))
    ests$cutoff = cutoffs[j]
    ests$names = c("lower", "med", "upper")
    cis[[j]] = data.table(ests)
  }

  quant_list = lapply(1:length(cutoffs), function(x) melt(cis[[x]], id.vars=c("cutoff", "names")))
  quant_list_df = do.call(rbind, quant_list)

  # Format for ggplot
  cnames = unique(quant_list_df$names)
  base_dt = quant_list_df[names == "lower", .(cutoff, variable)]

  for(nm in cnames){
    base_dt[[nm]] = quant_list_df[names == nm]$value
  }

  # Make pretty names
  name_map = c("beta_omega_density"="\u03C9, density", 
                "beta_phi_density"="\u03C6, density", 
                "beta_omega_temp"="\u03C9, temperature", 
                "beta_phi_temp"="\u03C6, temperature",
                "beta_omega_snow"="\u03C9, winter severity",
                "beta_phi_snow"="\u03C6, winter severity",
                "beta_omega_tempsnow"="\u03C9, temp. x winter",
                "beta_phi_tempsnow"="\u03C6, temp. x winter",
                "beta_omega_tadpole"="\u03C9, tadpole",
                "beta_phi_tadpole"="\u03C6, tadpole")
  base_dt$pretty_variable = plyr::revalue(base_dt$variable, name_map)

  return(base_dt)
}


prepare_stan_data = function(num_inf_mat, sample_size_mat, abund_missing_mat, cutoff,
                             imputed_cov_mats, a_or_d, num_to_include=3){
  # Format observed Sierra sampling data for HMM analysis
  #
  # Parameters
  # ----------
  # num_inf_mat : matrix
  #   Data.table (site by years) with number of infected swabs in each cell
  # sample_size_mat : matrix
  #   Data.table (site by years) with number of sample hosts in each cell
  # abund_missing_mat : matrix
  #   Data.table (site by years) specify where abundance data is missing
  # cutoff : float
  #   The prevalence value that delineates a high and low prevalence state
  # imputed_cov_mats : list
  #  List of imputed covariates including 'density', 'abund', 'bd', 'maxtemp'
  #
  # Returns
  # -------
  # : list
  #   Data formatted for 

  # 2:17 gives 2004:2019 columns
  # TODO: make this robust to different nubmer of years
  state_mat = convert_to_states(as.matrix(num_inf_mat[, 2:17, with=F]), 
                                as.matrix(sample_size_mat[, 2:17, with=F]),
                                as.matrix(abund_missing_mat[, 2:17, with=F]), 
                                cutoff)
  state_mat = as.data.table(cbind(num_inf_mat$site_id, state_mat))
  colnames(state_mat) = colnames(num_inf_mat)

  # Extract the lake ids that should be fit for this analysis
  anal1_dat = fread("../data/candidate_datasets/analysisI_failedinvasions_excluded.csv")

  # Extract the sites with greather than num_include samples
  num_include = num_to_include
  greaterthan_ind = apply(as.matrix(anal1_dat[, 2:17, with=F]), 1, function(x) sum(x != -1)) >= num_include

  include_sites = anal1_dat$site_id[greaterthan_ind]
  ind = state_mat$site_id %in% include_sites
  n = ncol(state_mat)

  # Print density range
  # tempd = as.matrix(imputed_cov_mats[['density']][lakes %in% include_sites][order(lakes)][, 2:n, with=F])
  # print(c(min(tempd[tempd != 0]), max(tempd)))

  state_mat_red = as.matrix(state_mat[ind][order(site_id)][, 2:n, with=F])
  num_inf_mat_red = as.matrix(num_inf_mat[ind][order(site_id)][, 2:n, with=F])
  sample_size_mat_red = as.matrix(sample_size_mat[ind][order(site_id)][, 2:n, with=F])
  density_mat_red = flat_scale_unflat(log10(as.matrix(imputed_cov_mats[['density']][lakes %in% include_sites][order(lakes)][, 2:n, with=F]) + 1))
  abund_mat_red = flat_scale_unflat(log10(as.matrix(imputed_cov_mats[['abund']][lakes %in% include_sites][order(lakes)][, 2:n, with=F]) + 1))
  intensity_mat_red = flat_scale_unflat(as.matrix(imputed_cov_mats[['bd']][lakes %in% include_sites][order(lakes)][, 2:n, with=F]))
  temp_mat_red = flat_scale_unflat(as.matrix(imputed_cov_mats[['maxtemp']][lakes %in% include_sites][order(lakes)][, 2:n, with=F]))
  snow_mat_red = flat_scale_unflat(as.matrix(imputed_cov_mats[['snowdepth']][lakes %in% include_sites][order(lakes)][, 2:n, with=F]))
  abund_missing_mat_red = as.matrix(abund_missing_mat[ind][order(site_id)][, 2:n, with=F])

  # Check if tadpole abundance is being used 
  max_tad = max(imputed_cov_mats[['tad']][, 2:n, with=F], na.rm=T)
  if(max_tad > 1){
    # Using tadpole abundance
    tadpole_mat_red_unscaled = as.matrix(imputed_cov_mats[['tad']][lakes %in% include_sites][order(lakes)][, 2:n, with=F])
    tadpole_mat_red = flat_scale_unflat(as.matrix(imputed_cov_mats[['tad']][lakes %in% include_sites][order(lakes)][, 2:n, with=F]))
    tadpole_mat_red[tadpole_mat_red_unscaled == -1] = 0 # Set missing values to mean 
  }
  else{
    tadpole_mat_red = as.matrix(imputed_cov_mats[['tad']][lakes %in% include_sites][order(lakes)][, 2:n, with=F])
    tadpole_mat_red[tadpole_mat_red == -1] = 0.5 # Set missing values to mean
  }


  lake_vect = state_mat[ind][order(site_id)]$site_id

  # Specify where intensity data or abundance data are missing
  missing_density = abund_missing_mat_red
  missing_density[abund_missing_mat_red == -1] = 1 # Density/abundance/VES data are missing
  missing_density[abund_missing_mat_red > -1] = 0 # Data are not missing

  # Select appropriate dataset
  if(a_or_d == "abund"){
    adat = abund_mat_red
  } else{
    adat = density_mat_red
  }

  stan_data = list(L=nrow(sample_size_mat_red),
                   Y=ncol(sample_size_mat_red),
                   cutoff=cutoff,
                   cutoff_inv_int=floor(1 / cutoff),
                   obs_mat=state_mat_red,
                   num_inf=num_inf_mat_red,
                   sample_size=sample_size_mat_red,
                   density=adat,
                   intensity=intensity_mat_red,
                   temperature=temp_mat_red,
                   snow=snow_mat_red,
                   missing_density=missing_density,
                   tadpole=tadpole_mat_red,
                   included_sites=lake_vect)

  return(stan_data)
}


max_double = function(x1, x2){
  # Function for returning specific maximum values from two different
  # vectors when either one or the other might have all missing values
  # Either take the max of one vector, the other, or both.
  #
  # Parameters
  # ----------
  # x1 : vector
  # x2 : vector
  # 
  # Returns
  # -------
  # : Summed maximum of two vectors

  if(all(is.na(x1))){
    maxx1 = as.integer(-1)
  } else{
    maxx1 = max(x1, na.rm=T)
  }

  if(all(is.na(x2))){
    maxx2 = as.integer(-1)
  } else{
    maxx2 = max(x2, na.rm=T)
  }

  if(maxx1 == -1){
    return(maxx2)
  } else if(maxx2 == -1){
    return(maxx2)
  } else{
    return(maxx1 + maxx2)
  }
}

tadpole_pa = function(tadpole){
  # Helper function to compute tadpole P/A
  return(ifelse(all(is.na(tadpole)), as.integer(-1), as.integer(any(tadpole > 0, na.rm=T))))
}

tadpole_log_abund = function(tadpole){
  # Helper function to compute tadpole P/A
  return(ifelse(all(is.na(tadpole)), as.numeric(-1), log10(max(tadpole, na.rm=T) + 1)))
}

unpack_looic = function(fit, cutoffs, mod, criteria="ic"){
  # Take in fit and unpack LOOIC values
  # 
  # fit : list of lists where each item is a stan fit
  # cutoffs : the different cutoffs used for the model
  # mod : list of model specifications used to id the model
  #
  # Returns
  # -------
  # : data.table of LOOIC values

  if(criteria == 'diag'){
    res = lapply(fit, function(x) lapply(x, function(y) sum(loo(y)$diagnostics$pareto_k > 0.7)))
  } else {
    res = lapply(fit, function(x) lapply(x, function(y) loo(y)$estimates[3, 1]))
  }

  looic_vals = array(NA, dim=c(length(res[[1]]), length(res)))
  for(i in 1:ncol(looic_vals)){
    looic_vals[, i] = unlist(res[[i]])
  }

  looic_dt = data.table(looic_vals)
  colnames(looic_dt) = names(res)
  looic_dt$rho = c(0.25, 0.333, 0.5) # This should be stored after each fit
  looic_dt$model = paste(paste0(names(unlist(mod)), unlist(mod)), collapse="_")
  return(looic_dt)

}


set_up_response_and_covariate_data = function(params, full_data, temperature_data, snow_data){

  # Function builds data matrices and covariate matrices needed for the HMM analysis.
  # Returns a list of data.frames and covariates after processing
  #
  # Parameters
  # ----------
  # params : list
  #   The boolean choices used to specify the model type
  # full_data : data.frame
  #   The full Bd dataset from Sierra Nevada system
  # temperature : data.frame
  #   Summer temperature data
  # snow_data : data.frame
  #   CDEC snow depth data from April 1
  #
  # Returns
  # -------
  # : list
  #   list of data.frames and covariates after processing

  include_subadults = params$include_subadults # Include subadults into abundance vs density measures
  rarify = params$rarify # Reduce sample size in some lakes to remove correlation between density and sample size
  subsample = params$subsample # Maximum number of samples to randomly rarify
  non_linear_temp = params$non_linear_temp # Include a non-linear effect of temperature in the model
  include_snow_depth = params$include_snow_depth # If false, don't include snow depth data
  log_tadpole_abundance = params$log_tadpole_abundance # Should we log10 tadpole abundance rather than presence absence? 

  # Set the tadpole function
  if(log_tadpole_abundance){
    tad_fxn = tadpole_log_abund
  } else{
    tad_fxn = tadpole_pa
  }

  if(include_subadults){

    red_data = full_data[, .(bd_present=ifelse(all(is.na(bd_load)), as.integer(-1), as.integer(any(bd_load > 0))), 
                            num_inf=ifelse(all(is.na(bd_load)), as.integer(0), sum(bd_load > 0)), 
                            num_samps=ifelse(all(is.na(bd_load)), as.integer(0), length(bd_load)),
                            prev=ifelse(all(is.na(bd_load)), as.numeric(-1), sum(bd_load > 0) / length(bd_load)),
                            tadpole=tad_fxn(tadpole),
                            mean_bd=ifelse(all(is.na(bd_load)), as.numeric(-1), ifelse(any(bd_load > 0), mean(log10(bd_load[bd_load > 0])), -1)), 
                            frog_abundance=max_double(adult, subadult)), by=.(site_id, year)][year >= 2004]


  } else {

    red_data = full_data[, .(bd_present=ifelse(all(is.na(bd_load)), as.integer(-1), as.integer(any(bd_load > 0))), 
                            num_inf=ifelse(all(is.na(bd_load)), as.integer(0), sum(bd_load > 0)), 
                            num_samps=ifelse(all(is.na(bd_load)), as.integer(0), length(bd_load)),
                            prev=ifelse(all(is.na(bd_load)), as.numeric(-1), sum(bd_load > 0) / length(bd_load)),
                            tadpole=tad_fxn(tadpole), 
                            mean_bd=ifelse(all(is.na(bd_load)), as.numeric(-1), ifelse(any(bd_load > 0), mean(log10(bd_load[bd_load > 0])), -1)), 
                            frog_abundance=ifelse(all(is.na(adult)), as.integer(-1), max(adult, na.rm=T))), by=.(site_id, year)][year >= 2004]

  }

  # Rarify data to reduce correlation
  if(rarify) { 
    red_data = rarify_dat(red_data, subsample)
  } else{
    subsample = NULL
  }

  # Adjust frog abundance to be the max between sample size a visual encounter survey
  # The rationale is that there has to be at least the number of frogs sampled in the pond.
  # We set this for now, but infer these missing values later.
  max_abund = apply(red_data[, .(frog_abundance, num_samps)], 1, max_na)
  red_data[, frog_abundance_combined:=(max_abund)]

  # Visualize relationship between num samps and frog_abundance
  # plot(log(red_data$frog_abundance + 1), log(red_data$num_samps + 1))
  # vals = seq(min(log(red_data$num_samps + 1)), max(log(red_data$num_samps + 1)), len=50)
  # lines(vals, vals)

  # Cast into wide form for analysis
  obs_bd_mat = dcast(red_data, site_id ~ year,  value.var="bd_present", fill=-1)
  sample_size_mat = dcast(red_data, site_id ~ year,  value.var="num_samps", fill=-1)
  num_inf_mat = dcast(red_data, site_id ~ year,  value.var="num_inf", fill=-1)
  bd_mat = dcast(red_data, site_id ~ year,  value.var="mean_bd", fill=-1)
  abund_mat = dcast(red_data, site_id ~ year,  value.var="frog_abundance_combined", fill=-1)
  abund_missing_mat = dcast(red_data, site_id ~ year,  value.var="frog_abundance", fill=-1)
  tad_mat = dcast(red_data, site_id ~ year,  value.var="tadpole", fill=-1)

  ################################
  ##### Build Covariate Data #####
  ################################

  # Extract relevant lakes and years
  years = as.integer(colnames(obs_bd_mat)[2:length (colnames(obs_bd_mat))])
  year_cols = as.character(years)
  red_temperature = temperature_data[lake_id %in% obs_bd_mat$site_id][year %in% years]
  red_snowpack = snow_data[lake_id %in% obs_bd_mat$site_id][year %in% years]
  colnames(red_temperature)[1] = "site_id"
  colnames(red_snowpack)[3] = "site_id"

  # Get temperature and snow depth covariates
  mintemp_mat = dcast(red_temperature, site_id ~ year,  value.var="min_summer_temp_C", fill=-1)
  maxtemp_mat = dcast(red_temperature, site_id ~ year,  value.var="max_summer_temp_C", fill=-1)
  snowdepth_mat = dcast(red_snowpack, site_id ~ year,  value.var="SWEmm", fill=-1)

  # Check correlation between snow depth + 1 year and temperature within lakes
  melted_temp = melt(maxtemp_mat, id.var="site_id", variable.name="year", value.name="temperature")
  melted_snow =  melt(snowdepth_mat, id.var="site_id", variable.name="year", value.name="snow")
  tempsnow_dt = merge(melted_snow, melted_temp, by=c("site_id", "year"))
  tempsnow_offset = tempsnow_dt[order(site_id, year)][, .(snow=snow[2:length(snow)], 
                                        temperature=temperature[1:(length(temperature) - 1)],
                                        year=year[1:(length(year) - 1)]), 
                                    by=.(site_id)]
  # plot(tempsnow_offset$snow, tempsnow_offset$temperature)

  # Site specific correlations to check for Simpsons paradox
  tempsnow_offset = tempsnow_offset[complete.cases(tempsnow_offset)]
  corvals = tempsnow_offset[, .(cor_value=cor(temperature, snow)), by=.(site_id)]
  # hist(corvals$cor_value)
  # median(corvals$cor_value)
  # mean(corvals$cor_value)
  # ggplot(tempsnow_offset) + geom_point(aes(x=temperature, y=snow)) + stat_smooth(aes(x=temperature, y=snow, group=site_id), method="lm", se=FALSE)

  # Impute missing covariate data

  # Make a list of the covariate matrices
  cov_mats = list(bd=bd_mat,
                  abund=abund_mat,
                  tad=tad_mat,
                  mintemp=mintemp_mat,
                  maxtemp=maxtemp_mat,
                  snowdepth=snowdepth_mat)

  # Loop through and impute the mean between two observations. See impute_vals() for approach
  imputed_cov_mats = list()
  for(nm in names(cov_mats)){

    cmat = copy(cov_mats[[nm]])
    lakes = cmat$site_id

    tmat = as.matrix(cmat[, (year_cols), with=F])

    # IMPUTE by row. Impute both -1 and NA values
    for(i in 1:nrow(tmat)){
      tmat[i, ] = impute_vals(tmat[i, ], c(-1, NA))
    }

    imputed_cov_mats[[nm]] = as.data.table(cbind(lakes, as.data.table(tmat)))

  }

  # Compute density data
  lake_dat = fread("../data/archival/sites_nov2020.csv")
  # Make this approximate perimeter assuming a circle
  areas = sqrt(lake_dat[id %in% imputed_cov_mats[['abund']]$lakes][order(id)]$area / pi)*pi*2
  host_density = apply(as.matrix(imputed_cov_mats[['abund']][, (year_cols), with=F]), 2, function(x) x / areas)
  host_density = as.data.table(cbind(imputed_cov_mats[['abund']]$lakes, host_density))
  colnames(host_density) = c("lakes", year_cols)
  imputed_cov_mats[['density']] = host_density

  all_return = list(imputed_cov_mats=imputed_cov_mats,
                    obs_bd_mat=obs_bd_mat,
                    sample_size_mat=sample_size_mat,
                    num_inf_mat=num_inf_mat,
                    bd_mat=bd_mat,
                    abund_mat=abund_mat,
                    abund_missing_mat=abund_missing_mat,
                    tad_mat=tad_mat)
  return(all_return)

}