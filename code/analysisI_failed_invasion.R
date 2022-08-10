library(magrittr)
library(data.table)
library(rstan)
library(ggplot2)
library(lubridate)
library(loo)
library(patchwork)
library(bridgesampling)
source("analysis_functions.R")

## This script analyzes whether host density affects failed Bd invasion probability
## It does this through three steps
##
## 1. Prepares abundance, intensity, and covariate data for inclusion in the model
## 2. Fits the HMM model to the data under different conditions
## 3. Compares models
## 4. Plots covariates with uncertainty
## 5. Explores the lack of density effect using model comparison
## 6. Plots example site trajectories

#############################
### Data set preparation ####
#############################

# Load and format data for analysis
full_data = fread("../data/formatted/combined_Bd_data_for_analysis_allVES.csv")
temperature_data = fread("../data/environmental_covariates/summer_lake_temperatures.csv")
snow_data = fread("../data/environmental_covariates/cdec_snowpack_april1.csv")

# NOTE: Results aren't sensitive to adult vs. adult + subadult distinction

# We explore the following combinations in the manuscript
# include_subadults = TRUE; rarify=FALSE; subsample=15; non_linear_temp=FALSE; log_tadpole_abundance=FALSE; include_snow_depth=TRUE
# include_subadults = TRUE; rarify=TRUE; subsample=15; non_linear_temp=FALSE; log_tadpole_abundance=FALSE; include_snow_depth=TRUE
# include_subadults = TRUE; rarify=FALSE; subsample=15; non_linear_temp=FALSE; log_tadpole_abundance=TRUE; include_snow_depth=TRUE
# include_subadults = TRUE; rarify=FALSE; subsample=15; non_linear_temp=TRUE; log_tadpole_abundance=FALSE; include_snow_depth=FALSE
# include_subadults = TRUE; rarify=FALSE; subsample=15; non_linear_temp=FALSE; log_tadpole_abundance=FALSE; include_snow_depth=FALSE

model_params = list(base=list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                              non_linear_temp=FALSE, log_tadpole_abundance=FALSE,
                              include_snow_depth=TRUE),
                    logtad=list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                                non_linear_temp=FALSE, log_tadpole_abundance=TRUE,
                                include_snow_depth=TRUE),
                    nonlintemp=list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                                    non_linear_temp=TRUE, log_tadpole_abundance=FALSE,
                                    include_snow_depth=FALSE),
                    lintemp=list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                                    non_linear_temp=FALSE, log_tadpole_abundance=FALSE,
                                    include_snow_depth=FALSE),
                    rarify_mod=list(rarify=TRUE, subsample=15, include_subadults=TRUE, 
                                    non_linear_temp=FALSE, log_tadpole_abundance=FALSE,
                                    include_snow_depth=TRUE))

# Set all to TRUE to run all analyses
full_fitting = TRUE # Fit model across all rhos
compare_models = TRUE # Generate LOOIC table for comparison among models
example_plots = TRUE # Plot example trajectories
explore_density_effect = TRUE # Examine predictive power of null density effect
make_tadpole_plots = TRUE # Make tadpole plots
plot_all = TRUE # Plot all results together
make_effect_plots = TRUE # Make effect plots for abundance

###################################
#### Fit candidate models      ####
###################################

if(full_fitting){

  for(mp in model_params){

    include_subadults = mp$include_subadults # Include subadults into abundance vs density measures
    rarify = mp$rarify # Reduce sample size in some lakes to remove correlation between density and sample size
    subsample = mp$subsample # Maximum number of samples to randomly rarify
    non_linear_temp = mp$non_linear_temp # Include a non-linear effect of temperature in the model
    include_snow_depth = mp$include_snow_depth # If false, don't include snow depth data
    log_tadpole_abundance = mp$log_tadpole_abundance # Should we log10 tadpole abundance rather than presence absence? 

    # Prepare data
    all_res = set_up_response_and_covariate_data(mp, full_data, temperature_data, snow_data)
    list2env(all_res, .GlobalEnv)

    # Convert to 4 states of HMM: -1: Missing; 0, No bd observed; 1, Low prev observed; 2, High prev observed; 3, extinction
    cutoffs = c(0.25, 0.3333, 0.5)

    # Specify which covariate to use: abundance or density
    abund_or_density = c("abund", "density")
    param_estimates_both = list()
    fitted_data_both = list()

    for(a_or_d in abund_or_density){

      param_estimates = list()
      fitted_data = list()

      for(i in 1:length(cutoffs)){

        cutoff = cutoffs[i]

        stan_data = prepare_stan_data(num_inf_mat, sample_size_mat, abund_missing_mat, cutoff,
                                      imputed_cov_mats, a_or_d)

        if(!non_linear_temp){

          if(include_snow_depth){

            if(log_tadpole_abundance){
              mod2 = stan_model(file="stan_files/hmm_occupancy_fourstates_covariates_imputation_logtadpole.stan")
            } else{
              mod2 = stan_model(file="stan_files/hmm_occupancy_fourstates_covariates_imputation.stan")
            }
            fitfull = sampling(mod2, data=stan_data, iter=2000, warmup=1000, chains=3, cores=3,
                               include=TRUE, pars=c("a", "ext",
                                        'beta_omega0', 'beta_omega_density','beta_omega_temp', 'beta_omega_snow', 'beta_omega_tempsnow', 'beta_omega_tadpole',
                                        'beta_phi0', 'beta_phi_density','beta_phi_temp', 'beta_phi_snow', 'beta_phi_tempsnow', 'beta_phi_tadpole',
                                        'sigma_density', "log_lik"))
          } else{

            # Fit a linear temp effect with no snow depth data
            mod2 = stan_model(file="stan_files/hmm_occupancy_fourstates_covariates_imputation_linear_temp.stan")
            fitfull = sampling(mod2, data=stan_data, iter=2000, warmup=1000, chains=3, cores=3,
                               include=TRUE, pars=c("a", "ext",
                                        'beta_omega0', 'beta_omega_density','beta_omega_temp','beta_omega_tadpole',
                                        'beta_phi0', 'beta_phi_density','beta_phi_temp','beta_phi_tadpole',
                                        'sigma_density', "log_lik"))
          }
        } else{

          include_snow_depth = FALSE # By default, don't include snowdepth here
          cat("Working on cutoff", cutoff, "\n")
          mod2 = stan_model(file="stan_files/hmm_occupancy_fourstates_covariates_imputation_nonlinear_temp.stan")
          fitfull = sampling(mod2, data=stan_data, iter=2000, warmup=1000, chains=2, cores=2,
                             include=TRUE, pars=c("a", "ext",
                                      'beta_omega0', 'beta_omega_density','beta_omega_temp', 'beta_omega_temp2','beta_omega_tadpole',
                                      'beta_phi0', 'beta_phi_density','beta_phi_temp', 'beta_phi_temp2', 'beta_phi_tadpole',
                                      'sigma_density', "log_lik"))

        }

        param_estimates[[i]] = fitfull
        fitted_data[[i]] = stan_data

      }

      param_estimates_both[[a_or_d]] = param_estimates
      fitted_data_both[[a_or_d]] = fitted_data

    }

    # Save and plot values
    saveRDS(fitted_data_both, paste0("../results/fitted_data_rarify_", 
                                        rarify, subsample, "_include_subadults", 
                                        include_subadults, 
                                        "_nonlinear_temp", non_linear_temp,
                                        "_tadpolelog10", log_tadpole_abundance,
                                        "_snowdepth", include_snow_depth, ".rds"))
    saveRDS(param_estimates_both, paste0("../results/model_parameter_estimates_rarify_", 
                                        rarify, subsample, "_include_subadults", 
                                        include_subadults,
                                        "_nonlinear_temp", non_linear_temp,
                                        "_tadpolelog10", log_tadpole_abundance,
                                        "_snowdepth", include_snow_depth, ".rds"))

  } # End model loop

} # End fitting if statement


###################################
#### Model comparison          ####
###################################


if(compare_models){

  model_params = list(base=list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                                non_linear_temp=FALSE, log_tadpole_abundance=FALSE,
                                include_snow_depth=TRUE),
                      logtad=list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                                  non_linear_temp=FALSE, log_tadpole_abundance=TRUE,
                                  include_snow_depth=TRUE),
                      nonlintemp=list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                                      non_linear_temp=TRUE, log_tadpole_abundance=FALSE,
                                      include_snow_depth=FALSE),
                      lintemp=list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                                      non_linear_temp=FALSE, log_tadpole_abundance=FALSE,
                                      include_snow_depth=FALSE))

  all_looic = list()
  all_looic_diag = list()
  for(i in 1:length(model_params)){
    mod = model_params[[i]]
    fit = readRDS(paste0("../results/model_parameter_estimates_rarify_", 
                                      mod$rarify, mod$subsample, "_include_subadults", 
                                      mod$include_subadults, 
                                      "_nonlinear_temp", mod$non_linear_temp,
                                      "_tadpolelog10", mod$log_tadpole_abundance,
                                      "_snowdepth", mod$include_snow_depth, ".rds"))   

    # Unpack LOOIC
    all_looic[[i]] = unpack_looic(fit, cutoffs, mod, criteria="ic")
    all_looic_diag[[i]] = unpack_looic(fit, cutoffs, mod, criteria="diag")

  }

  looic_dt = do.call(rbind, all_looic)
  looic_diag_dt = do.call(rbind, all_looic_diag)
  fwrite(looic_dt, "../results/looic_model_comparison_results.csv")

}


###################################
#### Plot all results together ####
###################################

if(plot_all){

  cutoffs = c(0.25, 0.3333, 0.5) 
  combos = list(list(rarify=TRUE, subsample=15, include_subadults=TRUE, non_linear_temp=FALSE, log_tadpole_abundance=FALSE, include_snow_depth=TRUE), 
                list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, non_linear_temp=FALSE, log_tadpole_abundance=FALSE, include_snow_depth=TRUE))
  # combos = list(list(rarify=FALSE, subsample=NULL, include_subadults=TRUE))
  abund_or_density = list('density', 'abund')
  all_base_a_and_d = list()

  for(j in 1:length(abund_or_density)){

    a_or_d = abund_or_density[[j]]
    all_base = list()
    for(i in 1:length(combos)){

      combo = combos[[i]]
      param_estimates_temp = readRDS(paste0("../results/model_parameter_estimates_rarify_", 
                                        combo$rarify, combo$subsample, "_include_subadults", 
                                        combo$include_subadults,
                                        "_nonlinear_temp", combo$non_linear_temp,
                                        "_tadpolelog10", combo$log_tadpole_abundance, 
                                        "_snowdepth", combo$include_snow_depth, ".rds"))
      fitted_data_temp = readRDS(paste0("../results/fitted_data_rarify_", 
                                        combo$rarify, combo$subsample, "_include_subadults", 
                                        combo$include_subadults, 
                                        "_nonlinear_temp", combo$non_linear_temp,
                                        "_tadpolelog10", combo$log_tadpole_abundance, 
                                        "_snowdepth", combo$include_snow_depth, ".rds"))

      fitted_data = fitted_data_temp[[a_or_d]]
      naive_dt = format_for_plotting_naive(fitted_data, cutoffs)
      naive_dt$bayes = "naive"
      naive_dt$rarify = combo$rarify

      param_estimates = param_estimates_temp[[a_or_d]] 
      base_dt = format_for_plotting(param_estimates, cutoffs)
      base_dt$rarify = combo$rarify
      base_dt$bayes = "bayes"
      all_base[[i]] = rbind(base_dt, naive_dt)

    }

    all_base_dt = do.call(rbind, all_base)
    all_base_dt$rarify_bayes = paste0(all_base_dt$rarify, "_", all_base_dt$bayes)

    all_base_dt$rarify_bayes = factor(all_base_dt$rarify_bayes, levels=c("FALSE_bayes", "TRUE_bayes", "FALSE_naive", "TRUE_naive"))
    all_base_dt$abund_type = ifelse(a_or_d == "density", "density", "abundance")
    all_base_a_and_d[[j]] = all_base_dt
  }

  all_base_a_and_d_dt = do.call(rbind, all_base_a_and_d)



  # Plot rarified and not rarified results
  for(rare in c(FALSE, TRUE)){

    unique(all_base_a_and_d_dt[bayes != "naive" & rarify != rare]$pretty_variable)

    facet_bounds = data.frame(pretty_variable=c("\u03c9, density", 
                                                "\u03c6, density",
                                                "\u03c9, temperature",
                                                "\u03c6, temperature",
                                                "\u03c9, winter severity",
                                                "\u03c6, winter severity",
                                                "\u03c9, temp. x winter", 
                                                "\u03c6, temp. x winter",
                                                "\u03c9, tadpole",
                                                "\u03c6, tadpole" ),
                             ymin=c(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -6.3, -6.3),
                             ymax=c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1))
    ff <- with(facet_bounds,
           data.frame(med=c(ymin,ymax),
                      pretty_variable=c(pretty_variable,pretty_variable)))

    tplot = ggplot(all_base_a_and_d_dt[bayes != "naive" & rarify != rare]) + geom_point(aes(x=as.factor(round(cutoff, 2)), y=med, color=abund_type, shape=abund_type), position=position_dodge(width=0.25)) +
                      geom_errorbar(aes(x=as.factor(round(cutoff, 2)), ymin=lower, ymax=upper, color=abund_type), position=position_dodge(width=0.25), width=0.01) +
                      geom_hline(aes(yintercept=0), linetype="dashed") + 
                      # geom_point(data=ff, aes(x=NA, y=med)) +
                      scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3')) +
                      facet_wrap(~pretty_variable, scales='free_y', nrow=5) + theme_classic() + ylab("Effect size") + 
                      xlab("\u03C1 \n (cutoff between low prev. and high prev. states)") + 
                      guides(color=guide_legend(title=NULL), shape=guide_legend(title=NULL))

    ggsave(paste0("../results/effect_sizes_across_rhos_rarify", !rare, ".jpg"), width=6, height=8)

  }

}


#######################################
##### Make effect plots ###############
#######################################

# Plot effect plots
if(make_effect_plots){

  
  cutoffs = c(0.25, 0.3333, 0.5)
  combo = list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, non_linear_temp=FALSE, log_tadpole_abundance=FALSE, include_snow_depth=TRUE)
  a_or_d = "abund"

  # Extract mean_sd
  all_res = set_up_response_and_covariate_data(combo, full_data, temperature_data, snow_data)
  anal1_dat = fread("../data/candidate_datasets/analysisI_failedinvasions_excluded.csv")
  num_include = 3
  greaterthan_ind = apply(as.matrix(anal1_dat[, 2:17, with=F]), 1, function(x) sum(x != -1)) >= num_include
  include_sites = anal1_dat$site_id[greaterthan_ind]
  imputed_cov_mats = all_res$imputed_cov_mats

  mean_sd = flat_scale_unflat_mean_sd(log10(as.matrix(imputed_cov_mats[[a_or_d]][lakes %in% include_sites][order(lakes)][, 2:17, with=F]) + 1))
  mean_sd_temp = flat_scale_unflat_mean_sd(as.matrix(imputed_cov_mats[['maxtemp']][lakes %in% include_sites][order(lakes)][, 2:17, with=F]))

  # Extract model fits and data used to fit
  param_estimates_temp = readRDS(paste0("../results/model_parameter_estimates_rarify_", 
                                    combo$rarify, combo$subsample, "_include_subadults", 
                                    combo$include_subadults,
                                    "_nonlinear_temp", combo$non_linear_temp,
                                    "_tadpolelog10", combo$log_tadpole_abundance, 
                                    "_snowdepth", combo$include_snow_depth, ".rds"))
  fitted_data_temp = readRDS(paste0("../results/fitted_data_rarify_", 
                                    combo$rarify, combo$subsample, "_include_subadults", 
                                    combo$include_subadults, 
                                    "_nonlinear_temp", combo$non_linear_temp,
                                    "_tadpolelog10", combo$log_tadpole_abundance, 
                                    "_snowdepth", combo$include_snow_depth, ".rds"))


  pest = param_estimates_temp[[a_or_d]]
  cutoff_ind = 1 # Use data from 0.25 cutoff between LP and HP
  psummary = summary(pest[[cutoff_ind]])$summary

  trans_types = c("omega", "phi")

  effect_plots = list()
  effect_plots_temp = list()
  for(t in trans_types){

    # Param plot
    params_nms = row.names(psummary)[row.names(psummary) %like% t]
    params = list()
    for(p in params_nms){
      params[[p]] = rstan::extract(pest[[cutoff_ind]], par=p)[[p]]
    }

    # Now I just need to translate standardized densty back to true density.
    density_vals = seq(-2, 3, len=20)
    temp_vals = seq(-2, 2, len=20)

    dens_effects = sapply(density_vals, function(d) {pnorm(params[[paste0('beta_', t, '0')]] + 
                                                           params[[paste0('beta_', t, '_density')]]*d + 
                                                           params[[paste0('beta_', t, '_tadpole')]] +
                                                           params[[paste0('beta_', t, '_temp')]]*0)})

    temp_effects = sapply(temp_vals, function(d) {pnorm(params[[paste0('beta_', t, '0')]] + 
                                                           params[[paste0('beta_', t, '_density')]]*0 + 
                                                           params[[paste0('beta_', t, '_tadpole')]] +
                                                           params[[paste0('beta_', t, '_temp')]]*d)})

    effects_ul = data.table(t(apply(dens_effects, 2, quantile, c(0.025, 0.5, 0.975))))
    effects_temp_ul = data.table(t(apply(temp_effects, 2, quantile, c(0.025, 0.5, 0.975))))

    colnames(effects_ul) = c("lower", "med", "upper")
    effects_ul$density = density_vals
    effects_ul$param = t
    effect_plots[[t]] = effects_ul

    colnames(effects_temp_ul) = c("lower", "med", "upper")
    effects_temp_ul$temperature = temp_vals
    effects_temp_ul$param = t
    effect_plots_temp[[t]] = effects_temp_ul


  }

  effect_plots_dt = do.call(rbind, effect_plots) 
  effect_plots_dt$density_natural = 10^(effect_plots_dt$density*mean_sd[2] + mean_sd[1])

  effect_plots_temp_dt = do.call(rbind, effect_plots_temp) 
  effect_plots_temp_dt$temperature_natural = effect_plots_temp_dt$temperature*mean_sd_temp[2] + mean_sd_temp[1]

  #scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))

  col = '#1b9e77'
  p1 = ggplot(effect_plots_dt[param == 'omega']) + geom_line(aes(x=density_natural, y=med), color=col) + 
                            geom_ribbon(aes(x=density_natural, ymin=lower, ymax=upper), fill=col, alpha=0.3) +
                            scale_x_log10() +
                            ylim(c(0, 1)) + theme_classic() + ylab("Yearly probability of remaining\nin low prevalence state (\u03C9)") + 
                            xlab("Host abundance") + annotate("text", x = 10, y = .25, size=2.5, label = "Weak effect,\nopposite direction of\na priori prediction")

  p2 = ggplot(effect_plots_dt[param == 'phi']) + geom_line(aes(x=density_natural, y=med), color=col) + 
                            geom_ribbon(aes(x=density_natural, ymin=lower, ymax=upper), fill=col, alpha=0.3) +
                            scale_x_log10() +
                            ylim(c(0, 1)) + theme_classic() + ylab("Yearly probability of\ntransitioning to Bd-free state (\u03C6)") + 
                            xlab("Host abundance") + annotate("text", x = 10, y = .25, size=2.5, label = "Weak effect,\nsame direction as\na priori prediction")

  p3 = ggplot(effect_plots_temp_dt[param == 'omega']) + geom_line(aes(x=temperature_natural, y=med), color=col) + 
                            geom_ribbon(aes(x=temperature_natural, ymin=lower, ymax=upper), alpha=0.3, fill=col) +
                            ylim(c(0, 1)) + theme_classic() + ylab("Yearly probability of remaining\nin low prevalence state (\u03C9)") + 
                            xlab("Maximum summer air temperature (C)")

  p4 = ggplot(effect_plots_temp_dt[param == 'phi']) + geom_line(aes(x=temperature_natural, y=med), color=col) + 
                            geom_ribbon(aes(x=temperature_natural, ymin=lower, ymax=upper), alpha=0.3, fill=col) +
                            ylim(c(0, 1)) + theme_classic() + ylab("Yearly probability of \ntransitioning to Bd-free state (\u03C6)") + 
                            xlab("Maximum summer air temperature (C)")
  eplot = (p1 + p2) / (p3 + p4) + plot_annotation(tag_levels="A", tag_suffix=".")
  ggsave("../results/effect_plots.jpg", width=7, height=5)

  
}


#######################################
#### Tadpole and temperature plots ####
#######################################

# Plot tadpole and temperature plot
if(make_tadpole_plots){

  
  params = list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                non_linear_temp=FALSE, log_tadpole_abundance=FALSE,
                include_snow_depth=TRUE)
  all_res = set_up_response_and_covariate_data(params, full_data, temperature_data, snow_data)
  list2env(all_res, .GlobalEnv)
  cutoff = 0.5
  a_or_d = "abund"
  pname = "phi"

  stan_data = prepare_stan_data(num_inf_mat, sample_size_mat, abund_missing_mat, cutoff,
                                imputed_cov_mats, a_or_d)
  saveRDS(stan_data$included_sites, "../results/analysisI_sites.rds")

  loss_dt = naive_probs_with_logistic_regression(stan_data$obs_mat, stan_data$density, pname, 
                                                 temp_mat=stan_data$temperature, 
                                                 tad_mat=stan_data$tadpole,
                                                 snow_mat=stan_data$snow)$data


  loss_dt$site_id = stan_data$included_sites[loss_dt$lid]
  trun_dt = loss_dt[loss_dt$tad %in% c(0), ]
  counts = as.data.table(loss_dt)[, .(prop=sum(y) / length(y), count=length(y)), by=.(tad)]
  tad_plot = ggplot(loss_dt, aes(x=tad, y=y)) + geom_jitter(width=0.05, height=0.05, color='gray', alpha=0.5) + stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE) + 
                                     # geom_label(data=trun_dt, aes(x=tad, y=seq(0.2, 0.8, len=nrow(trun_dt)), label=site_id)) +
                                     geom_point(data=counts, aes(x=tad, y=prop, size=count + 10), color='black', shape=17) +
                                     guides(size='none') +
                                     theme_classic() + xlab("Tadpole absence/presence") + 
                                     ylab("Observed loss of Bd in time step\nphi: 1=loss, 0=no loss")
  dens_plot = ggplot(loss_dt, aes(x=density, y=y)) + geom_jitter(width=0.05, height=0.05) + stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)                               

  # Look at tadpole sites
  tad_sites = loss_dt[loss_dt$tad %in% c(0), ]$site_id
  # tad_mat[site_id %in% tad_sites]
  # obs_bd_mat[site_id %in% tad_sites]
  # fit = glm(y ~ density + temp + tad + snow + temp:snow, data=loss_dt, family="binomial")

  # Make temperature plots

  loss_dt = naive_probs_with_logistic_regression(stan_data$obs_mat, stan_data$density, "omega", 
                                                 temp_mat=stan_data$temperature, 
                                                 tad_mat=stan_data$tadpole,
                                                 snow_mat=stan_data$snow)$data

  loss_dt$site_id = stan_data$included_sites[loss_dt$lid]
  loss_dt = as.data.table(loss_dt)

  lake_dat = fread("../data/archival/sites_nov2020.csv")
  temp = merge(loss_dt, lake_dat[, .(id, elevation)], by.x="site_id", by.y="id")

  # Show interaction for severe and less severe winter
  loss_dt$site_id = stan_data$included_sites[loss_dt$lid]
  temp_plot = ggplot() + geom_jitter(data=loss_dt[loss_dt$snow < 0, ], aes(x=temp, y=y, color="Less severe winter"), height=0.05) + 
                         geom_jitter(data=loss_dt[loss_dt$snow > 0 & loss_dt$snow < 1, ], aes(x=temp, y=y, color="More severe winter"), height=0.05) +
                                     stat_smooth(data=loss_dt[loss_dt$snow < 0, ], aes(x=temp, y=y, color="Less severe winter", fill="Less severe winter"), method="glm", method.args=list(family="binomial"), se=TRUE) + 
                                     stat_smooth(data=loss_dt[loss_dt$snow > 0 & loss_dt$snow < 1, ], aes(x=temp, y=y, color="More severe winter", fill="More severe winter"), method="glm", method.args=list(family="binomial"), se=TRUE) + 
                                     theme_classic() + xlab("Maximum summer temperature (standardized)") + 
                                     ylab("Probability of remaining in low prevalence state\nomega: 1=remain, 0=to high prevalence") +
                                     guides(fill="none") +
                                     theme(legend.position=c(0.7,0.25),
                                           legend.title = element_blank(),
                                           legend.background=element_rect(fill="transparent"))
      

  # Combine tadpole and temperature plots
  tplot = temp_plot + tad_plot + plot_annotation(tag_levels="A", tag_suffix=".")
  ggsave("../results/temp_tadpole_plot_example.pdf", plot=tplot, width=8, height=4)

}



#############################################################
#### Explore the no density effect with model comparison ####
#############################################################

if(explore_density_effect){

  # Run this for density and abundance and across rhos
  params = list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
                non_linear_temp=FALSE, log_tadpole_abundance=FALSE,
                include_snow_depth=TRUE)
  all_res = set_up_response_and_covariate_data(params, full_data, temperature_data, snow_data)
  list2env(all_res, .GlobalEnv)

  cutoffs = c(0.25, 0.3333, 0.5)
  looic_vals = array(NA, dim=c(length(cutoffs), 2))
  bf_vals = array(NA, dim=c(length(cutoffs)))
  for(i in 1:length(cutoffs)){


    cutoff = cutoffs[i]
    a_or_d = "density"
    stan_data = prepare_stan_data(num_inf_mat, sample_size_mat, abund_missing_mat, cutoff,
                                  imputed_cov_mats, a_or_d)
    # Density models

    # Fit null model with no density effect
    null_mod = stan_model(file="stan_files/hmm_occupancy_fourstates_covariates_imputation_no_density.stan")
    fit_null_mod = sampling(null_mod, data=stan_data, iter=2000, warmup=1000, chains=3, cores=3,
                            include=TRUE)#, pars=c("log_lik"))
    fit_null_mod_marginal = bridge_sampler(fit_null_mod)
    error_measures(fit_null_mod_marginal)

    # Fit model with negative density effect
    neg_mod = stan_model(file="stan_files/hmm_occupancy_fourstates_covariates_imputation_negative_density.stan")
    fit_neg_mod = sampling(neg_mod, data=stan_data, iter=2000, warmup=1000, chains=3, cores=3,
                           include=TRUE)#, pars=c("log_lik"))
    fit_neg_mod_marginal = bridge_sampler(fit_neg_mod)

    # Compute baye factor
    bayes_factor = bf(fit_null_mod_marginal, fit_neg_mod_marginal)
    bf_vals[i] = bayes_factor$bf

    # loovals_null = loo(fit_null_mod)
    # loovals_neg = loo(fit_neg_mod)
    # looic_vals[i, 1] = loovals_null$estimates[3, 1]
    # looic_vals[i, 2] = loovals_neg$estimates[3, 1]
  }

  # Save LOO and Bayes Factor results
  # saveRDS(looic_vals, "../results/looic_values_density_models_no_rarify_all_cutoffs.rds")
  saveRDS(bf_vals, "../results/bayesfactor_values_density_models_no_rarify_all_cutoffs.rds")
}


###################################
#### Plot state trajectories   ####
###################################

if(example_plots){

  combo = list(rarify=FALSE, subsample=NULL, include_subadults=TRUE, 
               non_linear_temp=FALSE, log_tadpole_abundance=FALSE, 
               include_snow_depth=TRUE)
  a_or_d = "density" 
  fitted_data_temp = readRDS(paste0("../results/fitted_data_rarify_", 
                                    combo$rarify, combo$subsample, "_include_subadults", 
                                    combo$include_subadults, 
                                    "_nonlinear_temp", combo$non_linear_temp,
                                    "_tadpolelog10", combo$log_tadpole_abundance, 
                                    "_snowdepth", combo$include_snow_depth, ".rds"))

  fitted_data = fitted_data_temp[[a_or_d]][[3]] # Cutoff 0.5
  sites = fitted_data$included_sites

  # Set up data for plotting
  obs_dat = fitted_data$obs_mat
  obs_dat[obs_dat == -1] = NA
  mdat = as.data.table(data.table::melt(t(obs_dat)))
  colnames(mdat) = c("year", "id", "state")
  mdat = mdat[!is.na(state)]
  id_dt = as.data.table(data.frame(id=1:nrow(obs_dat), lake_id=sites))
  mdat = merge(mdat, id_dt, key="id")

  # Plot the state trajectories
  ggplot(mdat[id %in% 10:25], aes(x=year, y=state)) + geom_path() + geom_point() + facet_wrap(~lake_id) + theme_bw() +   
                theme(strip.text.x = element_text(size=9),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                scale_y_continuous(labels=c("U", "LP", "HP", "E"), 
                                   breaks=0:3, 
                                   limits=c(0, 3), 
                                   minor_breaks=0:3) +
                scale_x_continuous(breaks=seq(2005, 2020, by=2), limits=c(2004, 2019)) +
                xlab("Year") + ylab("Observed state")

  ggsave("../results/state_trajectories.pdf", width=8, height=7)

}


