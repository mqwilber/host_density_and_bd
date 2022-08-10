library(data.table)
library(rstan)
library(ggplot2)
# install.packages("patchwork", version="1.1.1")
library(patchwork)

## The script addressed the question: 
## Does the initial abundance or density of the amphibian population affect the magnitude of disease-induced decline?

candidate_sites = function(decline_dat){
  # Generate data set with sites that have shown some decline
  #
  # decline_dat : data.table
  #
  # Returns
  # -------
  # : data.table with same columns as decline_dat but with sites removed based on criteria

  decline_dat = decline_dat[(year_max < year_min) & 
                          (prev_overall_max >= 0.25) & 
                          (year_prev_max >= year_max) & 
                          (year_prev_max <= year_min)]

  # Drop sites identified by Roland
  exclude_sites = c(10055, 10276, 10314, 10315, 10316, 12102, 21111)
  decline_dat = decline_dat[!(site_id %in% exclude_sites)]
  decline_dat = decline_dat[abund_max > 0]
  return(decline_dat)

}

check_convergence = function(fit){
  # Check basic convergence diagnostics of model
  #
  # fit: stan model


  fit_summary = summary(fit1)
  print(all((fit_summary$summary[, "Rhat"] < 1.01), na.rm=T))
  
  # Which Rhat aren't converging well?
  Rhat = fit_summary$summary[, "Rhat", drop=T]
  print(Rhat[(Rhat > 1.01) & !(is.na(Rhat))])
  
  # Check effective sample sizes of parameters
  n_eff = fit_summary$summary[, 'n_eff', drop=T]
  print(n_eff[(n_eff < 400) & !(is.na(n_eff))])
}


plot_fits = function(fit, name, adult_nm, standata, legend){
  # Plot the best fit results
  #
  # fit : stan fit
  # name : str, "abundance" or "density"
  # adult_nm : str: "adults" or "adults + subadults"
  # standata : list of data passed into stan model
  # legend : bool, TRUE include legend, FALSE no legend
  #
  # Returns
  # -------
  # ggplot object
  
  # Check out predicted relationship
  beta0 = extract(fit, pars='beta0')$beta0
  beta1 = as.vector(extract(fit, pars='beta')$beta[, 1])

  # Density effect
  lower_beta1 = round(quantile(beta1, 0.025), 2)
  upper_beta1 = round(quantile(beta1, 0.975), 2)
  
  mean_x = standata$mean_x
  sd_x = standata$sd_x
  lower = min(standata$X[, 1])
  upper = max(standata$X[, 1])
  vals = seq(lower, upper, len=50)
  
  pred = sapply(vals, function(v) beta0 + beta1*v)
  prob = 1 / (1 + exp(-pred))
  probs = apply(prob, 2, median)
  probs_lower = apply(prob, 2, quantile, 0.025)
  probs_upper = apply(prob, 2, quantile, 0.975)
  
  xvals = as.vector(standata$X[, 1])*sd_x + mean_x
  yvals = 1 - (standata$end_abund / standata$start_abund)

  txvals = vals*sd_x + mean_x
  print(txvals)
  tplot = ggplot() + geom_line(aes(x=txvals, y=probs)) + 
    geom_ribbon(aes(x=txvals, ymin=probs_lower, ymax=probs_upper), alpha=0.25) +
    geom_point(aes(x=xvals, y=yvals)) + 
    annotate("text", x = quantile(txvals, 0.7), y = 0.05, 
             label = paste0("Density/abundance effect\n", "95% CI: [", lower_beta1, ", ", upper_beta1, "]"),
             size=3) +
    theme_classic() + xlab(paste("log10", adult_nm, name)) + 
    ylab("Magnitude of decline") + ylim(c(0, 1))
  
  if(legend){
    tplot = tplot + theme(legend.position = c(0.7, 0.2))
  } else{
    tplot = tplot + theme(legend.position = "none")
  }
  return(tplot)

}


#########################
###### Load data ########
#########################


decline_dat_wsub = fread("../data/candidate_datasets/analysisII_magnitudeofdecline_withsubadults.csv")
decline_dat_nosub = fread("../data/candidate_datasets/analysisII_magnitudeofdecline_withoutsubadults.csv")

nrow(decline_dat_wsub)
nrow(decline_dat_nosub)

# Set of sites that don't overlap
setdiff(decline_dat_wsub$site_id, decline_dat_nosub$site_id)

# Set up datasets to compare
datasets = list(wsub=decline_dat_wsub, nosub=decline_dat_nosub)
ddats = list()
for(nm in names(datasets)){

  ddat = datasets[[nm]]
  ddats[[nm]] = candidate_sites(ddat)

}

nrow(ddats[["wsub"]]) # 82 sites
nrow(ddats[["nosub"]]) # 52 sites

######################
#### Fit model #######
######################

# Set up model
bbmod = stan_model("stan_files/betabinomial_decline.stan")

mod_fits = list()
all_standata = list()
for(nm in names(ddats)){

  decline_dat = ddats[[nm]]

  # Set up design matrices 
  X1_abund = model.matrix(lm(abund_min ~ scale(log10(abund_max)) - 1, data=decline_dat))
  X1_density = model.matrix(lm(abund_min ~ scale(log10(density)) - 1, data=decline_dat))

  standata1 = list(N=nrow(decline_dat),
                  P=1,
                  X=X1_abund,
                  mean_x=mean(log10(decline_dat$abund_max)),
                  sd_x=sd(log10(decline_dat$abund_max)),
                  start_abund=as.integer(decline_dat$abund_max),
                  end_abund=as.integer(decline_dat$abund_min))

  standata2 = list(N=nrow(decline_dat),
                  P=1,
                  X=X1_density,
                  mean_x=mean(log10(decline_dat$density)),
                  sd_x=sd(log10(decline_dat$density)),
                  start_abund=as.integer(decline_dat$abund_max),
                  end_abund=as.integer(decline_dat$abund_min))

  fit1 = sampling(bbmod, data=standata1, chains=3, cores=3)
  fit2 = sampling(bbmod, data=standata2, chains=3, cores=3)

  cat("Checking converged of abundance model", "\n")
  check_convergence(fit1)
  cat("Checking converged of density model", "\n")
  check_convergence(fit2)

  mod_fits[[nm]] = list(abund=fit1, density=fit2)
  all_standata[[nm]] = list(abund=standata1, density=standata2)
}

##########################
#### Plot results  #######
##########################

all_plots = list()
count = 1
for(nm in names(mod_fits)){

  fits = mod_fits[[nm]]
  data = all_standata[[nm]]
  if(nm == "wsub"){
    adult_nm = "adult + subadult"
  } else{
    adult_nm = "adult"
  }

  abund_plot = plot_fits(fits[['abund']], "abundance", adult_nm, data[['abund']], TRUE)
  density_plot = plot_fits(fits[['density']], "density (per m)", adult_nm, data[['density']], FALSE)
  all_plots[[count]] = abund_plot
  all_plots[[count + 1]] = density_plot
  count = count + 2
}

tp = (all_plots[[3]] + all_plots[[4]]) / (all_plots[[1]] + all_plots[[2]]) + plot_annotation(tag_levels = "A", tag_suffix = ".")
tp
ggsave("../results/magnitude_of_decline.pdf", width=7, height=6)


#########################################
#### Predict decline percentages ########
#########################################


num_frogs = 7
standata = all_standata[['wsub']]$abund
lower_abund = (log10(num_frogs) - standata$mean_x) / standata$sd_x
beta0 = extract(mod_fits[['wsub']]$abund, pars="beta0")$beta0
beta1 = extract(mod_fits[['wsub']]$abund, pars="beta")$beta[, 1]

# Prediction at 7 frogs
pred = beta0 + beta1*lower_abund
probs = (1 / (1 + exp(-pred)))
quantile(probs, c(0.025, 0.5, .975))






