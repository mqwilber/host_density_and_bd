library(data.table)
library(expm)
library(rstan)

naive_probs = function(sim_res){
	# Calculate the naive probabilities of a transition by summing
	# transitions

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

	# Arrival
	arrival_df = all_df[first == 0]
	a = 1 - (nrow(arrival_df[second == 0]) / nrow(arrival_df))
	arrival_df2 = arrival_df[second != 0]
	zeta = nrow(arrival_df2[second == 1]) / nrow(arrival_df2)

	# Loss df
	loss_df = all_df[first == 1]
	phi = nrow(loss_df[second == 0]) / nrow(loss_df)
	loss_df2 = loss_df[second != 0]
	omega = nrow(loss_df2[second == 1]) / nrow(loss_df2)

	# Extinction df
	ext_df = all_df[first == 2]
	ext = nrow(ext_df[second == 3]) / nrow(ext_df)
	params = c(a=a, zeta=zeta, phi=phi, omega=omega, e=ext)
	return(params)

}


## Building a toy model to test our inference

# Parameters
a = 0.4
zeta = 0.5
phi = 0.4
omega = 0.4
e = 0.2

# Set up initial condition and transition matrix
init = c(1, 0, 0, 0)
Pmat = matrix(c(1 - a, a*(zeta), a*(1 - zeta), 0,
			   phi, (omega)*(1 - phi), (1 - omega)*(1 - phi), 0,
			   0, 0, (1 - e), e, 
			   0, 0, 0, 1), nrow=4, ncol=4, byrow=TRUE)

num_reps = 1
bias_results_new = array(NA, dim=c(num_reps, 5))

# Loop over multiple simulations
for(b in 1:num_reps){

	tsteps = 15
	nsims = 20

	all_results = array(NA, dim=c(nsims, tsteps + 1))

	# Loop over multiple "lakes"
	for(n in 1:nsims){

		results = array(NA, dim=c(4, tsteps + 1))
		results[, 1] = init

		# Run stochastic simulation
		for(i in 2:(tsteps + 1)){

			old_val = results[, i - 1]
			probs = old_val %*% Pmat
			new_val = rmultinom(1, 1, probs)
			results[, i] = new_val

		}

		# Process which items are ones and zeros
		collapsed_res = apply(results, 2, function(x) which.max(x)) - 1
		all_results[n, ] = collapsed_res
	}

	# Assign missing at random
	total_missing = (tsteps + 1)*(nsims)
	num_missing = 1 #floor((tsteps - 1) * nsims * 0.5)
	missing_inds = sample(31:total_missing, num_missing, replace=F)
	all_results[missing_inds] = -1

	# Pass to stan for processing
	mod1 = stan_model(file="stan_files/hmm_occupancy_fourstates_test.stan")

	stan_data = list(Y=tsteps + 1,
					 L=nsims,
					 obs_mat=all_results)

	fit1 = sampling(mod1, data=stan_data, iter=1000, warmup=500, chains=4, cores=4,
				   include=TRUE, pars=c("a", "zeta", "phi", "omega", "ext"))

	# That looks good.  Recovering the true parameters
	bias_results_new[b, ] = summary(fit1)$summary[1:5, 1]
}

naive_probs(all_results)

## Let's add in additional observation error and try to account for it

cutoff = 0.5

# Assign sample sizes...number of frogs sampled
sample_size = array(sample(19:20, total_missing, replace=T), dim=dim(all_results))
sample_size[all_results == c(3)] = 0
sample_size[all_results == -1] = -1

# Assign true prevalences
true_prev = array(NA, dim=dim(sample_size))
low_prev = all_results == 1
true_prev[low_prev] = runif(sum(low_prev), min=0, max=cutoff) 
high_prev = all_results == 2
true_prev[high_prev] = runif(sum(high_prev), min=cutoff, max=1)

# Draw observed numer of infected hosts sampled
num_inf_sampled = array(NA, dim=dim(sample_size))
bd_present_ind = !is.na(true_prev)
num_inf_sampled[bd_present_ind] = rbinom(sum(bd_present_ind), 
										 sample_size[bd_present_ind], 
										 true_prev[bd_present_ind])
num_inf_sampled[all_results == -1] = -1
num_inf_sampled[all_results == 0] = 0
num_inf_sampled[all_results == 3] = 0

obs_prev = num_inf_sampled / sample_size

# Classify the "Observed" state data
all_results_obs = array(NA, dim=dim(sample_size))
all_results_obs[obs_prev == 0] = 0
all_results_obs[obs_prev > 0 & obs_prev <= cutoff] = 1
all_results_obs[obs_prev > cutoff] = 2
all_results_obs[all_results == 3] = 3
all_results_obs[all_results == -1] = -1

# Pass to stan for processing
mod2 = stan_model(file="stan_files/hmm_occupancy_fourstates_test_obs_process.stan")

stan_data = list(Y=tsteps + 1,
				 L=nsims,
				 cutoff=cutoff,
				 obs_mat=all_results_obs,
				 sample_size=sample_size,
				 num_inf=num_inf_sampled)

fit3 = sampling(mod2, data=stan_data, iter=1000, warmup=500, chains=4, cores=4,
			   include=TRUE, pars=c("a", "zeta", "phi", "omega", "ext"))

# Close
fit3









