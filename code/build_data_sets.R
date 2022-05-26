library(dplyr)
library(data.table)
library(lubridate)

max_double = function(x1, x2){
  # Function for returning specific maximum values from two different
  # vectors when either one or the other might have all missing values
  # Either take the max of one vector, the other, or both.
  #
  # Parameters
  # ----------
  # x1 : vector
  # x2 : vector

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


max_double_na = function(x1, x2){
  # Same as max_double, but missing value is NA instead of -1

  # Function for 
  if(all(is.na(x1))){
    maxx1 = as.integer(NA)
  } else{
    maxx1 = max(x1, na.rm=T)
  }

  if(all(is.na(x2))){
    maxx2 = as.integer(NA)
  } else{
    maxx2 = max(x2, na.rm=T)
  }

  if(is.na(maxx1)){
    return(maxx2)
  } else if(is.na(maxx2)){
    return(maxx2)
  } else{
    return(maxx1 + maxx2)
  }
}

#### Generate candidate data sets for the three analyses used in manuscript ####

# Load raw data provided by Roland
count_dat = fread("../data/archival/ramu_counts_nov2020.csv")
load_dat = fread("../data/archival/ramu_bdload_nov2020.csv")
site_dat = fread("../data/archival/sites_nov2020.csv")

# Set specific year column
count_dat[, datetime:=as.POSIXct(visit_date, format="%Y-%m-%d")][, year:=year(datetime)]
load_dat[, datetime:=as.POSIXct(visit_date, format="%Y-%m-%d")][, year:=year(datetime)]

# Get maximum counts per site per year for each lifestage.
# NOTE: This is where we account for multiple surveys per year
cast_counts = count_dat[, .(amphib_max=max(stage_count)), by=.(site_id, year, visual_life_stage)][
                        visual_life_stage != ""] %>%  
                        dcast(site_id + year ~ visual_life_stage,  value.var="amphib_max", fill=0)

# Check on complete extinction...only two sites 
cast_counts[adult == 0 & eggmass == 0 & subadult == 0 & tadpole == 0]

# Merge load data and count data.
# Only include data with swabs
roland_full1 = merge(load_dat, cast_counts, by=c("site_id", "year"), all.x=T)

# Include data with swabs OR VES, but not necessarily both
roland_full2 = merge(load_dat, cast_counts, by=c("site_id", "year"), all.x=T, all.y=T)

# Save data for analysis
fwrite(roland_full1, "../data/formatted/combined_Bd_data_for_analysis.csv")
fwrite(roland_full2, "../data/formatted/combined_Bd_data_for_analysis_allVES.csv")

######################################################
#### Analysis I: Building failed invasion dataset ####
######################################################

full_data = fread("../data/formatted/combined_Bd_data_for_analysis_allVES.csv")
lake_dat = fread("../data/archival/sites_nov2020.csv")
colnames(lake_dat)[1] = "site_id"

# Dataset that identifies all lakes that experienced candidate failed invasions
# This data set has been manually processed by Roland to identify which lakes
# likely experienced a failed invasion. See make_candidate_failed_invasions.Rmd
rolands_priors = fread("../data/formatted/candidate_failed_invasion_sites_reduced_roland_additions.csv")[in_roland == TRUE]

# Should we include or exclude tadpoles when computing mean loads?
# Including tadpoles in the load calculation for this analysis
red_data = full_data[, .(bd_present=ifelse(all(is.na(bd_load)), as.integer(-1), as.integer(any(bd_load > 0))), 
                        num_inf=ifelse(all(is.na(bd_load)), as.integer(0), sum(bd_load > 0)), 
                        num_samps=ifelse(all(is.na(bd_load)), as.integer(0), length(bd_load)),
                        prev=ifelse(all(is.na(bd_load)), as.numeric(-1), sum(bd_load > 0) / length(bd_load)),
                        tadpole=ifelse(all(is.na(tadpole)), as.integer(-1), max(tadpole, na.rm=T)),
                        mean_bd=ifelse(all(is.na(bd_load)), as.numeric(-1), ifelse(any(bd_load > 0), mean(log10(bd_load[bd_load > 0])), -1)), # Keep this on the natural scale?
                        frog_abundance=max_double(adult, subadult)), by=.(site_id, year)][year >= 2004]

# Join datasets
red_data = merge(red_data, lake_dat[, .(site_id, jurisdiction)], by="site_id", all.x=TRUE)
obs_bd_mat = dcast(red_data, site_id ~ year,  value.var="bd_present", fill=-1)
obs_bd_mat = merge(obs_bd_mat, lake_dat[, .(site_id, jurisdiction)], by="site_id", all.x=TRUE)

# Which sites had the first sample as a 0 and a later sample as 1
num_cols = ncol(obs_bd_mat)
good_lakes = list()
for(i in 1:nrow(obs_bd_mat)){

	rowvals = obs_bd_mat[i, 2:(num_cols - 1)]
	first_present = which(rowvals == 1)
	first_absent = which(rowvals == 0)

	# Was there ever any Bd absence?
	if(length(first_absent) > 0){

		# Was there ever any Bd presence?
		if(length(first_present > 0)){

			if(first_absent[1] < first_present[1]){
				good_lakes = append(good_lakes, obs_bd_mat[i, 1])
			}
		}

	}
}

candidate_lakes = unlist(good_lakes)

# Add lakes that had putative failed invasions but that didn't explicitly fit first sample 0 and later
# sample 1 criteria. See 
potential_lake_ids = unique(c(candidate_lakes, rolands_priors$site_id))
potential_lakes = obs_bd_mat[site_id %in% potential_lake_ids]
potential_lakes = merge(potential_lakes, rolands_priors[, .(site_id, loss_years, failed_invasion_ranking)], all.x=TRUE)

# Additional lakes to exclude based on Roland's assessments in his email.
# See docs/other_docs/rolands_notes_analysisI.txt
exclude = c(10277, 11506, 10315, 21520, 52103, 
            84303, 84312, 84313, 84326, 84327, 
            84329, 70215, 70224, 70508, 70629, 
            70634, 71967, 72336, 72554, 72708, 
            72990, 72998, 74062)

 # Drop all yosemite lakes based on Roland's suggestions
potential_lakes = potential_lakes[!(site_id %in% exclude)][jurisdiction != "yosemite"]

# Analysis I dataset after exclusions
fwrite(potential_lakes, "../data/candidate_datasets/analysisI_failedinvasions_excluded.csv")


####################################################
#### Analysis II: Magnitude of decline analysis ####
####################################################

full_data = fread("../data/formatted/combined_Bd_data_for_analysis_allVES.csv")

# Group data by lake and year
include_subadults = TRUE
include_only_adult_swabs = FALSE

if(include_subadults){

  red_data = full_data[, .(bd_present=ifelse(all(is.na(bd_load)), as.integer(NA), as.integer(any(bd_load > 0, na.rm=TRUE))), 
                          num_inf=ifelse(all(is.na(bd_load)), as.integer(0), sum(bd_load > 0, na.rm=TRUE)), 
                          num_samps=ifelse(all(is.na(bd_load)), as.integer(0), length(bd_load)),
                          prev=ifelse(all(is.na(bd_load)), as.numeric(NA), sum(bd_load > 0, na.rm=TRUE) / length(bd_load)),
                          tadpole=ifelse(all(is.na(tadpole)), as.integer(NA), as.integer(any(tadpole > 0, na.rm=TRUE))),
                          mean_bd=ifelse(all(is.na(bd_load)), as.numeric(NA), ifelse(any(bd_load > 0, na.rm=TRUE), mean(log10(bd_load[bd_load > 0] + 1), na.rm=TRUE), 0)), # Keep this on the natural scale?
                          frog_abundance=max_double_na(adult, subadult)), by=.(site_id, year)][year >= 2004]


} else {

  if(include_only_adult_swabs){

      full_data$bd_load[full_data$capture_life_stage != "adult"] = as.numeric(NA)
      red_data = full_data[, .(bd_present=ifelse(all(is.na(bd_load)), as.integer(NA), as.integer(any(bd_load > 0, na.rm=TRUE))), 
                          num_inf=ifelse(all(is.na(bd_load)), as.integer(0), sum(bd_load > 0, na.rm=TRUE)), 
                          num_samps=ifelse(all(is.na(bd_load)), as.integer(0), length(bd_load)),
                          prev=ifelse(all(is.na(bd_load)), as.numeric(NA), sum(bd_load > 0, na.rm=TRUE) / length(bd_load)),
                          tadpole=ifelse(all(is.na(tadpole)), as.integer(NA), as.integer(any(tadpole > 0, na.rm=TRUE))),
                          mean_bd=ifelse(all(is.na(bd_load)), as.numeric(NA), ifelse(any(bd_load > 0, na.rm=TRUE), mean(log10(bd_load[bd_load > 0] + 1), na.rm=TRUE), 0)), # Keep this on the natural scale?
                          frog_abundance=ifelse(all(is.na(adult)), as.integer(NA), max(adult, na.rm=T))), by=.(site_id, year)][year >= 2004]

  } else {
    red_data = full_data[, .(bd_present=ifelse(all(is.na(bd_load)), as.integer(NA), as.integer(any(bd_load > 0, na.rm=TRUE))), 
                            num_inf=ifelse(all(is.na(bd_load)), as.integer(0), sum(bd_load > 0, na.rm=TRUE)), 
                            num_samps=ifelse(all(is.na(bd_load)), as.integer(0), length(bd_load)),
                            prev=ifelse(all(is.na(bd_load)), as.numeric(NA), sum(bd_load > 0, na.rm=TRUE) / length(bd_load)),
                            tadpole=ifelse(all(is.na(tadpole)), as.integer(NA), as.integer(any(tadpole > 0, na.rm=TRUE))),
                            mean_bd=ifelse(all(is.na(bd_load)), as.numeric(NA), ifelse(any(bd_load > 0, na.rm=TRUE), mean(log10(bd_load[bd_load > 0] + 1), na.rm=TRUE), 0)), # Keep this on the natural scale?
                            frog_abundance=ifelse(all(is.na(adult)), as.integer(NA), max(adult, na.rm=T))), by=.(site_id, year)][year >= 2004]
  }

}

# Cast into wide form for analysis and exploration
obs_bd_mat = dcast(red_data, site_id ~ year,  value.var="bd_present", fill=-1)
sample_size_mat = dcast(red_data, site_id ~ year,  value.var="num_samps", fill=-1)
num_inf_mat = dcast(red_data, site_id ~ year,  value.var="num_inf", fill=-1)
bd_mat = dcast(red_data, site_id ~ year,  value.var="mean_bd", fill=-1)
abund_mat = dcast(red_data, site_id ~ year,  value.var="frog_abundance", fill=-1)
tad_mat = dcast(red_data, site_id ~ year,  value.var="tadpole", fill=-1)

num_yrs = 3 # Number of samples needed for inclusion
unq_lakes = unique(red_data$site_id)

all_res = list()
count_greater_than = 0

for(id in unq_lakes){
  
  tdat = red_data[site_id == id]
  tdat = tdat[year >= 2004]

  # Include lakes with enough samples
  if((sum(!is.na(tdat$frog_abundance)) >= num_yrs) & (sum(!is.na(tdat$bd_present)) >= num_yrs)){


    if((length(tdat$frog_abundance) != sum(is.na(tdat$frog_abundance))) & (length(tdat$bd_present) != sum(is.na(tdat$bd_present)))){

      ind_max = which(tdat$frog_abundance == max(tdat$frog_abundance, na.rm=T))
      ind_min =  which(tdat$frog_abundance == min(tdat$frog_abundance, na.rm=T))
      ind_max_bd = which(tdat$mean_bd == max(tdat$mean_bd, na.rm=T))

      # Make sure some frog abundance was available
      # Account for multiple mins and maxs
      if(length(ind_max) >= 1){
        
        print(c(length(ind_max), length(ind_min)))
        count_greater_than = count_greater_than + 1
        year_max = min(tdat$year[ind_max]) # Take the minimum year
        year_max_bd = min(tdat$year[ind_max_bd])
        year_vals =  year_max_bd:max(tdat$year)
        abund_min_vals = tdat[year %in% year_vals]$frog_abundance

        # Look at years after Bd max. If any non-NA abundances, return year of minimum
        # abundance. Otherwise, return the most recent year with minimum abundance for all
        # surveys
        year_min = ifelse(all(is.na(abund_min_vals)), max(tdat$year[ind_min]), 
                          tdat[year %in% year_vals][frog_abundance == min(frog_abundance, na.rm=T)]$year)
        bd_abund_max = tdat$mean_bd[which(tdat$year == year_max)]
        prev_max = tdat$prev[which(tdat$year == year_max)]
        prev_overall_max = max(tdat$prev, na.rm=T)
        year_prev_overall_max = tdat$year[which.max(tdat$prev)]
        bd_max = tdat$mean_bd[ind_max_bd][1]

        abund_max = tdat$frog_abundance[ind_max][1]
        abund_min = tdat$frog_abundance[which(tdat$year == year_min)]
        tabund = tdat$frog_abundance[tdat$frog_abundance != -1]
        abund_cv = mean(tabund) / sd(tabund)
        
        tres = data.frame(site_id=id, year_max=year_max, year_min=year_min, 
                          year_max_bd=year_max_bd, 
                          bd_max_when_abund_max_log10=bd_abund_max,
                          bd_max=bd_max,
                          prev_when_abund_max=prev_max,
                          prev_overall_max=prev_overall_max,
                          year_prev_max=year_prev_overall_max,
                          abund_max=abund_max, abund_min=abund_min)
        
        all_res[[id]] = tres
        
      }
    }
    
  }
}

decline_dat = as.data.table(do.call(rbind, all_res))
decline_dat = decline_dat[abund_max >= 1]

# Extract lake data to calculate perimeter
lake_dat = fread("../data/archival/sites_nov2020.csv")
colnames(lake_dat)[1] = "site_id"
decline_dat = merge(decline_dat, lake_dat[, .(site_id, area, jurisdiction, elevation)], key="site_id")
decline_dat[, density:=(abund_max / sqrt(area / pi)*pi*2)]

# Check certain test sites are present 
check_sites = c(11506, 10081, 10100, 10101, 10102, 11858, 10090, 10091,
                10206, 10222, 10284, 10538, 10537, 10571, 11506, 11525, 11526,
                20198, 21091)

for(cs in check_sites){
  cat(cs, any(decline_dat$site_id == cs), "\n")
}

# Make sure site had some Bd
decline_dat = decline_dat[bd_max != 0]

# Only include analysis I sites
decline_dat = decline_dat[jurisdiction %in% c("kings_canyon", "sequoia")]

if(include_subadults){
  fwrite(decline_dat, "../data/candidate_datasets/analysisII_magnitudeofdecline_withsubadults.csv")
} else{
  fwrite(decline_dat, "../data/candidate_datasets/analysisII_magnitudeofdecline_withoutsubadults.csv")
}

############################################################
#### Analysis III: Effects of density on infection load ####
############################################################

full_data = fread("../data/formatted/combined_Bd_data_for_analysis_allVES.csv")

# Do not include tadpole loads in this analysis
full_data$bd_load[full_data$capture_life_stage == 'tadpole'] = NA
red_data = full_data[, .(bd_present=ifelse(all(is.na(bd_load)), as.integer(NA), as.integer(any(bd_load > 0, na.rm=TRUE))), 
                         num_inf=ifelse(all(is.na(bd_load)), as.integer(0), as.integer(sum(bd_load > 0, na.rm=TRUE))),
                         num_samps=ifelse(all(is.na(bd_load)), as.integer(0), length(bd_load)),
                         mean_bd=ifelse(all(is.na(bd_load)), as.numeric(NA), ifelse(any(bd_load > 0, na.rm=TRUE), as.numeric(mean(log(bd_load[bd_load > 0] + 1), na.rm=TRUE)), as.numeric(NA))), # Keep this on the natural scale?
                         frog_abundance=max_double_na(adult, subadult)), by=.(site_id, year)][year >= 2004]

# Cast into wide form for analysis
obs_bd_mat = dcast(red_data, site_id ~ year,  value.var="bd_present", fill=-1)
sample_size_mat = dcast(red_data, site_id ~ year,  value.var="num_samps", fill=-1)
num_inf_mat = dcast(red_data, site_id ~ year,  value.var="num_inf", fill=-1)
bd_mat = dcast(red_data, site_id ~ year,  value.var="mean_bd", fill=-1)
abund_mat = dcast(red_data, site_id ~ year,  value.var="frog_abundance", fill=-1)

# Set up to view data as time-step transitions
sites = unique(red_data$site_id)
ipm_dat = list()
for(site in sites){
  
  tdat = red_data[site_id == site][order(year)]
  n = nrow(tdat)
  newdat = list()
  newdat$site_id = site
  newdat$year1 = tdat$year[1:(n - 1)]
  newdat$year2 = tdat$year[2:(n)]
  newdat$bd_occupied_past = tdat$bd_present[1:(n - 1)]
  newdat$bd_occupied_present = tdat$bd_present[2:(n)]
  newdat$mean_bd_past = tdat$mean_bd[1:(n - 1)]
  newdat$mean_bd_present = tdat$mean_bd[2:(n)]
  newdat$frog_abund_present = tdat$frog_abundance[2:(n)]
  newdat$num_inf_present = tdat$num_inf[2:n]
  ipm_dat[[site]] = as.data.frame(newdat)
  
}

ipm_dat = as.data.table(do.call(rbind, ipm_dat))


ipm_dat$diff = ipm_dat$year2 - ipm_dat$year1
anal_dat = ipm_dat[(diff == 1) & 
                   (bd_occupied_past == 1) & 
                   (bd_occupied_present == 1) & 
                   (!is.na(frog_abund_present)) & 
                   (frog_abund_present > 0)]

# Re-transform log loads
anal_dat[, mean_bd_past:=log10(exp(mean_bd_past) - 1)]
anal_dat[, mean_bd_present:=log10(exp(mean_bd_present) - 1)]

# Merge data with lake data
lake_dat = fread("../data/archival/sites_nov2020.csv")
colnames(lake_dat)[1] = "site_id"
anal_dat = merge(anal_dat, lake_dat[, .(site_id, area, elevation, jurisdiction)], key="site_id")
anal_dat[, density:=(frog_abund_present / sqrt(area / pi)*pi*2)]

# Delineate a unique identifier for region
anal_dat[, location:=(as.numeric(jurisdiction == "yosemite"))]
anal_dat_yose = anal_dat[jurisdiction == "yosemite"][order(site_id, year1)]
fwrite(anal_dat, "../data/candidate_datasets/analysisIII_loadanddensity.csv")

