Scripts used in the analysis for "Host density has limited effects on pathogen invasion, disease-induced declines, and within-host infection dynamics across a landscape of disease"

- `analysis_functions.R`:  Functions used primarily for preparing data and fitting the Hidden Markov Model to Bd state trajectories.

- `build_data_sets.R`: Script that applies specific criteria (described in the manuscript) to identify the frog populations included in each analysis

- `analysisI_failed_invasion.R`: Script for fitting HMM models to test whether host density affects failed Bd invasion

- `analysisII_decline.R`: Analysis of whether host density prior to decline is predictive of disease-induced declines.

- `analysisIII_infection_load.*`: Analysis of whether host density affects within-host infection intensity.

- `extract_snowdepth.R`: Extracts and formats CDEC snow-water equivalent data for each lake in the database from 2004-2019.

- `extract_temperature.R`: Extracts maximum and minimum daily temperatures for all lakes from 2004-2019.

- `make_candidate_failed_invasions.Rmd`:  Identifies sites that likely experiences failed invasions for human inspection.

- `test_fourstate_HMM.R`:  Tests that the HMM model with and without observation error can recover known quantities.

- `runall.R`: Specifies the order that R scripts need to to be run in order to reproduce the full analysis.
