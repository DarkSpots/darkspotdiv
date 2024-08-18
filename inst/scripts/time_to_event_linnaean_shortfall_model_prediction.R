#' Predicting the number of species remaining to be described in a particular region using STAN.
#'
#' The script allows to predict the number of species remaining to be described in a given region and the median probability of description
#' for an "average" species in the region
#'
#' EXPECTED INPUTS:
#'  - `stan_model`: path to the (.rds) file containing the `stan_surv` model fit of the time-to-description model
#'  - `start_date`: date at which starts the description process (should be consistent with `stan_model`)
#'  - `end_date`: date at which to make inference about the time to species description.
#'  - `nsamp`: number of samples to draw from the posterior predictive
#'
#' EXAMPLE OF COMMAND LINE:
#'  Rscript time_to_event_linnaean_shortfall_model_prediction.R --stan_model=modelObjects/FRA.rds --start_date=1754 --end_date=2021 --nsamp=2000
#'
#' EXAMPLE OF SOURCE COMMAND:
#'  Fill in `time_to_event_linnaean_shortfall_model_prediction.R` with required information and parameters
#'  source("time_to_event_linnaean_shortfall_model_prediction.R")
#'
#'
#'
suppressWarnings(
suppressMessages({
library(dplyr)
library(magrittr)
library(data.table)
library(rstanarm)
})
)

# File path to the stan_surv model object file (.rds)
STAN_MODEL = ""

# Date at which to start the description process (should be the same as in the stan_model)
START_DATE = 1754

# Date at which to make the prediction about the time to species description
END_DATE = 2021

# Path to the directory/folder where to save the predictions
OUTPUTDIR = dirname(dirname(STAN_MODEL))

# Draw NSAMP samples from the posterior predictive
NSAMP = 2000
SEED = 1234

# check command line arguments if any
if (sys.nframe() == 0L) {

  #--------------------------------------------------------
  #= Read command line argument and check input parameters
  #--------------------------------------------------------
  cli::cli_alert_info("Reading command line arguments and checking parameters")

  # default arguments. Change if needed.
  default_args <- list(
    start_date = 1754,
    end_date = 2021
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)

  STAN_MODEL <- args$stan_model
  if(is.null(STAN_MODEL)) cli::cli_abort("The `stan_model` file is missing !")

  START_DATE <- args$start_date
  END_DATE <- args$end_date
  OUTPUTDIR <- dirname(dirname(STAN_MODEL))
  NSAMP <- args$nsamp
}

if(!file.exists(STAN_MODEL)) cli::cli_abort("Unable to find file: {.file {STAN_MODEL}}")
if(!dir.exists(OUTPUTDIR)) cli::cli_abort("Unable to find output directory: {.path {OUTPUTDIR}}")
if(END_DATE <= START_DATE) cli::cli_abort("END_DATE must be > START_DATE")

#====================================================================
#== 1. Load time-to-description model
#====================================================================

  # Get the TDWG level3 geographic code from the file name
	tdwg_code <- stringr::str_extract(basename(STAN_MODEL),".+?(?=\\.rds$)")
	cli::cli_h2("Load time-to-description model for the geographic region: {tdwg_code}")

	# load the model
	cli::cli_progress_step("read in {.file {STAN_MODEL}}")
	stan_fit <- readRDS(STAN_MODEL)

	# fix for 'predict' method
	if(stan_fit$basehaz$type_name=="ms"){
	  class(stan_fit$basehaz$basis) <- c("MSpline","splines2","matrix")
	  cli::cli_alert_success("fixed the {.code stan_surv} object to allow the {.code predict} method")
	}

	# convert family random effect into numeric factor
	fam_num <- match(stan_fit$data$family, stan_fit$xlevs[["family"]])
	fam_num[is.na(fam_num)] <- as.numeric(as.factor(stan_fit$data$family[is.na(fam_num)])) + dplyr::n_distinct(stan_fit$xlevs[["family"]])
	stan_fit$data$family <- fam_num

	# time parameter
	cli::cli_progress_step("compute time period for prediction")
	period_length = min(END_DATE - START_DATE, max(stan_fit$eventtime))

	#====================================================================
	#== 2. Sample the probability of remaining undescribed
	#====================================================================
	cli::cli_h2("Sample the probability of remaining undescribed in {tdwg_code}")

  ps_surv <- stan_fit %>%
    rstanarm::posterior_survfit(type="surv",            # => survival probability
                                newdata = stan_fit$data,# re-use calibration data
                                times=period_length,    # duration of the description process
                                last_time="time",       # last known or censoring time of species not yet described
                                extrapolate=FALSE,      # do not extrapolate beyond `period_length`
                                condition=FALSE,        # do not compute conditional survival probability
                                return_matrix=TRUE,     # return matrix of the samples (i.e. do not aggregate results)
                                draws=NSAMP, seed=SEED)

  # compute median probabilities +/- 95% credible interval (i.e. median of standardised predictions)
  cli::cli_progress_step("compute median probabilities +/- 95% credible interval")
  median_prob <- ps_surv %>%
    `[[`(1) %>%
  {
     data.frame(prob=quantile(rowMeans(.,na.rm=T), probs=0.5, na.rm=T)['50%'],
     ci_ll=quantile(rowMeans(.,na.rm=T), probs=0.025,na.rm=T)['2.5%'],
     ci_ul=quantile(rowMeans(.,na.rm=T), probs=0.975,na.rm=T)['97.5%'], row.names=NULL)
  }


	dn <- file.path(OUTPUTDIR,"modelPredictions")
	if(!dir.exists(dn)) dir.create(dn, recursive=TRUE)

	# get the species richness in the botanical country
	cli::cli_progress_step("get the observed species richness in the botanical country")

	bot_country_species_number = rWCVP::wcvp_checklist(taxon_rank="species",
	                                            area_codes = tdwg_code[1],
	                                            native=TRUE,
	                                            introduced=FALSE,
	                                            extinct=FALSE,
	                                            location_doubtful = FALSE,
	                                            hybrids=FALSE,
	                                            synonyms = FALSE) %>%
	  dplyr::distinct(plant_name_id) %>%
	  nrow()

	# combine predictions
	cli::cli_progress_step("combine predictions")

	out = data.frame(LEVEL3_COD=tdwg_code,
	                 SR=bot_country_species_number,
	                 median_prob,
	                 # compute number of species not described yet (1-median_prob$prob is the probability of being described)
	                 SRdesc=(bot_country_species_number/(1-median_prob$prob))-bot_country_species_number,
	                 # compute standard deviation
	                 SRdesc_sd = sd((bot_country_species_number/(1-rowMeans(ps_surv[[1]],na.rm=T)))-bot_country_species_number,na.rm=T),
	                 # compute lower interval bound
	                 SRdesc_ll=(bot_country_species_number/(1-median_prob$ci_ll))-bot_country_species_number,
	                 # compute upper interval bound
	                 SRdesc_ul=(bot_country_species_number/(1-median_prob$ci_ul))-bot_country_species_number)

	# save the results
	cli::cli_progress_step("save results in {.path {dn}}")

  write.csv(out, file.path(dn,paste0(tdwg_code[1],".csv")), row.names=FALSE)

  cli::cli_progress_done()
