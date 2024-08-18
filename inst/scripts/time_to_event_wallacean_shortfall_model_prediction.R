#' Predicting the number of species remaining to be geolocated in a particular region using STAN.
#'
#' The script allows to predict the number of species remaining to be geolocated in a given region and the median probability of geolocation
#' for an "average" species in the region
#'
#' EXPECTED INPUTS:
#'  - `stan_model`: path to the (.rds) file containing the `stan_surv` model object of the time-to-geolocation model
#'  - `no_geoloc`: number of species in the region (botanical country) without valid geographic coordinates
#'  - `nsamp`: number of samples to draw from the posterior predictive
#'
#' EXAMPLE OF COMMAND LINE:
#'  Rscript time_to_event_wallacean_shortfall_model_prediction.R --stan_model=modelObjects/FRA.rds --no_geoloc=1000 --nsamp=2000
#'
#' EXAMPLE OF SOURCE COMMAND:
#'  Fill in `time_to_event_wallacean_shortfall_model_prediction.R` with required information and parameters
#'  source("time_to_event_wallacean_shortfall_model_prediction.R")
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

# Number of species in the region (botanical country) without valid geographic coordinates
NUM_SPECIES_NOT_GEOLOCATED = 0 # should be changed to a number > 0, otherwise the prediction will be equal to zero

# Path to the directory where to save model prediction file (.csv)
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
    nsamp = 2000,
    no_geoloc = 0
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)

  STAN_MODEL <- args$stan_model
  if(is.null(STAN_MODEL)) cli::cli_abort("the path to the {.code stan_surv} object is missing !")
  NUM_SPECIES_NOT_GEOLOCATED <- args$no_geoloc
  OUTPUTDIR <- dirname(dirname(STAN_MODEL))
  NSAMP <- args$nsamp

}

if(!nzchar(STAN_MODEL)) cli::cli_abort("The {.code stan_model} file is missing !")
if(!file.exists(STAN_MODEL)) cli::cli_abort("Unable to find file: {.file {STAN_MODEL}}")
if(!dir.exists(OUTPUTDIR)) cli::cli_abort("Unable to find output directory: {.path {OUTPUTDIR}}")

#====================================================================
#== 1. Load time-to-geolocation model
#====================================================================

  # Get the TDWG level3 geographic code from the file name
  tdwg_code <- stringr::str_extract(basename(STAN_MODEL),".+?(?=\\.rds$)")
  cli::cli_h2("Load time-to-geolocation model for the geographic region: {tdwg_code}")

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

  # make prediction for species not geolocated yet
  nd = stan_fit$data %>%
    dplyr::filter(event==0)

  # if none is found in the data re-use the whole dataset
  if(nrow(nd)==0L){
    cli::cli_alert_warning("found no species remaining to be geolocated in the training data")
    nd = stan_fit$data
  }

  #====================================================================
  #== 2. Sample the probability of remaining non geolocated
  #====================================================================

  cli::cli_h2("Sample the probability of remaining non geolocated in {tdwg_code}")

  ps_surv <- stan_fit %>%
    rstanarm::posterior_survfit(type="surv",            # => survival probability
                                newdata = nd,           # prediction data
                                times="time",           # duration of the description process
                                last_time="time",       # last known or censoring time of species not yet described
                                extrapolate=FALSE,      # do not extrapolate
                                return_matrix=TRUE,     # return matrix of the samples (i.e. do not aggregate results)
                                draws=NSAMP, seed=SEED)

  # compute median probabilities +/- 95% credible interval (i.e. median of standardised predictions)
  cli::cli_progress_step("compute median probabilities +/- 95% credible interval")
  median_prob <- ps_surv %>%
    `[[`(1) %>%
    {
      data.frame(prob_nogeoloc=quantile(rowMeans(.,na.rm=T), probs=0.5, na.rm=T)['50%'],
                 ci_ll=quantile(rowMeans(.,na.rm=T), probs=0.025,na.rm=T)['2.5%'],
                 ci_ul=quantile(rowMeans(.,na.rm=T), probs=0.975,na.rm=T)['97.5%'], row.names=NULL)
    }

  dn <- file.path(OUTPUTDIR,"modelPredictions")
  if(!dir.exists(dn)) dir.create(dn, recursive=TRUE)

  # combine predictions
  cli::cli_progress_step("combine predictions")

  out = data.frame(LEVEL3_COD=tdwg_code,
                   SR=NUM_SPECIES_NOT_GEOLOCATED,
                   median_prob,
                   # compute number of species not geolocated yet
                   SRnoloc=NUM_SPECIES_NOT_GEOLOCATED * median_prob$prob_nogeoloc,
                   # compute standard deviation
                   SRnoloc_sd = sd(NUM_SPECIES_NOT_GEOLOCATED * rowMeans(ps_surv[[1]],na.rm=T)),
                   # compute lower interval bound
                   SRnoloc_ll= NUM_SPECIES_NOT_GEOLOCATED * median_prob$ci_ll,
                   # compute upper interval bound
                   SRnoloc_ul= NUM_SPECIES_NOT_GEOLOCATED * median_prob$ci_ul)

  # save the results
  cli::cli_progress_step("save results in {.path {dn}}")

	write.csv(out, file.path(dn,paste0(tdwg_code[1],".csv")), row.names=FALSE)

	cli::cli_progress_done()
