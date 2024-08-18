#' Modelling the time to species description using STAN.
#'
#' The script allows to:
#' * fit a time-to-event model using the time to the description of an assemblage of species within a given region
#'  as a function of species attributes/characteristics.
#' * assess the model fit using the Brier score
#' * predict the number of species remaining to be described in the region
#'
#' EXPECTED INPUTS:
#'  - `occ_data`: path to the file (.parquet) containing the time-to-event data for species in a given region with valid occurrence records
#'  - `gap_data`: path to the file (.parquet) containing the time-to-event data for species in the same region as `occ_data`, but with no
#'   valid occurrence records
#'  - `output_dir`: path to a directory to save modelling outputs. Default is current directory.
#'  - `end_date`: date at which to make inference about the time to species description.
#'  - `evaluate`: should the model calibration be evaluated ? Default is TRUE.
#'  - `predict`: should the model make (new) predictions ? Default is TRUE.
#'
#' EXAMPLE OF COMMAND LINE:
#'  Rscript time_to_event_linnaean_shortfall_model_fitting.R --occ_data=~/data/FRA.parquet --gap_data=~/gap/FRA.parquet --evaluate=FALSE --predict=TRUE
#'
#' EXAMPLE OF SOURCE COMMAND:
#'  Fill in `time_to_event_linnaean_shortfall_model_fitting.R` with required information and parameters
#'  source("time_to_event_linnaean_shortfall_model_fitting.R")
#'
#'
#'
suppressWarnings(
suppressMessages({
library(arrow)
library(dplyr)
library(magrittr)
library(data.table)
library(rstanarm)
})
)

# Path to the file (.parquet) containing the time-to-event data for species with valid occurrence records
OCC_DATA = ""

# Path to the file (.parquet) containing the time-to-event data for species with no valid occurrence records
GAP_DATA = ""

# Path to the directory/folder where to save all modelling outputs
OUTPUT_DIR = getwd()

# Date at which to make inference about the time to species description
END_DATE = 2021

# Should the model calibration be evaluated ?
EVALUATE = TRUE

# Should the model make (new) predictions ?
PREDICT = TRUE

# check command line arguments if any
if (sys.nframe() == 0L) {

  #--------------------------------------------------------
  #= Read command line argument and check input parameters
  #--------------------------------------------------------
  cli::cli_alert_info("Reading command line arguments and checking parameters")

  # default arguments. Change if needed.
  default_args <- list(
    output_dir = getwd(),
    end_date = 2021,
    evaluate = TRUE,
    predict = TRUE
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)

  OCC_DATA <- args$occ_data
  GAP_DATA <- args$gap_data
  OUTPUT_DIR <- args$output_dir
  END_DATE <- args$end_date
  EVALUATE <- as.logical(args$evaluate)
  PREDICT <- as.logical(args$predict)
}


if(!file.exists(OCC_DATA)) cli::cli_abort("Unable to find {.code OCC_DATA} file: {.file {OCC_DATA}}")
if(!file.exists(GAP_DATA)) cli::cli_abort("Unable to find {.code GAP_DATA} file: {.file {GAP_DATA}}")
if(!dir.exists(OUTPUT_DIR)) cli::cli_abort("Unable to find {.code OUTPUT_DIR} directory: {.file {OUTPUT_DIR}}")


#====================================================================
#== 1. Build time-to-event data
#====================================================================
  cli::cli_h1("Building time-to-description data")

	# Get the TDWG level3 geographic code from the file name
	tdwg_code <- stringr::str_extract(basename(OCC_DATA),".+?(?=\\.)")

	# time-to-event data for species in tdwg_code with occurrence records
	cli::cli_progress_step("read in data for species with occurrence records in botanical country: {tdwg_code}")
	time_to_event <- OCC_DATA %>%
	 arrow::read_parquet(as_data_frame = FALSE) %>%
	  # keep species described after the starting date (1753) up to END_DATE
	 dplyr::filter(time_todescription > 0, time_ofdescription <= END_DATE) %>%
	 dplyr::collect() %>%
	 data.table::setDT()

	# time-to-event data for species in `tdwg_code` with no occurrence records
	cli::cli_progress_step("read in data for species with no occurrence records in {tdwg_code}")
	time_to_eventt <- arrow::read_parquet(GAP_DATA, as_data_frame = FALSE) %>%
	  # keep species described after 1753 and before END_DATE
	 dplyr::filter(time_todescription > 0, time_ofdescription <= END_DATE) %>%
	 dplyr::collect() %>%
	 data.table::setDT() %>%
	 data.table::setcolorder(names(time_to_event))

   time_to_event %<>%
	 dplyr::bind_rows(time_to_eventt)

	# pull description times
	time_obs <- time_to_event %>%
	 dplyr::pull(time_todescription)

	# explanatory variables
	cli::cli_progress_step("select and scale continuous predictors")
	Xdata <- time_to_event %>%
	 dplyr::select(CHELSA_bio_1:num_uses,
				   lifeform_description,
				   family) %>%
	  # scale the continuous variables (substract their mean and divide by their standard deviation)
	 dplyr::mutate_if(is.numeric, function(x) as.vector(scale(x)))

	# set the number of uses for species with no known use to 0
	if(all(is.na(Xdata$num_uses))){
	 Xdata[,num_uses:=0]
	}

	# create the time-to-event dataset for stan_surv
	cli::cli_progress_step("create dataset for {.code stan_surv} function")
	stan_survdata = data.frame(time=time_obs, event=1, Xdata) %>%
        dplyr::mutate(lifeform_description=as.factor(lifeform_description))

	# define model formula
	EXCLUDE_VARS = "CHELSA_bio_1" # exclude CHELSA_bio_1 for collinearity issue (set to "" if unnecessary)

	if(length(EXCLUDE_VARS)>0L){
	  cli::cli_alert_info("exclude {.vars {EXCLUDE_VARS}} from the list of predictors")
	  stan_survdata %<>%
	    dplyr::select(-dplyr::all_of(EXCLUDE_VARS))
	}
	form0 <- formula(stats::model.frame(survival::Surv(time, event)~ ., data=stan_survdata))
	form1 <- update(form0, ~. - family -1 + (1|family))

#---------------------------------
#== 1.1 Detect time-varying effects
#---------------------------------
	cli::cli_progress_step("detect time-varying effects")

	# update model formula for Cox model
	form1b <- update.formula(form1, ~ . -(1|family) + cluster(family))

	# get model data
	stan_fit_data <- stan_survdata

	# extract time-varying terms
	tve_terms = tryCatch({
  	# fit a Cox model
  	survival_fit <- survival::coxph(form1b, data=stan_fit_data)

  	# test proportional hazard assumption
  	test_PH <- survival::cox.zph(survival_fit,terms=FALSE,global=FALSE)

  	# detect terms with time-varying effects (tve) i.e. non-proportional hazards
  	test_PH[[1]] %>%
  	as.data.frame() %>%
  	tibble::rownames_to_column(var="Variables") %>%
  	dplyr::filter(p<0.05, !grepl("lifeform_description",Variables)) %>% dplyr::pull(Variables)
  		},error=function(err){
  		setdiff(attr(terms(form1),"term.labels"),c("lifeform_description","1 | family"))
	})

#---------------------------------------------------------
#= 1.2 Update model formula to include time-varying effect
#---------------------------------------------------------

	if(length(tve_terms)>0L){
	  cli::cli_alert_info("detected {.strong {length(tve_terms)}} time-varying effect{?s} for {.vars tve_terms}")
	  cli::cli_progress_step("update model formula")
		tfe_terms <- setdiff(attr(terms(form1),"term.labels"),c(tve_terms,"1 | family","lifeform_description"))
    formula_tfe_part <- paste(tfe_terms, collapse=" + ")
		formula_tve_part <- paste(paste0("tve(",tve_terms,")"), collapse = " + ")
		formula_rest <- paste(tve_terms, collapse="-")
		formula_tve <- as.formula(paste("~",formula_tve_part,ifelse(nzchar(formula_tfe_part)," + ",""),formula_tfe_part,"-",formula_rest,"+ lifeform_description -1 + (1|family)"))

		form1 <- update.formula(form1,formula_tve)
	}else{
	  cli::cli_alert_warning(" No time-varying effect detected")
	}

#====================================================================
#== 2. Run the time-to-event model
#====================================================================

	cli::cli_h1("Building time-to-description model")

	#-----------------------------------
	#== 2.1 Set up modelling parameters
	#-----------------------------------

	cli::cli_progress_step("set up modelling parameters")
	# MCMC stan setup parameters
	NCHAINS = 2
	NCORES = parallel::detectCores() - 1
	NITER = 3000
	WARMUP = 1000
	SEED = 1234
	options(mc.cores = NCORES)
	# model baseline (+options)
	BASEHAZ = "ms" # m-splines
	BASEHAZ_OPS = list(degree=3, df=10) # cubic m-spline model (with df = 10)

	cli::cli_bullets(c("i"="Monte-Carlo Markov Chain (MCMC) parameters:",
	                   "*"="{.strong {NCHAINS}} chain{?s}",
	                   "*"="{.strong {NITER}} iterations",
	                   "*"="{.strong {WARMUP}} burn-in iteration{?s}")
	)
	cli::cli_alert_info("Baseline hazard function: {.strong {BASEHAZ} model}")

	cli::cli_progress_step("run {.code stan_surv}")
	stan_fit <- rstanarm::stan_surv(form1,
								   data = stan_survdata,
								   basehaz =  BASEHAZ,
								   basehaz_ops = BASEHAZ_OPS,
								   chains = NCHAINS,
								   cores = NCORES,
								   iter = NITER,
                   warmup = WARMUP,
								   seed = SEED)

#====================================================================
#== 3. Save the model outputs
#====================================================================
	dn <- file.path(OUTPUT_DIR,"time_to_event","stan_models","Linnaean","modelObjects")
	if(!dir.exists(dn)) dir.create(dn, recursive=TRUE)
	cli::cli_progress_step("save the {.code stan_surv} model object in {.path {dn}}")
	saveRDS(stan_fit, file.path(dn,paste0(tdwg_code[1],".rds")))
	cli::cli_progress_done()

#====================================================================
#== 4. Run model assessment script
#====================================================================
	if(EVALUATE){
	  cli::cli_h1("Evaluating model performance")
	  model_assessment_script = system.file("scripts","time_to_event_linnaean_shortfall_model_assessment.R", package="darkspotdiv")
  	cmd <- sprintf("%s %s --stan_model=%s",
  	               Sys.which("Rscript"),
  	               model_assessment_script,
  	               file.path(dn,paste0(tdwg_code[1],".rds")))
  	system(cmd)
	}

#====================================================================
#== 5. Run model prediction script
#====================================================================
	if(PREDICT){
	  cli::cli_h1("Making model prediction")
	  model_prediction_script = system.file("scripts","time_to_event_linnaean_shortfall_model_prediction.R", package="darkspotdiv")
  	cmd <- sprintf("%s %s --stan_model=%s --no_geoloc=%g",
  	               Sys.which("Rscript"),
  	               model_prediction_script,
  	               file.path(dn,paste0(tdwg_code[1],".rds")),
  	               nrow(time_to_eventt))
  	system(cmd)
	}

	cli::cli_progress_done()
