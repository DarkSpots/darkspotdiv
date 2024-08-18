#' Assessing the model fit of the time to species description in a particular region using the (integrated) Brier score.
#'
#' The script allows to assess the quality of the fit of a time to species description model fitted with an assemblage of species occurring
#' in a given geographic region and using a time-dependent estimate of the Brier score
#'
#' EXPECTED INPUTS:
#'  - `stan_model`: path to the (.rds) file containing the `stan_surv` model object of the time-to-description model
#'  - `nsamp`: number of samples to draw from the posterior predictive
#'
#' EXAMPLE OF COMMAND LINE:
#'  Rscript time_to_event_wallacean_shortfall_model_assessment.R --stan_model=modelObjects/FRA.rds --nsamp=100
#'
#' EXAMPLE OF SOURCE COMMAND:
#'  Fill in `time_to_event_wallacean_shortfall_model_assessment.R` with required information and parameters
#'  source("time_to_event_wallacean_shortfall_model_assessment.R")
#'
#'
#'

suppressWarnings(
suppressMessages({
library(dplyr)
library(magrittr)
library(data.table)
library(rstanarm)
library(pec)
library(survival)
})
)

# function that draws survival probabilities from a Stan model fit object and returns a matrix of survival probabilities as a function of time
surv_mat <- function(fit, ndraws=500) {

  ps <- rstanarm::posterior_survfit(fit,
                          times = 0,
                          extrapolate = TRUE,
                          control = list(epoints = 180),
                          return_matrix=TRUE,
                          draws = ndraws)

  ps_mat <- vector(mode="list",length=ndraws)
  ps_mat <- ps %>%
    purrr::map(function(x) data.frame(iterations=1:nrow(x), times=unique(attributes(x)$times), x, row.names=NULL)) %>%
    data.table::rbindlist() %>%
    split(by='iterations') %>%
    purrr::map(function(x) x %>% dplyr::select(-c(1:2)) %>% t())

  attr(ps_mat, "times") <- ps %>% purrr::map_dbl(function(x) unique(attributes(x)$times))
  names(ps_mat) <- paste0("iterations",1:ndraws)
  return(ps_mat)
}

# File path to the stan_surv model object file (.rds)
STAN_MODEL = ""

# Path to the directory/folder where to save the evaluations
OUTPUTDIR = dirname(dirname(STAN_MODEL))

# Draw NSAMP samples from the posterior predictive
NSAMP = 500

# check command line arguments if any
if (sys.nframe() == 0L) {

  #--------------------------------------------------------
  #= Read command line argument and check input parameters
  #--------------------------------------------------------
  cli::cli_alert_info("Reading command line arguments and checking parameters")

  # default arguments. Change if needed.
  default_args <- list(
    nsamp = 500
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)

  STAN_MODEL <- args$stan_model
  if(is.null(STAN_MODEL)) cli::cli_abort("The {.code stan_model} file is missing !")
  OUTPUTDIR <- dirname(dirname(STAN_MODEL))
  NSAMP <- args$nsamp
}

if(!nzchar(STAN_MODEL)) cli::cli_abort("The {.code stan_model} file is missing !")
if(!file.exists(STAN_MODEL)) cli::cli_abort("Unable to find file: {.file {STAN_MODEL}}")
if(!dir.exists(OUTPUTDIR)) cli::cli_abort("Unable to find output directory: {.path {OUTPUTDIR}}")

#====================================================================
#== 1. Get the time-to-event model fit
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

	# compute survival probabilities
  cli::cli_progress_step("compute survival probabilities from {NSAMP} posterior samples")
	surv_mat_obj <- surv_mat(stan_fit, ndraws=NSAMP)

	# release memory
	gc()
#====================================================================
#== 2. Compute Brier score
#====================================================================

	# extract the range of times used to compute the survival probability (i.e. probability of remaining undescribed)
	times_pec <- attr(surv_mat_obj, "times")

	# convert family random effect into numeric factor
	fam_num <- match(stan_fit$data$family, stan_fit$xlevs[["family"]])
	fam_num[is.na(fam_num)] <- as.numeric(as.factor(stan_fit$data$family[is.na(fam_num)])) + dplyr::n_distinct(stan_fit$xlevs[["family"]])
	stan_fit$data$family <- fam_num

	# extract model formula
	form1 <- update.formula(stan_fit$formula$fe_form, Surv(time, event) ~ .)

	# compute the prediction error curves (pec) aka the time-dependent Brier scores
	cli::cli_progress_step("compute the prediction error curves i.e. the {.strong time-dependent Brier score}")
	pec_rslt <- pec::pec(surv_mat_obj,
					formula = form1,
					data = stan_fit$data,
					times = times_pec,
					reference = FALSE,
					exact = FALSE,
					cens.model="none",  # no model censoring logic becaseu there is no censored observation.
					splitMethod = "none")

	# compute the cumulative prediction error curves (crps), aka integrated Brier scores.
	cli::cli_progress_step("compute the cumulative prediction error curves (crps), i.e. {.strong integrated Brier scores}")
	ibs <- pec::crps(pec_rslt, start=0, times=range(times_pec)[2], what="AppErr")

	# compute the average integrated Brier score +/- 95% interval
	cli::cli_progress_step("compute the average integrated Brier score +/- 95% interval")
	ibs.out <- data.frame(LEVEL3_COD=tdwg_code[1],
	                      bs=mean(ibs),
	                      bs_ci_ll = unname(quantile(ibs, probs=0.025)['2.5%']),
	                      bs_ci_ul=unname(quantile(ibs, probs=0.975)['97.5%']))

	# extract and format Brier scores for plotting
	cli::cli_progress_step("extract and format Brier scores for plotting")

	res <- pec_rslt$AppErr %>%
	  purrr::map_dfc(function(x) x) %>%
	  data.table::setDT()

	# mean Brier scores by time
	mbs <- res[, list(bs=mean(unlist(.SD)),
					  bs_ci_ll=quantile(unlist(.SD), probs=0.025)['2.5%'],
					  bs_ci_ul=quantile(unlist(.SD), probs=0.975)['97.5%']), by=seq_len(nrow(res))]
	mbs[, times:=times_pec]

#====================================================================
#== 4. Save the model outputs
#====================================================================

	dn <- file.path(OUTPUTDIR,"modelEvaluations","testing","Brier_score")
	if(!dir.exists(dn)) dir.create(dn, recursive=TRUE)

	cli::cli_progress_step("save results in {.path {dn}}")

	write.csv(ibs.out, file.path(dn,paste0(tdwg_code[1],".csv")), row.names=FALSE)

#====================================================================
#== 5. Plot Brier scores
#====================================================================

	Bstplt <- ggplot2::ggplot(data=mbs, ggplot2::aes(x=times, y=bs)) +
	  ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin = bs_ci_ll, ymax = bs_ci_ul),
	                       fill = "grey70", alpha=0.7) +
	  ggplot2::geom_step(color="black",
	                     linewidth=1,
	                     mapping=ggplot2::aes(x=times, y=bs),
	                     data=mbs, inherit.aes = FALSE) +
	  ggplot2::theme_classic() +
	  ggplot2::labs(title="Time-dependent model assessment",
					subtitle=paste0("Brier score calculated at WGSRPD-level3 code: ", tdwg_code[1]))+
	  ggplot2::theme(
		text=ggplot2::element_text(family="serif"),
		plot.title = ggplot2::element_text(face="bold", size=14, hjust = 0),
		axis.text=ggplot2::element_text(colour="black", size=12),
		axis.title=ggplot2::element_text(colour="black", size=16,face="bold")
	  ) +
	  ggplot2::xlab(expression('Years from 1754')) + ggplot2::ylab("Brier score") +
	  ggplot2::scale_x_continuous(expand=c(0,0)) +
	  ggplot2::scale_y_continuous(expand=c(0,0), limits=c(0,1))

	dn1 <- file.path(dirname(dn),"plots")
	if(!dir.exists(dn1)) dir.create(dn1, recursive=TRUE)

	cli::cli_progress_step("plot Brier scores and save results in {.path {dn1}}")

	pdf(file=file.path(dn1, paste0(tdwg_code[1],'.pdf')))
	# a main (background) plot
	print(Bstplt);
	dev.off()

	cli::cli_progress_done()
