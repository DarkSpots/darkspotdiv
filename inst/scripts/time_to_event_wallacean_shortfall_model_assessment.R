#' Assessing the model fit of the time to species geolocation in a particular region using the (integrated) Brier score.
#'
#' The script allows to assess the quality of the fit of a time to species geolocation model fitted with an assemblage of species occurring
#' in a given geographic region and using a time-dependent estimate of the Brier score
#'
#' EXPECTED INPUTS:
#'  - `stan_model`: path to the (.rds) file containing the `stan_surv` model object of the time-to-geolocation model
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
surv_mat <- function(fit, newdata_=NULL, times_=0, extrapolate_=FALSE, ndraws=100) {

  ps <- rstanarm::posterior_survfit(fit,
                                    newdata = newdata_,
                                    times = times_,
                                    extrapolate = extrapolate_,
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

# Path to the directory where to save model evaluations
OUTPUTDIR = dirname(dirname(STAN_MODEL))

# Draw NSAMP samples from the posterior predictive
NSAMP = 100

# check command line arguments if any
if (sys.nframe() == 0L) {

  #--------------------------------------------------------
  #= Read command line argument and check input parameters
  #--------------------------------------------------------
  cli::cli_alert_info("Reading command line arguments and checking parameters")

  # default arguments. Change if needed.
  default_args <- list(
    nsamp = 100
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

cli::cli_h2("Load time-to-geolocation model for the geographic region: {tdwg_code}")

# load the model
stan_fit <- readRDS(STAN_MODEL)

#====================================================================
#== 2. Compute Brier score
#====================================================================

if(stan_fit$nlcens>0L){  # if the model contains left-censored observations
  cli::cli_progress_step(" the {.code stan_model} contains {.strong {stan_fit$nlcens}} left-censored observations")

  # convert family random effect into numeric factor
  fam_num <- match(stan_fit$data$family, stan_fit$xlevs[["family"]])
  fam_num[is.na(fam_num)] <- as.numeric(as.factor(stan_fit$data$family[is.na(fam_num)])) + dplyr::n_distinct(stan_fit$xlevs[["family"]])
  stan_fit$data$family <- fam_num

  # convert the training data to interval censored data
  cli::cli_progress_step("convert the training data to interval censored data")
  stan_survdata <- stan_fit$data %>%
    dplyr::mutate(time2=dplyr::if_else(event==0,Inf,time2),
                  time=dplyr::if_else(event==2,0,time))

  # create a `Surv` object
  cli::cli_progress_step("create a {.code Surv} object with newly created censored intervals")
  sobject <- survival::Surv(stan_survdata$time, stan_survdata$time2, type="interval2")
  intL <- sobject[, 1]
  intR <- sobject[, 2]
  times_ <- c(intL,intR)

  # compute survival probabilities at each interval limit for each species
  cli::cli_progress_step("compute survival probabilities at each interval limit for each species from {NSAMP} posterior samples")
  newdata = data.frame(time=times_,rbind(stan_survdata[,-c(1:3)], stan_survdata[,-c(1:3)]))
  surv_mat <- rstanarm::posterior_survfit(stan_fit,
                                          newdata = newdata,
                                          times = 'time',
                                          extrapolate = FALSE,
                                          return_matrix=TRUE,
                                          draws = NSAMP)

  require(doParallel)
  cli::cli_progress_step("compute the average integrated Brier score +/- 95% interval", spinner=TRUE)
  ibs <- foreach::foreach(j=iterators::iter(surv_mat[[1]], by='row'), .combine="c") %do% {
    cli::cli_progress_update()
    survs_m <- matrix(rep(as.vector(j), nrow(sobject)), nrow = length(times_))
    ICcforest::sbrier_IC(sobject, survs_m, type="IBS")[1] |>  unname()
  }

  ibs.out <- data.frame(LEVEL3_COD=tdwg_code[1], bs=mean(ibs), bs_ci_ll = unname(quantile(ibs, probs=0.025)['2.5%']), bs_ci_ul=unname(quantile(ibs, probs=0.975)['97.5%']))

}else{ # If the model contains right-censored of no censored observations

  cli::cli_progress_step(" the {.code stan_model} contains {.strong {ifelse(stan_fit$nrcens>0L,stan_fit$nrcens,'no')}} {ifelse(stan_fit$nrcens>0L,'left-',''}censored observation{?s}")

  # compute survival probabilities
  cli::cli_progress_step("compute survival probabilities from {NSAMP} posterior samples")
  surv_mat_obj <- surv_mat(stan_fit, extrapolate_=TRUE, ndraws=NSAMP)

  # release memory
  gc()

  times_pec <- attr(surv_mat_obj, "times")

  # convert family random effect into numeric factor
  fam_num <- match(stan_fit$data$family, stan_fit$xlevs[["family"]])
  fam_num[is.na(fam_num)] <- as.numeric(as.factor(stan_fit$data$family[is.na(fam_num)])) + dplyr::n_distinct(stan_fit$xlevs[["family"]])
  stan_fit$data$family <- fam_num

  # extract model formula
  form1 <- update.formula(stan_fit$formula$formula, Surv(time, event) ~ .)
  cli::cli_progress_step("compute the prediction error curves i.e. the {.strong time-dependent Brier score}")

  pec_rslt <- pec::pec(surv_mat_obj,
                       formula = form1,
                       data = stan_fit$data,
                       times = times_pec,
                       reference = FALSE,
                       exact = FALSE,
                       cens.model="marginal",
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


  #-----------------------
  #=2.1 Plot Brier scores
  #-----------------------

  Bstplt <- ggplot2::ggplot(data=mbs, ggplot2::aes(x=times, y=bs)) +
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin = bs_ci_ll, ymax = bs_ci_ul), fill = "grey70", alpha=0.7) +
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
    ggplot2::xlab(expression('Years from description')) + ggplot2::ylab("Brier score") +
    ggplot2::scale_x_continuous(expand=c(0,0)) +
    ggplot2::scale_y_continuous(expand=c(0,0), limits=c(0,1))

  dn1 <- file.path(dirname(dn),"plots")
  if(!dir.exists(dn1)) dir.create(dn1, recursive=TRUE)

  cli::cli_progress_step("plot Brier scores and save results in {.path {dn1}}")

  pdf(file=file.path(dn1, paste0(tdwg_code[1],'.pdf')))
  # a main (background) plot
  print(Bstplt);
  dev.off()
}

  #====================================================================
  #== 3. Save the model outputs
  #====================================================================

  dn <- file.path(OUTPUTDIR,"modelEvaluations","testing","Brier_score")
  if(!dir.exists(dn)) dir.create(dn, recursive=TRUE)

  cli::cli_progress_step("save results in {.path {dn}}")

  write.csv(ibs.out, file.path(dn,paste0(tdwg_code[1],".csv")), row.names=FALSE)

  cli::cli_progress_done()

