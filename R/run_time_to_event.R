#'                                 Run a time-to-event model
#'
#' Utility function to run a time-to-event model
#'
#' @param occ_data A path to the file (.parquet) containing the time-to-event data for species in a given region with valid occurrence records.
#' @param gap_data A path to the file (.parquet) containing the time-to-event data for species in the same region as `occ_data`, but with no
#'   valid occurrence record
#' @param outputdir A path to a directory to save modelling outputs. Default is current directory
#' @param event A string specifying the event to model. It must be either `description` or `geolocation`.
#' @param end_date A numeric specifying the year at which to make inference about the time to species description.
#' @param evaluate A logical. Should the model calibration be evaluated ? Default is `TRUE`
#' @param predict A logical. Should the model make (new) predictions ?  Default is `TRUE`.
#' @return An integer number. `0` in case of success otherwise `-1`
#' @export
run_time_to_event_model <- function(occ_data,
                                    gap_data,
                                    outputdir = NULL,
                                    event = "description",
                                    end_date=NULL,
                                    evaluate=TRUE,
                                    predict=TRUE){

  # type of event to model
  event = match.arg(event, choices=c("description","geolocation"))

  # find the script to run the model fitting given the event
  event_model <- switch(event,
         description = system.file("scripts","time_to_event_linnaean_shortfall_model_fitting.R", package="darkspotdiv"),
         geolocation = system.file("scripts","time_to_event_wallacean_shortfall_model_fitting.R", package="darkspotdiv")
         )

  # find the Rscript executable file
  cmd <- Sys.which("Rscript")
  if(!nzchar(cmd)) stop("Cannot locate Rscript.exe")

  # create the base commande line
  cmd_args <- sprintf('%s --occ_data=%s --gap_data=%s --evaluate=%s --predict=%s', event_model, occ_data, gap_data, evaluate, predict)

  # update the command line with optional arguments
  if(!is.null(end_date)){
    cmd_args <- paste0(cmd_args,' --end_date=',end_date)
  }
  if(!is.null(outputdir)){
    cmd_args <- paste0(cmd_args,' --output_dir=',outputdir)
  }

  # run the command line
  tryCatch({
    out <- suppressWarnings(system2(cmd, cmd_args, stdout=TRUE, stderr=TRUE))
    status <- attr(out,"status")
    # check for error
    if(!is.null(status)){
      err <- attr(out,"errmsg")
      if(!is.null(err)){
        stop(err)
      }
      if(status==1) return(-1)
    }
    return(0)
  }, error=function(err) {message(err); return(-1);})

  on.exit({cli::cli_progress_done()})

}


#'                                 Make time-to-event predictions
#'
#' Utility function to make prediction from a time-to-event model
#'
#' @param stan_model A path to the file (.parquet) containing the time-to-event data for species in a given region with valid occurrence records.
#' @param start_date A numeric specifying the year of reference the time to species description follows from. Default is 1754. Only used if `event='description'`.
#' @param end_date A numeric specifying the year at which to make inference about the time to species description. Only used if `event='description'`.
#' @param no_geoloc A numeric specifying the number of species in the region (e.g. botanical country) without valid geographic coordinates. Only used if `event='geolocation'`.
#' @param event A string specifying the event to model. It must be either `description` or `geolocation`.
#' @param outputdir A path to a directory to save the predictions. Default is parent of the current directory
#' @return An integer number. `0` in case of success otherwise `-1`
#' @export
predict_time_to_event_model <- function(stan_model,
                                        start_date=1754,
                                        end_date=2021,
                                        no_geoloc = NULL,
                                        outputdir = NULL,
                                        event = "description",
                                        nsamp=NULL){

  # type of event to model
  event = match.arg(event, choices=c("description","geolocation"))

  if(event=="geolocation" & is.null(no_geoloc)){
    stop("For the 'time-to-geolocation' model the number of species missing geolocation must be provided.")
  }
  # find the script to make predictions given the event
  event_model <- switch(event,
                        description = system.file("scripts","time_to_event_linnaean_shortfall_model_prediction.R", package="darkspotdiv"),
                        geolocation = system.file("scripts","time_to_event_wallacean_shortfall_model_prediction.R", package="darkspotdiv")
  )

  # find the Rscript executable file
  cmd <- Sys.which("Rscript")
  if(!nzchar(cmd)) stop("Cannot locate Rscript.exe")

  # create the base commande line
  cmd_args <- switch(event,
                     description = sprintf('%s --stan_model=%s --start_date=%s --end_date=%s', event_model, stan_model, start_date, end_date),
                     geolocation = sprintf('%s --stan_model=%s --no_geoloc=%s', event_model, stan_model, no_geoloc)
                     )

  # update the command line with optional arguments
  if(!is.null(nsamp)){
    cmd_args <- paste0(cmd_args,' --nsamp=',nsamp)
  }
  if(!is.null(outputdir)){
    cmd_args <- paste0(cmd_args,' --output_dir=',outputdir)
  }

  # run the command line
  tryCatch({
    out <- suppressWarnings(system2(cmd, cmd_args, stdout=TRUE, stderr=TRUE))
    status <- attr(out,"status")
    # check for error
    if(!is.null(status)){
      err <- attr(out,"errmsg")
      if(!is.null(err)){
         stop(err)
      }
      if(status==1) return(-1)
    }
    return(0)
    }, error=function(err) {message(err); return(-1);})

}
