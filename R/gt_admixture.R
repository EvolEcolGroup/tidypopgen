#' Run ADMIXTURE from R
#'
#' This function runs ADMIXTURE, taking either a `gen_tibble` or a file as an input.
#'
#' @details This is a wrapper for the command line program ADMIXTURE. It can either
#' use a binary present in the main environment, or use a copy installed in a conda
#' environment. If the package `tidygenclust` is installed, and its conda environment
#' has been set up with `ADMIXTURE` in it (the default), it will automatically
#' use that version unless you change `conda_env` to "none".
#' @param x a `gen_tibble` or a character giving the path of the input file
#' @param k an integer giving the number of clusters
#' @param crossval boolean, should cross validation be used to assess the fit (defaults ot FALSE)
#' @param n_cors number of cores (defaults to 1)
#' @param conda_env the name of the conda environment to use. If set to "auto", the default,
#' a copy from `tidygenclust` will be preferred if available, otherwise a local copy will
#' be used. "none" forces the use of a local copy, whilst any other string will
#' direct the function to use a custom conda environment.
#' @return a list with the following elements:
#' - `Q` a matrix with the admixture proportions
#' - `loglik` the log likelihood of the model
#' - `cv` the cross validation error (if `crossval` is TRUE)
#' @export


gt_admixture <- function (x, k, crossval = FALSE, n_cores = 1, conda_env="auto"){
  # if x is a character, check that it is file that exists
  if (is.character(x)) {
    if (!file.exists(x)) {
      stop("The file ", x, " does not exist")
    }
    input_file <- x
  } else { # if x is a gen_tibble
    if (!is(x, "gen_tibble")) {
      stop("x must be a gen_tibble or a character")
    }
    # write the gen_tibble to a temp file
    input_file <- gt_as_plink(x, file = tempfile(fileext = "."))
  }

  # cast k as an integer
  k <- as.integer(k)

  # set the arguments for admixture
  admixture_args <- paste(input_file, k)
  if (crossval) {
    admixture_args <- paste(admixture_args, "--cv")
  }
  if (n_cores > 1) {
    admixture_args <- paste(admixture_args, paste0("-j", n_cores))
  }
  # set and check the conda environemnt
  # if there is no reticulate
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    # but we have not set conda_env to none
    if (conda_env!="none"){
      # if we have auto, set to none
      if (conda_env == "auto"){
        conda_env <- "none"
      } else {
        stop("The package reticulate is needed to use a conda environemnt in R.")
      }
    }
  } else { # if reticulate is available
    # if we have "auto" and the conda environment rtidygenclust does not exist, or
    if (conda_env =="auto") {
      if (("rtidygenclust" %in% reticulate::conda_list()[["name"]])){
        conda_env <- "rtidygenclust"
      } else {
        conda_env <- "none"
      }
      # check that the conda environment does exist
      if ((conda_env !="none") && (!conda_env %in% reticulate::conda_list()[["name"]])){
        stop("The conda environment ", conda_env, " does not exist.")
      }
    }
  }
  # if we conda_env is none
  if (conda_env == "none") {
    # use the system version of admixture (if it is installed)
    if (system2("which", args = "admixture") == 0) {
      system2("admixture", args = admixture_args)
    } else {
      stop("admixture is not installed in this path")
    }
  } else {
    reticulate::conda_run2("admixture", args = admixture_args,
                           conda = conda_env)
  }

  # read the output
  Q <- read.table(paste0(input_file, ".Q"), header = FALSE)
  loglik <- read.table(paste0(input_file, ".log"), header = FALSE)
  if (crossval) {
    cv <- read.table(paste0(input_file, ".cv"), header = FALSE)
    return(list(Q = Q, loglik = loglik, cv = cv))
  } else {
    return(list(Q = Q, loglik = loglik))
  }

}


