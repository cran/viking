# Copyright 2019-2022 EDF, Sorbonne Université and CNRS.
# Author : Joseph de Vilmarest (EDF, Sorbonne Université)
# The package Viking is distributed under the terms of the license GNU LGPL 3.

#' Select time-invariant variances of a State-Space Model
#'
#' \code{select_Kalman_variances} is a function to choose hyper-parameters of the
#' linear Gaussian State-Space Model with time-invariant variances. It relies on the
#' functions \code{iterative_grid_search} and \code{expectation_maximization}.
#'
#' @param ssm the statespace object
#' @param X explanatory variables
#' @param y time series
#' @param method (optional, default \code{'igd'}) it can be either
#' \describe{
#' \item{\code{'igd'}}{\code{iterative_grid_search} is called}
#' \item{\code{'em'}}{\code{expectation_maximization} is called}
#' }
#' @param ... additional parameters
#'
#' @return a new statespace object with new values in \code{kalman_params}
#' @export
select_Kalman_variances <- function(ssm, X, y, method = 'igd', ...) {

  ### 1. Estimation of the hyper-parameters
  if (method == 'igd')
    l <- iterative_grid_search(X, y, ...)
  else if (method == 'em')
    l <- expectation_maximization(X, y, ...)
  else
    stop("Selection method is not recognized. It should be either 'igd' or 'em'.")

  ### 2. Creation of the new model
  ssm_new <- ssm
  ssm_new$kalman_params <- l
  ssm_new$kalman_params$opt_Kalman_call <- match.call()

  # ### 3. Compute the state-space inference for designed parameters
  # # Kalman Filtering and Smoothing
  # ssm_new <- predict(ssm_new, X, newy = y, type = 'model')

  ssm_new
}
