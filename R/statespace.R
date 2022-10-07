# Copyright 2019-2022 EDF, Sorbonne Université and CNRS.
# Author : Joseph de Vilmarest (EDF, Sorbonne Université)
# The package Viking is distributed under the terms of the license GNU LGPL 3.

#' Design a State-Space Model
#'
#' The function \code{statespace} builds a state-space model, with known or unknown variances.
#' By default, this function builds a state-space model in the static setting, with a constant
#' state (zero state noise covariance matrix) and constant observation noise variance.
#'
#' @param X design matrix.
#' @param y variable of interest.
#' @param kalman_params (default \code{NULL}) list containing initial values for \code{theta,P}
#' as well as the variances (\code{Q,sig}). If it is not specified, the state-space model is
#' constructed in the static setting (\code{theta=0, P=I, Q=0, sig=1}).
#' @param viking_params (default \code{NULL}) list of parameters for the Viking algorithm.
#' @param ... additional parameters
#'
#' @return a statespace object.
#'
#' @importFrom stats predict
#' @export
#'
#' @example tests/example_statespace.R
statespace <- function(X, y, kalman_params = NULL, viking_params = NULL, ...) {

  ### Create state-space object
  d <- if (is.null(ncol(X))) 1 else ncol(X)

  kp <- list(theta = filter_null(kalman_params$theta, matrix(0,d,1)),
             P = filter_null(kalman_params$P, diag(d)),
             Q = filter_null(kalman_params$Q, diag(d, x=0)),
             sig = filter_null(kalman_params$sig, 1),
             opt_Kalman_call = kalman_params$opt_Kalman_call)

  vp <- viking_params
  if (!is.null(viking_params)) {
    vp$theta <- filter_null(vp$theta, matrix(0,d,1))
    vp$P <- filter_null(vp$P, diag(d))
    vp$hata <- filter_null(vp$hata, 0)
    vp$s <- filter_null(vp$s, 0)
    viking_params$mode <- filter_null(viking_params$mode, 'diagonal')
    vp$hatb <- filter_null(vp$hatb, if (viking_params$mode == 'scalar') 0 else matrix(0, d, 1))
    vp$Sigma <- filter_null(vp$Sigma, if (viking_params$mode == 'scalar') 0 else diag(d, x=0))
  }

  ssm <- list(kalman_params = kp, viking_params = vp)
  class(ssm) <- 'statespace'

  ### Compute the state-space inference for designed parameters
  # Using either Kalman Filtering and Smoothing or Viking
  ssm <- predict(ssm, X, newy = y, type = 'model', ...)

  ssm
}

filter_null <- function(x, y) {if (is.null(x)) y else x}
