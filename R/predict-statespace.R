# Copyright 2019-2022 EDF, Sorbonne Université and CNRS.
# Author : Joseph de Vilmarest (EDF, Sorbonne Université)
# The package Viking is distributed under the terms of the license GNU LGPL 3.

#' Predict using a statespace object
#'
#' \code{predict.statespace} makes a prediction for a statespace object, in the offline or online
#' setting.
#'
#' @param object the statespace object
#' @param newX the design matrix in the prediction set
#' @param newy (default \code{NULL}) the variable of interest in the prediction set. If specified
#' it allows to use the state-space model in the online setting. Otherwise the prediction is
#' offline.
#' @param online (default \code{TRUE}) specifies if the prediction is made online, that is if
#' the observation at time t-1 is used to update the model before predicting at time t.
#' @param compute_smooth (default \code{FALSE}) specifies if Kalman Smoothing is also computed.
#' @param type type of prediction. Can be either
#' \describe{
#' \item{mean}{return the mean forecast.}
#' \item{proba}{return a probabilistic forecast (list containing estimation of the mean and
#' standard deviation).}
#' \item{model}{return the updated statespace object (containing also the forecasts).}
#' }
#' @param ... additional parameters
#'
#' @return Depending on the type specified, the result is \cr
#' - a vector of mean forecast if \code{type='mean'}
#' - a list of two vectors, mean forecast and standard deviations if \code{type='proba'}
#' - a statespace object if \code{type='model'}
#'
#' @export
predict.statespace <- function(object, newX, newy = NULL, online = TRUE, compute_smooth = FALSE,
                               type = c('mean', 'proba', 'model'), ...) {
  type <- match.arg(type)

  object_new <- object
  X <- if (is.null(ncol(newX))) matrix(newX, length(newX), 1) else newX
  object_new$X <- X
  object_new$y <- newy

  if (online) {
    if (is.null(newy))
      stop('newy must be given when online = TRUE')

    if (is.null(object$viking_params)) {
      if (compute_smooth) {
        object_new$ks <- kalman_smoothing(X, newy, object$kalman_params$theta, object$kalman_params$P,
                                       Q = object$kalman_params$Q, sig = object$kalman_params$sig)
      }

      object_new$kf <- kalman_filtering(X, newy, object$kalman_params$theta, object$kalman_params$P,
                                     Q = object$kalman_params$Q, sig = object$kalman_params$sig)

      object_new$kalman_params$theta <- object_new$kf$theta
      object_new$kalman_params$P <- object_new$kf$P
      object_new$pred_mean <- sapply(1:nrow(X), function(t) {crossprod(object_new$kf$theta_arr[t,], X[t,])[1]})
      object_new$pred_sd <- sapply(1:nrow(X), function(t) {
        sqrt(crossprod(X[t,], object_new$kf$P_arr[t,,] %*% X[t,]) + object_new$kalman_params$sig^2)
      })
    } else {
      object_new$vik <- do.call(viking, c(list(X = newX, y = newy), object$viking_params))
      object_new$pred_mean <- sapply(1:nrow(X), function(t) {crossprod(object_new$vik$theta_arr[t,], X[t,])[1]})
      object_new$pred_sd <- sapply(1:nrow(X), function(t) {
        sqrt(crossprod(X[t,], object_new$vik$P_arr[t,,] %*% X[t,]))
      }) + exp(object_new$vik$hata + object_new$vik$s / 2)
    }
  } else {
    theta <- if (is.null(object$viking_params)) object$kalman_params$theta else object$viking_params$theta
    object_new$pred_mean <- sapply(1:nrow(X), function(t) {
      crossprod(X[t,], theta)
    })
    P <- if (is.null(object$viking_params)) object$kalman_params$P else object$viking_params$P
    object_new$pred_sd <- sapply(1:nrow(X), function(t) {
      sqrt(crossprod(X[t,], (P + (t-1) * object$kalman_params$Q) %*% X[t,]) + object$kalman_params$sig^2)
    })
  }

  switch(type, mean = object_new$pred_mean,
         proba = list(pred_mean = object_new$pred_mean, pred_sd = object_new$pred_sd),
         model = object_new)
}
