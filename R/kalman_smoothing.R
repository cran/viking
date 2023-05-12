# Copyright 2019-2022 EDF, Sorbonne Université and CNRS.
# Author : Joseph de Vilmarest (EDF, Sorbonne Université)
# The package Viking is distributed under the terms of the license GNU LGPL 3.

#' @title Kalman Smoothing
#' @description Compute the smoothed estimation of the parameters \code{theta} and \code{P}.
#'
#' @param X the explanatory variables
#' @param y the time series
#' @param theta1 initial \code{theta}
#' @param P1 initial \code{P}
#' @param Q (optional, default \code{0}) covariance matrix of the state noise
#' @param sig (optional, default \code{1}) variance of the spate noise
#'
#' @return a list containing \code{theta_arr} and \code{P_arr}, the smoothed estimation of
#' the parameters \code{theta} and \code{P}.
#' @export
kalman_smoothing <- function(X, y, theta1, P1, Q=0, sig=1) {
  n <- nrow(X)
  d <- ncol(X)

  theta_arr <- array(0,dim=c(n,d))
  theta <- theta1
  P_arr <- array(0,dim=c(n,d,d))
  P <- P1

  # filtering
  for (t in 1:n) {
    Xt <- X[t,]
    if (sum(is.na(c(X[t,],y[t]))) == 0)
      P <- P - tcrossprod(P %*% Xt) / (sig^2 + crossprod(Xt, P %*% Xt)[1])
    theta_arr[t,] <- theta
    P_arr[t,,] <- P
    if (sum(is.na(c(X[t,],y[t]))) == 0) {
      theta <- theta + P %*% Xt * (y[t] - crossprod(theta, Xt)[1]) / sig^2
      P <- P + Q
    }
  }
  Pn <- P - Q

  # smoothing
  for (t in (n-1):1) {
    Pt <- P_arr[t,,]
    P_inv <- solve(Pt+Q)
    thetat <- matrix(theta_arr[t,],d,1)
    theta_arr[t,] <- thetat + Pt %*% P_inv %*% (theta_arr[t+1,] - thetat)
    Pn <- Pt + Pt %*% P_inv %*% (Pn - Pt - Q) %*% P_inv %*% Pt
    P_arr[t,,] <- Pn
  }

  list(theta_arr = theta_arr, P_arr = P_arr)
}
