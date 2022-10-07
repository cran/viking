# Copyright 2019-2022 EDF, Sorbonne Université and CNRS.
# Author : Joseph de Vilmarest (EDF, Sorbonne Université)
# The package Viking is distributed under the terms of the license GNU LGPL 3.

parameters_star <- function(X, y, Qstar, p1 = 0) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  thetastar <- matrix(0,d,1)
  Pstar <- diag(d,x=p1)
  C <- diag(d,x=1)
  thetastar_arr <- array(0,dim=c(n+1,d))
  Pstar_arr <- array(0,dim=c(n+1,d,d))
  C_arr <- array(0,dim=c(n+1,d,d))
  C_arr[1,,] <- C

  for (t in 1:n) {
    Xt <- X[t,]
    err <- y[t] - crossprod(thetastar, Xt)[1]
    inv <- 1 / (1 + crossprod(Xt, Pstar %*% Xt)[1])

    thetastar_new <- thetastar + Pstar %*% Xt * inv * err
    Pstar_new <- Pstar + Qstar - tcrossprod(Pstar %*% Xt) * inv
    C_new <- (diag(d) - tcrossprod(Pstar %*% Xt, Xt) * inv) %*% C

    thetastar <- thetastar_new
    Pstar <- Pstar_new
    C <- C_new

    thetastar_arr[t+1,] <- thetastar
    Pstar_arr[t+1,,] <- Pstar
    C_arr[t+1,,] <- C
  }
  list(thetastar_arr=thetastar_arr, Pstar_arr=Pstar_arr, C_arr=C_arr)
}

get_theta1 <- function(X, y, par, Qstar, use, mode = 'gaussian') {
  n <- dim(X)[1]
  d <- dim(X)[2]
  ##########################################################
  A <- diag(d,x=0)
  b <- matrix(0,d,1)

  for (t in 1:n) {
    Xt <- X[t,]
    err <- y[t] - crossprod(par$thetastar_arr[use[t]+1,], Xt)[1]

    inv <- 1
    if (mode == 'gaussian') {
      inv <- 1 / (1 + crossprod(Xt, (par$Pstar_arr[use[t]+1,,] +
                                       max(0,t-use[t]-1) * Qstar) %*% Xt)[1])
    }

    A <- A + tcrossprod(crossprod(par$C_arr[use[t]+1,,],Xt)) * inv
    b <- b + crossprod(par$C_arr[use[t]+1,,],Xt) * err * inv
  }
  solve(A,b)
}

get_sig <- function(X, y, par, Qstar, use) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  theta1 <- get_theta1(X, y, par, Qstar, use)
  sig2 <- 0
  for (t in 1:n) {
    Xt <- matrix(X[t,],d,1)
    P <- par$Pstar_arr[use[t]+1,,] + max(t-use[t]-1,0) * Qstar
    sig2 <- sig2 + (y[t] - crossprod(par$thetastar_arr[use[t]+1,] +
                                       par$C_arr[use[t]+1,,]%*%theta1,Xt)[1])^2 /
      (1 + crossprod(Xt, P%*%Xt)[1]) / n
  }
  sqrt(sig2)
}

#' Log-likelihood
#'
#' \code{loglik} computes the log-likelihood of a state-space model of specified
#' \code{Q/sig^2, P1/sig^2, theta1}.
#'
#' @param X explanatory variables
#' @param y time series
#' @param Qstar the ratio \code{Q/sig^2}
#' @param use the availability variable
#' @param p1 coefficient for \code{P1/sig^2 = p1 I}
#' @param train_theta1 training set for estimation of \code{theta1}
#' @param train_Q time steps on which the log-likelihood is computed
#' @param mode (optional, default \code{gaussian})
#'
#' @return a numeric value for the log-likelihood.
#' @export
loglik <- function(X, y, Qstar, use, p1, train_theta1, train_Q, mode = 'gaussian') {
  par <- parameters_star(X, y, Qstar, p1)
  theta1 <- get_theta1(X[train_theta1,], y[train_theta1], par, Qstar, use,
                       mode = mode)
  sum1 <- 0
  sum2 <- 0
  n <- length(train_Q)
  for (t in train_Q) {
    Xt <- X[t,]
    P <- par$Pstar_arr[use[t]+1,,] + max(0,t-use[t]-1) * Qstar
    sum1 <- sum1 + log(1 + crossprod(Xt, P%*%Xt)[1]) / n
    err2 <- (y[t] - crossprod(par$thetastar_arr[use[t]+1,] +
                                par$C_arr[use[t]+1,,]%*%theta1,Xt)[1])^2
    if (mode == 'gaussian')
      sum2 <- sum2 + err2 / (1 + crossprod(Xt, P%*%Xt)[1]) / n
    else if (mode == 'rmse')
      sum2 <- sum2 + err2 / n
  }
  if (mode == 'gaussian')
    return(- 0.5 * sum1 - 0.5 - 0.5 * log(2 * pi * sum2))
  else if (mode == 'rmse')
    return(- sqrt(sum2))
}
