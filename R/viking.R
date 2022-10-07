# Copyright 2019-2022 EDF, Sorbonne Université and CNRS.
# Author : Joseph de Vilmarest (EDF, Sorbonne Université)
# The package Viking is distributed under the terms of the license GNU LGPL 3.

#' Viking: Variational bayesIan variance tracKING
#'
#' \code{viking} is the state-space estimation realized by Viking,
#' generalizing the Kalman Filter to variance uncertainty.
#'
#' @param X the explanatory variables
#' @param y the time series
#' @param theta0 initial \code{theta}
#' @param P0 initial \code{P}
#' @param hata0 initial \code{hata}
#' @param s0 initial \code{s}
#' @param hatb0 initial \code{hatb}
#' @param Sigma0 initial \code{Sigma}
#' @param n_iter (optional, default \code{2}) number of alternate steps
#' @param mc (optional, default \code{10}) number of Monte-Carlo samples
#' @param rho_a (optional, default \code{0}) learning rate of \code{a}
#' @param rho_b (optional, default \code{0}) learning rate of \code{b}
#' @param learn_sigma (optional, default \code{TRUE}) asserts the estimation of \code{a}
#' @param learn_Q (optional, default \code{TRUE}) asserts the estimation of \code{b}
#' @param K (optional, default \code{NULL}) if not \code{NULL} then it is a multiplicative
#' factor of the state in the state update
#' @param mode (optional, default \code{'diagonal'})
#' @param thresh (optional, default \code{TRUE})
#' @param phi (optional, default \code{logt})
#' @param phi1 (optional, default \code{logt1})
#' @param phi2 (optional, default \code{logt2})
#'
#' @importFrom stats rnorm
#'
#' @return a list composed of the evolving value of all the parameters:
#' \code{theta_arr, P_arr, q_arr, hata_arr, s_arr, hatb_arr, Sigma_arr}
#'
#' @references J. de Vilmarest, O. Wintenberger (2021), Viking: Variational Bayesian Variance
#' Tracking. <arXiv:2104.10777>
#'
#' @export
viking <- function(X, y, theta0, P0, hata0, s0, hatb0, Sigma0, n_iter = 2, mc = 10,
                   rho_a = 0, rho_b = 0, learn_sigma = TRUE, learn_Q = TRUE, K = NULL,
                   mode = 'diagonal', thresh = TRUE,
                   phi = logt, phi1 = logt1, phi2 = logt2) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  if (is.null(K))
    K <- diag(d)

  q_arr <- matrix(0,n,d)
  theta_arr <- matrix(0,n,d)
  P_arr <- array(0,dim=c(n,d,d))
  hata_arr <- numeric(n)
  s_arr <- numeric(n)
  hatb_arr <- matrix(0,n,d)
  Sigma_arr <- array(0,dim=c(n,d,d))

  theta <- theta0
  P <- P0
  hata <- hata0
  s <- s0
  hatb <- hatb0
  Sigma <- Sigma0

  for (t in 1:n) {

    hata_new <- hata
    s_new <- s + rho_a
    hatb_new <- hatb
    Sigma_new <- Sigma + rho_b * (if (mode == 'scalar') 1 else diag(d))

    if (mode == 'scalar')
      q_arr[t,] <- mean(phi( rnorm(10^3, mean=hatb_new, sd=Sigma_new) ))
    else if (mode == 'diagonal')
      q_arr[t,] <- colMeans(phi( mvtnorm::rmvnorm(10^3, mean=hatb_new, sigma=Sigma_new) ))
    theta_arr[t,] <- K %*% theta
    P_arr[t,,] <- K %*% P %*% t(K) + diag(diag(q_arr[t,]))

    hata_arr[t] <- hata_new
    s_arr[t] <- s_new
    hatb_arr[t,] <- hatb_new
    Sigma_arr[t,,] <- Sigma_new

    for (iter in 1:n_iter) {
      if (mode == 'scalar') {
        if (mc > 0) {
          bMC <- rnorm(mc, mean=hatb_new, sd=Sigma_new)
          A <- diag(d,x=0)
          for (b in bMC)
            A <- A + solve(K %*% P %*% t(K) + diag(d) * phi(b)) / length(bMC)
        } else {
          C_inv <- solve(K %*% P %*% t(K) + phi(hatb_new) * diag(d))
          A <- C_inv - phi2(hatb_new) * Sigma_new * C_inv %*% C_inv / 2 +
            phi1(hatb_new)^2 * Sigma_new * C_inv %*% C_inv %*% C_inv
        }
      }
      else if (mode == 'diagonal') {
        if (mc > 0) {
          bMC <- mvtnorm::rmvnorm(mc, mean=hatb_new, sigma=Sigma_new)
          A <- diag(d,x=0)
          for (i in 1:nrow(bMC))
            A <- A + solve(K %*% P %*% t(K) + diag(c( phi(bMC[i,]) ))) / nrow(bMC)
        } else {
          C_inv <- solve(K %*% P %*% t(K) + diag(c( phi(hatb_new) )))
          A <- C_inv - C_inv %*% diag(c( phi2(hatb_new) * diag(Sigma_new) )) %*% C_inv / 2 +
            C_inv %*% (C_inv * tcrossprod(phi1(hatb_new)) * Sigma_new) %*% C_inv
        }
      }
      A_inv <- solve((A+t(A))/2)

      # update of theta and P
      v <- exp(hata_new - s_new/2)
      P_new <- A_inv - tcrossprod(A_inv %*% X[t,]) / (crossprod(X[t,], A_inv %*% X[t,])[1] + v)
      theta_new <- K %*% theta + P_new %*% X[t,] / v * (y[t] - crossprod(K %*% theta, X[t,])[1])

      # update of hata and s
      if (learn_sigma) {
        c <- (y[t] - crossprod(theta_new, X[t,])[1])^2 + crossprod(X[t,], P_new %*% X[t,])[1]
        s_new <- ((s + rho_a)^-1 + 0.5 * c * exp(-hata_new))^-1
        Ma <- 3 * s
        diff <- 0.5 * ((s + rho_a)^-1 + 0.5 * c * exp(- hata + s_new / 2 + Ma))^-1 *
          (c * exp(- hata + s_new / 2) - 1)
        hata_new <- hata + max(min(diff, Ma), -Ma)
      }

      # update of hatb and Sigma
      if (learn_Q) {
        B <- P_new + tcrossprod(theta_new - K %*% theta)

        if (mode == 'scalar') {
          C_inv <- solve(K %*% P %*% t(K) + phi(hatb) * diag(d))
          g <- sum(diag( C_inv %*% (diag(d) - B %*% C_inv) )) * phi1(hatb)
          H <- - sum(diag(C_inv %*% B %*% C_inv)) * phi2(hatb) +
            2 * sum(diag(C_inv %*% C_inv %*% B %*% C_inv)) * phi1(hatb)^2
          Sigma_new <- ((Sigma + rho_b)^-1 + H / 2)^-1
          hatb_new <- hatb - Sigma_new * g / 2
        } else if (mode == 'diagonal') {
          C_inv <- solve(K %*% P %*% t(K) + diag(c( phi(hatb) )))
          g <- diag( C_inv %*% (diag(d) - B %*% C_inv) ) * phi1(hatb)
          H <- - (C_inv %*% B %*% C_inv %*% diag(c( phi2(hatb) ))) * diag(d) +
            2 * (C_inv %*% B %*% C_inv) * C_inv * tcrossprod(phi1(hatb))
          Sigma_new <- solve(solve(Sigma + rho_b * diag(d)) + H / 2)
          hatb_new <- hatb - Sigma_new %*% g / 2
        }
        if (thresh)
          hatb_new <- (hatb_new + abs(hatb_new)) / 2
      }
    }

    P <- P_new
    theta <- theta_new
    hata <- hata_new
    s <- s_new
    hatb <- hatb_new
    Sigma <- Sigma_new
  }
  list(theta_arr=theta_arr, P_arr=P_arr, q_arr=q_arr, hata_arr=hata_arr, s_arr=s_arr,
       hatb_arr=hatb_arr, Sigma_arr=Sigma_arr)
}

logt <- function(b) {log(1 + (b+abs(b)) / 2)}
logt1 <- function(b) {1 / (1 + (b+abs(b)) / 2) * (b >= 0)}
logt2 <- function(b) {-1 / (1 + (b+abs(b)) / 2)^2 * (b >= 0)}

