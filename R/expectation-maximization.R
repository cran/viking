# Copyright 2019-2022 EDF, Sorbonne Université and CNRS.
# Author : Joseph de Vilmarest (EDF, Sorbonne Université)
# The package Viking is distributed under the terms of the license GNU LGPL 3.

#' Expectation-Maximization
#'
#' \code{expectation_maximization} is a method to choose hyper-parameters of the
#' linear Gaussian State-Space Model with time-invariant variances relying on the
#' Expectation-Maximization algorithm.
#'
#' E-step is realized through recursive Kalman formulae (filtering then smoothing).\cr
#' M-step is the maximization of the expected complete likelihood with respect to the
#' hyper-parameters.\cr
#' We only have the guarantee of convergence to a LOCAL optimum.
#' We fix P1 = p1 I (by default p1 = 0). We optimize theta1, sig, Q.
#'
#' @param X explanatory variables
#' @param y time series
#' @param n_iter number of iterations of the EM algorithm
#' @param Q_init initial covariance matrix on the state noise
#' @param sig_init (optional, default \code{1}) initial value of the standard deviation
#' of the observation noise
#' @param verbose (optional, default \code{1000}) frequency for prints
#' @param lambda (optional, default \code{10^-9}) regularization parameter to avoid singularity
#' @param mode_diag (optional, default \code{FALSE}) if \code{TRUE} then we restrict the
#' search to diagonal matrices for \code{Q}
#' @param p1 (optional, default \code{0}) deterministic value of \code{P1 = p1 I}
#'
#' @return a list containing values for \code{P,theta,sig,Q}, and two vectors
#' \code{DIFF, LOGLIK} assessing the convergence of the algorithm.
#' @export
#'
#' @example tests/example_em.R
expectation_maximization <- function(X, y, n_iter, Q_init, sig_init = 1, verbose = 1000,
                                     lambda = 10^-9, mode_diag = FALSE, p1 = 0) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  sig2 <- sig_init^2
  theta1 <- matrix(0,d,1)
  Q <- Q_init
  DIFF <- numeric(n_iter)
  LOGLIK <- rep(- 0.5 * log(2 * pi), n_iter)
  flag_decrease <- FALSE
  init_time <- Sys.time()

  for (i in 1:n_iter) {
    theta_arr <- array(0,dim=c(n,d))
    theta_arr2 <- array(0,dim=c(n,d))
    theta <- matrix(0,d,1)
    P_arr <- array(0,dim=c(n,d,d))
    P <- diag(d,x=p1)
    C_arr <- array(0,dim=c(n,d,d))
    C_arr2 <- array(0,dim=c(n,d,d))
    C <- diag(d)

    ##################################################################################
    # E-step
    # filtering
    for (t in 1:n) {
      theta_arr2[t,] <- theta
      C_arr2[t,,] <- C
      Xt <- X[t,]
      LOGLIK[i] <- LOGLIK[i] - 0.5 * log(sig2 + crossprod(Xt, P %*% Xt)[1]) / n -
        0.5 * (y[t] - crossprod(theta + C %*% theta1, Xt)[1])^2 /
        (sig2 + crossprod(Xt, P %*% Xt)[1]) / n
      P <- P - tcrossprod(P%*%Xt) / (sig2 + crossprod(Xt, P%*%Xt)[1])
      P_arr[t,,] <- P
      theta <- theta + P %*% Xt * (y[t] - crossprod(theta, Xt)[1]) / sig2
      theta_arr[t,] <- theta
      C <- (diag(d) - tcrossprod(P%*%Xt, Xt) / sig2) %*% C
      C_arr[t,,] <- C
      P <- P + Q
    }
    Pn <- P - Q

    # smoothing
    Q_inv <- solve(Q)
    Q_new <- diag(d,x=0)
    sig2_new <- crossprod(X[n,], Pn %*% X[n,])[1] / n
    for (t in (n-1):1) {
      Pt <- P_arr[t,,]
      Pt_inv <- solve(Pt+Q)
      Q_new <- Q_new + (Pn - Pn %*% Pt_inv %*% Pt - t(Pn %*% Pt_inv %*% Pt)) / (n-1)
      thetat <- matrix(theta_arr[t,],d,1)
      theta_arr[t,] <- thetat + Pt %*% Pt_inv %*% (theta_arr[t+1,] - thetat)
      Pn <- Pt + Pt %*% Pt_inv %*% (Pn - Pt - Q) %*% Pt_inv %*% Pt
      Q_new <- Q_new + Pn / (n-1)
      sig2_new <- sig2_new + crossprod(X[t,], Pn %*% X[t,])[1] / n
      C <- C_arr[t,,] + Pt %*% Pt_inv %*% (C - C_arr[t,,])
      C_arr[t,,] <- C
    }

    ##################################################################################
    # M-step
    # update theta1
    A <- diag(d,x=lambda)
    b <- matrix(0,d,1)
    for (t in 1:n) {
      Xt <- X[t,]
      Ct <- C_arr2[t,,]
      A <- A + tcrossprod(crossprod(Ct, Xt))
      b <- b + (y[t] - crossprod(theta_arr2[t,], Xt)[1]) * crossprod(Ct, Xt)
    }
    theta1 <- solve(A,b)

    for (t in 1:n)
      theta_arr[t,] <- theta_arr[t,] + C_arr[t,,] %*% theta1
    # update Q
    Q_last <- Q
    Q <- Q_new
    for (t in 1:(n-1))
      Q <- Q + tcrossprod(theta_arr[t+1,] - theta_arr[t,]) / (n-1)
    Q <- (Q + t(Q)) / 2 # force symmetry to avoid approx error propagation
    if (mode_diag)
      Q <- diag(diag(Q))

    # update sig2
    sig2 <- sig2_new
    for (t in 1:n)
      sig2 <- sig2 + (y[t] - crossprod(theta_arr[t,], X[t,])[1])^2 / n

    DIFF[i] <- sqrt(sum(diag(crossprod(Q-Q_last))) / sum(diag(crossprod(Q_last))))
    if (verbose > 0) {
      if (i %% verbose == 0) {
        print(paste('Iteration:',i,
                    '| Relative variation of Q:', format(DIFF[i], scientific=TRUE, digits=3),
                    '| Log-likelihood:', round(LOGLIK[i], digits=6)))
        flag_decrease <- FALSE
      }
      if (sum(sum(abs(Q-t(Q)))>10^-6)>0)
        print('Q not symmetric')
      if (min(eigen(Q, symmetric=TRUE)$values) < 0)
        print('Q is not positive')
    }
    if (i > 1) {
      if (LOGLIK[i] < LOGLIK[i-1] & flag_decrease == FALSE) {
        warning('WARNING: log-likelihood decreased')
        flag_decrease <- TRUE
      }
    }
  }
  if (verbose > 0)
    print(c('Computation time of the expectation-maximization:', Sys.time()-init_time))
  list(P = diag(d,x=p1), theta = theta1, Q = Q, sig = sqrt(sig2), DIFF = DIFF, LOGLIK = LOGLIK)
}
