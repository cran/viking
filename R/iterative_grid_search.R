# Copyright 2019-2022 EDF, Sorbonne Université and CNRS.
# Author : Joseph de Vilmarest (EDF, Sorbonne Université)
# The package Viking is distributed under the terms of the license GNU LGPL 3.

#' Iterative Grid Search
#'
#' \code{iterative_grid_search} is an iterative method to choose hyper-parameters of
#' the linear Gaussian State-Space Model with time-invariant variances.
#'
#' We restrict ourselves to a diagonal matrix \code{Q} and we optimize \code{Q / sig^2} on
#' a grid. Each diagonal coefficient is assumed to belong to a pre-defined \code{q_list}.\cr
#' We maximize the log-likelihood on that region of search in an iterative fashion.
#' At each step, we change the diagonal coefficient improving the most the log-likelihood.
#' We stop when there is no possible improvement. This doesn't guarantee an optimal point
#' on the grid, but the computational time is much lower.
#'
#' @param X the explanatory variables
#' @param y the observations
#' @param q_list the possible values of \code{diag(Q) / sig^2}
#' @param Q_init (default \code{NULL}) initial value of \code{Q / sig^2},
#' if \code{NULL} it is set to 0
#' @param max_iter (optional 0) maximal number of iterations. If 0 then optimization is
#' done as long as we can improve the log-likelihood
#' @param delay (optional, default 1) to predict \code{y[t]} we have access to \code{y[t-delay]}
#' @param use (optional, default \code{NULL}) the availability variable
#' @param restrict (optional, default \code{NULL}) if not \code{NULL} it allows to specify the
#' indices of the diagonal coefficient to optimize
#' @param mode (optional, default \code{gaussian})
#' @param p1 (optional, default \code{0}) coefficient for \code{P1/sig^2 = p1 I}
#' @param ncores (optional, default \code{1}) number of available cores for computation
#' @param train_theta1 (optional, default \code{NULL}) training set for estimation of \code{theta1}
#' @param train_Q (optional, default \code{NULL}) time steps on which the log-likelihood is computed
#' @param verbose (optional, default \code{TRUE}) whether to print intermediate progress
#'
#' @return a list containing values for \code{P,theta,sig,Q}, as well as \code{LOGLIK},
#' the evolution of the log-likelihood during the search.
#' @export
#'
#' @example tests/example_igd.R
iterative_grid_search <- function(X, y, q_list, Q_init = NULL, max_iter = 0, delay = 1,
                                  use = NULL, restrict = NULL, mode = 'gaussian', p1 = 0,
                                  ncores = 1, train_theta1 = NULL, train_Q = NULL, verbose=TRUE) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  train_theta1 <- filter_null(train_theta1, 1:n)
  train_Q <- filter_null(train_Q, 1:n)
  Qstar <- if(is.null(Q_init)) diag(d,x=0) else diag(diag(Q_init))
  use <- filter_null(use, sapply(1:n - delay, function(x) max(x,0)))
  search_dimensions <- filter_null(restrict, 1:d)
  b <- TRUE
  l_opt <- loglik(X, y, Qstar, use, p1, train_theta1, train_Q, mode = mode)
  n_iter <- 0
  LOGLIK <- c()
  init_time <- Sys.time()

  while (b) {
    n_iter <- n_iter + 1
    b <- n_iter < max_iter | max_iter == 0
    i <- -1
    q_i <- 0
    for (k in search_dimensions) {
      q_prev <- Qstar[k,k]
      l_arr <- unlist(parallel::mclapply(q_list,function(q) {
        Qstar[k,k] <- q
        loglik(X, y, Qstar, use, p1, train_theta1, train_Q, mode = mode)
      }, mc.cores=ncores))
      if (max(l_arr) > l_opt) {
        l_opt <- max(l_arr)
        i <- k
        q_i <- q_list[which.max(l_arr)]
      }
      Qstar[k,k] <- q_prev
    }
    LOGLIK <- c(LOGLIK, l_opt)
    if (i == -1)
      b <- FALSE
    else {
      Qstar[i,i] <- q_i
      if (verbose) {
        if (q_i %in% c(min(q_list),max(q_list)))
          warning(c('Values may not be suited',q_i))
        print(paste('Iteration:', n_iter,
                    ' | Log-likelihood:', round(l_opt, digits=6),
                    ' | Diagonal of Q/sigma^2:'))
        print(diag(Qstar))
      }
    }
  }
  par <- parameters_star(X, y, Qstar, p1)
  sig <- get_sig(X, y, par, Qstar, use)
  if (verbose)
    print(paste('Computation time of the iterative grid search:',
                format(difftime(Sys.time(),init_time))))

  list(theta = get_theta1(X, y, par, Qstar, use, mode = mode), P = diag(d,x=p1*sig^2),
       sig = sig, Q = Qstar*sig^2, LOGLIK = LOGLIK)
}
