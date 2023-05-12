# Copyright 2019-2022 EDF, Sorbonne Université and CNRS.
# Author : Joseph de Vilmarest (EDF, Sorbonne Université)
# The package Viking is distributed under the terms of the license GNU LGPL 3.

#' Plot a statespace object
#'
#' \code{plot.statespace} displays different graphs expressing the behavior of the state-space
#' model:\cr
#' 1. Evolution of the Bias: rolling version of the error of the model.\cr
#' 2. Evolution of the RMSE: root-mean-square-error computed on a rolling window.\cr
#' 3. State Evolution: time-varying state coefficients, subtracted of the initial state vector.\cr
#' 4. Normal Q-Q Plot: we check if the observation follows the Gaussian distribution of estimated
#' mean and variance. To that end, we display a Q-Q plot of the residual divided by the estimated
#' standard deviation, against the standard normal distribution.\cr
#'
#' @param x the statespace object.
#' @param pause (default \code{FALSE}) if set to \code{FALSE} then the plots are displayed on a single
#' page, otherwise a new page is created for each plot.
#' @param window_size (default \code{7}) the window size of the rolling mean computed on the
#' error to display the bias, and on the mean squared error to display a rolling RMSE.
#' @param date (default \code{NULL}) defines the values for the x-axis.
#' @param sel (default \code{NULL}) defines a subset of the data on which we zoom.
#' For instance one can display the evolution of the SSM on a test set and not the whole data set.
#' @param ... additional parameters
#'
#' @importFrom graphics axis box legend lines mtext par plot points text
#' @importFrom stats qqnorm sd
#'
#' @return No return value, called to display plots.
#'
#' @export
plot.statespace <- function(x, pause = FALSE, window_size = 7, date = NULL, sel = NULL, ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (pause) par(ask = TRUE, mfrow = c(1,1))
  else par(ask = FALSE, mfrow = c(3,2), mai=c(0.4,0.4,0.4,0.4))
  if (is.null(date)) date <- 1:length(x$y)
  if (is.null(sel)) sel <- 1:length(x$y)

  # 1. Evolution of the bias
  res <- (x$y - x$pred_mean)[sel]
  th_fixed <- if (is.null(x$vik)) x$kf$theta_arr[sel[1],] else x$vik$theta_arr[sel[1],]
  res_fixed <- x$y[sel] - x$X[sel,] %*% th_fixed
  date_sel <- date[sel]
  lfixed <- sapply(window_size:length(res), function(t) {mean(res_fixed[(t-window_size+1):t], na.rm=T)})
  l <- sapply(window_size:length(res), function(t) {mean(res[(t-window_size+1):t], na.rm=T)})
  plot(date_sel[window_size:length(res)], lfixed, type='l', ylim=range(c(lfixed, l), na.rm=T),
       main='Evolution of the Bias', ylab='', xlab='', lwd=2)
  lines(date_sel[window_size:length(res)], l, col='darkgreen', lwd=2)
  mtext('SSM', side = 4, at = l[length(l)], col='darkgreen', las=2, cex=0.5*(1+pause), line=0.3)
  mtext('Fixed', side = 4, at = lfixed[length(lfixed)], col='black', las=2, cex=0.5*(1+pause), line=0.3)

  # 2. Evolution of the RMSE
  lfixed <- sapply(window_size:length(res), function(t) {sqrt(mean(res_fixed[(t-window_size+1):t]^2, na.rm=T))})
  l <- sapply(window_size:length(res), function(t) {sqrt(mean(res[(t-window_size+1):t]^2, na.rm=T))})
  plot(date_sel[window_size:length(res)], lfixed, ylim=c(0,max(c(lfixed, l), na.rm=T)), type='l',
       ylab='', xlab='', lwd=2,
       main=if (pause) 'Evolution of the Root-Mean-Square-Error' else 'Evolution of the RMSE')
  lines(date_sel[window_size:length(res)], l, col='darkgreen', lwd=2)
  mtext('SSM', side=4, at=l[length(l)], col='darkgreen', las=2, cex=0.5*(1+pause), line=0.3)
  mtext('Fixed', side=4, at=lfixed[length(lfixed)], col='black', las=2, cex=0.5*(1+pause), line=0.3)

  # 3. Evolution of the state coefficients
  # Kalman Filtering
  plot_evol <- function(evol, title, coln = NULL, logpar='') {
    c <- RColorBrewer::brewer.pal(max(ncol(evol), 11), name = "Spectral")[sapply(1:ncol(evol), function(i) min(i,11))]
    plot(date_sel, evol[,1], type='l', ylim=range(evol), col=c[1],
         main=title, ylab='', xlab='', lwd=2, log=logpar)
    for (j in 2:ncol(evol))
      lines(date_sel, evol[,j], col=c[j], lwd=2)
    if (!is.null(coln))
      mtext(side=4, text=coln, at=evol[nrow(evol),], las=2, col=c, cex=0.5, line=0.3)
  }
  if (!is.null(x$kf)) {
    evoltheta <- x$kf$theta_arr[sel,]
    evoltheta <- evoltheta - matrix(rep(evoltheta[1,], nrow(evoltheta)), nrow(evoltheta),
                                    ncol(evoltheta), byrow=TRUE)
    plot_evol(evoltheta, 'State Evolution: Kalman Filtering', coln = colnames(x$X))
  }
  # Kalman Smoothing
  if (!is.null(x$ks)) {
    evoltheta <- x$ks$theta_arr[sel,]
    evoltheta <- evoltheta - matrix(rep(evoltheta[1,], nrow(evoltheta)), nrow(evoltheta),
                                    ncol(evoltheta), byrow=TRUE)
    plot_evol(evoltheta, 'State Evolution: Kalman Smoothing', coln = colnames(x$X))
  }
  # Viking
  if (!is.null(x$vik)) {
    evoltheta <- x$vik$theta_arr[sel,]
    evoltheta <- evoltheta - matrix(rep(evoltheta[1,], nrow(evoltheta)), nrow(evoltheta),
                                    ncol(evoltheta), byrow=TRUE)
    plot_evol(evoltheta, 'State Evolution: Viking', coln = colnames(x$X))
  }

  # 4. Viking: evolution of the variances
  if (!is.null(x$vik)) {
    plot(date_sel, exp(x$vik$hata_arr[sel] + x$vik$s_arr[sel] / 2), type='l', lwd=2,
         main='Evolution of the observation noise variance', ylab='', xlab='')
    plot_evol(x$vik$q_arr[sel,], 'Evolution of the state noise covariance matrix',
              coln = colnames(x$X), logpar = 'y')
  }

  # 5. Plot the diagonal of Q
  if (!is.null(x$kf)) {
    diagQ <- diag(x$kalman_params$Q)
    # diagQ[diagQ == 0] <- NA
    c <- RColorBrewer::brewer.pal(max(length(diagQ), 11), name = "Spectral")[sapply(1:length(diagQ), function(i) min(i,11))]
    plot(c(1:length(diagQ)), diagQ, type='h', xlab='', ylab='', col=c, axes=FALSE,
         main='Diagonal of Q', ylim=c(0,max(diagQ)*1.4))
    text(c(1:length(diagQ)), diagQ, labels=format(diagQ, scientific=TRUE, digits=2),
         srt=90, adj=c(-0.15,0.5), col=c)
    points(c(1:length(diagQ)), diagQ, pch=20, col=c)
    axis(2)
    box()
    mtext(side=1, text=colnames(x$X), at=1:length(diagQ), las=2, col=c, cex=0.5)
  }

  # 6. Check of the Gaussian assumption
  qqnorm(res / x$pred_sd[sel], xlab='', ylab='', main='Normal Q-Q Plot of the Residuals')
  lines(-4:4, -4:4, col='darkred', lwd=2)
  lines(-4:4, -4:4 * sd(res / x$pred_sd[sel]), col='darkgreen', lwd=2)
  legend('topleft', legend=c('Standard Normal', 'Best Constant'), col=c('darkred','darkgreen'),
         lty=1, lwd=2)

  par(mfrow = c(1,1))
}
