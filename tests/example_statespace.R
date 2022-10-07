set.seed(1)
### Simulate data
n <- 1000
d <- 5
Q <- diag(c(0,0,0.25,0.25,0.25))
sig <- 1

X <- cbind(matrix(rnorm((d-1)*n,sd=1),n,d-1),1)
theta <- matrix(rnorm(d), d, 1)
theta_arr <- matrix(0, n, d)
for (t in 1:n) {
  theta_arr[t,] <- theta
  theta <- theta + matrix(mvtnorm::rmvnorm(1,matrix(0,d,1),Q),d,1)
}
y <- rowSums(X * theta_arr) + rnorm(n, sd=sig)

####################
### Kalman Filter
# Default Static Setting
ssm <- viking::statespace(X, y)
plot(ssm)

# Known variances
ssm <- viking::statespace(X, y, kalman_params = list(Q=Q, sig=sig))
plot(ssm)

