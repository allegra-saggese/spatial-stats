###############################################################################
# HW2 Problems 4-6: K-L representation of exponential correlation
# Range parameter phi = 1, domain [-L, L] with L = 5
###############################################################################

library(geoR)
library(fields)

set.seed(42)

phi <- 1
L   <- 5
s_grid <- seq(-L, L, length.out = 200)

###############################################################################
# P4: True K-L approximation
# Solve transcendental equations tan(wL) = phi/w  and  tan(wL) = -w/phi
###############################################################################

# Find roots of f(w) = tan(wL) - phi/w  in each branch (cosine family)
# The k-th branch of tan lies in ((k-1/2)pi/L, (k+1/2)pi/L); root sits inside.
find_w_cos <- function(k, L, phi) {
  lo <- (k - 0.5) * pi / L + 1e-6
  hi <- (k + 0.5) * pi / L - 1e-6
  if (k == 0) lo <- 1e-6  # first branch starts at 0
  uniroot(function(w) tan(w * L) - phi / w, c(lo, hi))$root
}

# Roots of tan(wL) = -w/phi  (sine family)
find_w_sin <- function(k, L, phi) {
  lo <- (k - 0.5) * pi / L + 1e-6
  hi <- (k + 0.5) * pi / L - 1e-6
  uniroot(function(w) tan(w * L) + w / phi, c(lo, hi))$root
}

# Build first N pairs (so 2N modes total)
build_KL <- function(N, L, phi) {
  w1 <- sapply(0:(N-1), find_w_cos, L = L, phi = phi)
  w2 <- sapply(1:N,     find_w_sin, L = L, phi = phi)
  
  lam1 <- 2 * phi / (w1^2 + phi^2)
  lam2 <- 2 * phi / (w2^2 + phi^2)
  
  list(w1 = w1, w2 = w2, lam1 = lam1, lam2 = lam2)
}

# Reconstruct rho(tau) = C(s, s+tau) with s = 0  (stationary => depends on tau)
# Use rho(tau) = sum_j lambda_j psi_j(0) psi_j(tau)
KL_corr <- function(tau, KL, L) {
  w1 <- KL$w1; w2 <- KL$w2
  lam1 <- KL$lam1; lam2 <- KL$lam2
  
  # Normalizers
  norm1 <- sqrt(L + sin(2 * w1 * L) / (2 * w1))
  norm2 <- sqrt(L - sin(2 * w2 * L) / (2 * w2))
  
  # cosine modes: psi(0) = 1/norm1; psi(tau) = cos(w1 tau)/norm1
  contrib1 <- sapply(seq_along(w1), function(j) {
    lam1[j] * (1 / norm1[j]) * (cos(w1[j] * tau) / norm1[j])
  })
  # sine modes: psi(0) = 0  -> contributes nothing to rho(tau) when s=0
  contrib2 <- sapply(seq_along(w2), function(j) {
    lam2[j] * (sin(w2[j] * 0) / norm2[j]) * (sin(w2[j] * tau) / norm2[j])
  })
  
  if (length(tau) == 1) sum(contrib1) + sum(contrib2)
  else rowSums(cbind(contrib1, contrib2))
}

# Above used s=0; but rho(tau) shouldn't depend on s
KL_cov_matrix <- function(s_pts, KL, L) { # set cov matrix
  w1 <- KL$w1; w2 <- KL$w2
  lam1 <- KL$lam1; lam2 <- KL$lam2
  norm1 <- sqrt(L + sin(2 * w1 * L) / (2 * w1))
  norm2 <- sqrt(L - sin(2 * w2 * L) / (2 * w2))
  
 
  m <- length(s_pts)
  C <- matrix(0, m, m)
  for (j in seq_along(w1)) {
    psi <- cos(w1[j] * s_pts) / norm1[j]
    C <- C + lam1[j] * outer(psi, psi)
  }
  for (j in seq_along(w2)) {
    psi <- sin(w2[j] * s_pts) / norm2[j]
    C <- C + lam2[j] * outer(psi, psi)
  }
  C
}

true_corr <- function(tau, phi = 1) exp(-phi * abs(tau)) # get corr

# Compare across a set of orders --- random 
orders <- c(5, 10, 25, 50)
tau_grid <- seq(0, 2*L, length.out = 300)

par(mfrow = c(2, 2))
for (N in orders) {
  KL <- build_KL(N, L, phi)
  # Use C(0, tau) from the covariance matrix evaluated on a grid containing 0
  pts <- c(0, tau_grid)
  Cmat <- KL_cov_matrix(pts, KL, L)
  approx_rho <- Cmat[1, -1]
  
  # plot approx for several orders
  plot(tau_grid, true_corr(tau_grid), type = "l", lwd = 2, col = "black",
       ylim = c(-0.05, 1), xlab = expression(tau), ylab = expression(rho(tau)),
       main = paste("K-L truncation: 2N =", 2*N, "modes"))
  lines(tau_grid, approx_rho, col = "red", lwd = 2, lty = 2)
  legend("topright", c("True exp(-tau)", "K-L approx"),
         col = c("black", "red"), lty = c(1, 2), lwd = 2, bty = "n")
}

###############################################################################
# P 5: form slides - Page 13 approximation  lambda_j ~ f(j pi/(2L))
# Spectral density of exponential correlation in 1D:
#   f(k) = (1/pi) * phi / (k^2 + phi^2)
# Eigenfunctions: psi_j(t) ~ c * exp(i j pi t / (2L))
###############################################################################

f_exp_spectrum <- function(k, phi = 1) (1/pi) * phi / (k^2 + phi^2)

KL_page13_corr <- function(tau, N, L, phi = 1) {
  # Slide 13 (LEC5): lambda_j = f(j*pi/(2L)), psi_j(t) = (1/sqrt(2L)) * exp(i*w_j*t)
  w <- (1:N) * pi / (2 * L)
  lam <- f_exp_spectrum(w, phi)
  
  # rho(tau) = sum_j lambda_j * psi_j(0) * psi_j(tau) = sum_j lambda_j * cos(w_j*tau) / (2L)
  approx <- sapply(tau, function(t) sum(lam * cos(w * t)) / (2 * L))
  
  # Normalize so rho(0) = 1
  approx / (sum(lam) / (2 * L))
}


par(mfrow = c(2, 2))
for (N in orders) {
  approx_rho <- KL_page13_corr(tau_grid, N, L, phi)
  plot(tau_grid, true_corr(tau_grid), type = "l", lwd = 2, col = "black",
       ylim = c(-0.1, 1), xlab = expression(tau), ylab = expression(rho(tau)),
       main = paste("stationary correlation function -  spectral approx: N =", N))
  lines(tau_grid, approx_rho, col = "blue", lwd = 2, lty = 2)
  legend("topright", c("True", "Spectral approx"),
         col = c("black", "blue"), lty = c(1, 2), lwd = 2, bty = "n")
}

###############################################################################
# P6: Simulate 100 realizations of GP with exponential correlation,
# compare empirical eigenvalues/eigenfunctions to K-L and spectral approximation
###############################################################################

# Discretize [-L, L]
m <- 100
s_disc <- seq(-L, L, length.out = m)
ds <- s_disc[2] - s_disc[1]

# True covariance matrix (sigma^2 = 1)
D <- as.matrix(dist(s_disc))
C_true <- exp(-phi * D)

# Cholesky for simulation
R <- chol(C_true + 1e-8 * diag(m))

# Generate 100 realizations
n_real <- 100
Z <- matrix(rnorm(n_real * m), nrow = m, ncol = n_real)
X <- t(R) %*% Z  # m x n_real, each column is a realization

# Empirical covariance
C_hat <- cov(t(X))   # m x m

# Empirical eigendecomposition
# Continuous K-L:  int C(s,s') psi(s') ds' = lambda psi(s)
# Discrete approximation: (C * ds) v = lambda v
eig_emp <- eigen(C_hat * ds, symmetric = TRUE)
lam_emp <- eig_emp$values

# Eigenfunctions normalized so int psi^2 ds = 1  =>  sum(psi^2) * ds = 1
psi_emp <- eig_emp$vectors / sqrt(ds)

# True K-L eigenvalues (analytic)
KL_full <- build_KL(20, L, phi)
lam_KL <- sort(c(KL_full$lam1, KL_full$lam2), decreasing = TRUE)

# Page-13 spectral approx eigenvalues
N_spec <- 20
w_spec <- (1:N_spec) * pi / (2 * L)
lam_spec <- sort(f_exp_spectrum(w_spec, phi) * 2 * pi, decreasing = TRUE) # factor, but depends on normalization convention 

# Plot top 15 eigenvalues
par(mfrow = c(1, 1))
K <- 15
plot(1:K, lam_emp[1:K], type = "b", pch = 19, col = "black",
     ylim = range(c(lam_emp[1:K], lam_KL[1:K], lam_spec[1:K])),
     xlab = "Index j", ylab = expression(lambda[j]),
     main = "Eigenvalue comparison: empirical vs K-L vs spectral")
lines(1:K, lam_KL[1:K], type = "b", pch = 17, col = "red")
lines(1:K, lam_spec[1:K], type = "b", pch = 15, col = "blue")
legend("topright",
       c("Empirical (100 reals)", "Analytic K-L", "Page-13 spectral"),
       col = c("black", "red", "blue"), pch = c(19, 17, 15), lty = 1, bty = "n")

# Compare top 4 eigenfunctions
# Build analytic K-L eigenfunctions on the same grid, ordered by eigenvalue
build_KL_eigfuns <- function(s_pts, KL, L) {
  w1 <- KL$w1; w2 <- KL$w2
  lam1 <- KL$lam1; lam2 <- KL$lam2
  norm1 <- sqrt(L + sin(2 * w1 * L) / (2 * w1))
  norm2 <- sqrt(L - sin(2 * w2 * L) / (2 * w2))
  
  funs <- cbind(
    sapply(seq_along(w1), function(j) cos(w1[j] * s_pts) / norm1[j]),
    sapply(seq_along(w2), function(j) sin(w2[j] * s_pts) / norm2[j])
  )
  lams <- c(lam1, lam2)
  ord <- order(lams, decreasing = TRUE)
  list(funs = funs[, ord], lams = lams[ord])
}

KL_funs 

<- build_KL_eigfuns(s_disc, KL_full, L)

par(mfrow = c(2, 2))
for (j in 1:4) {
  # plot eigenvalues
  emp <- psi_emp[, j]
  ana <- KL_funs$funs[, j]
  if (sum(emp * ana) < 0) emp <- -emp
  
  plot(s_disc, ana, type = "l", lwd = 2, col = "red",
       xlab = "s", ylab = expression(psi(s)),
       main = paste("Eigenfunction j =", j))
  lines(s_disc, emp, col = "black", lwd = 2, lty = 2)
  legend("topright", c("Analytic K-L", "Empirical"),
         col = c("red", "black"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
}