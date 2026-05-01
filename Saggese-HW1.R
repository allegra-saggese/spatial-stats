library(geoR)

# ============================================================
# Q3: COVARIOGRAMS AND VARIOGRAMS
# sigma^2 = 1, phi calibrated so correlation = 0.05, tau = 1 (except for powered family)
# ============================================================

sigma2 <- 1
tau  <- seq(0, 2, length.out = 500)  # distance grid

# PART 1: VAR 
# 1. Powered Exponential: C(tau) = exp(-(tau/phi)^nu), nu=1 (exponential)
# Solve: exp(-1/phi) = 0.05 => phi = -1/log(0.05)

nu_vals <- c(0.5, 1, 1.5, 2)

plot(tau, tau*0, type="n", ylim=c(0,1.1),
     xlab=expression(tau), ylab="C(tau)",
     main="Powered Exponential family")

for(nu in nu_vals){
  phi_nu <- 1/(-log(0.05))^(1/nu)
  lines(tau, exp(-(tau/phi_nu)^nu), lwd=2,
        col=which(nu_vals==nu))
}

legend("topright",
       legend=paste("nu =", nu_vals),
       col=1:4, lwd=2)

# 2. Spherical: C(tau) = 1 - 3tau/2phi + tau^3/2phi^3 for tau <= phi, 0 otherwise
# Solve numerically: need phi >= 1
f_sph <- function(phi) {
  if(phi < 1) return(Inf)  # phi must be >= 1 for C(1) to be defined
  1 - 1.5*(1/phi) + 0.5*(1/phi)^3 - 0.05
}
phi_sph <- uniroot(f_sph, c(1, 100))$root
1 - 1.5*(1/phi_sph) + 0.5*(1/phi_sph)^3 # should give .05 for correlation 


# 3. Rational Quadratic: C(tau) = 1 - tau^2/(phi^2 + tau^2)
# Solved from plugging in s.t. 1 - 1/(phi^2+1) = 0.05 => phi = sqrt(1/0.95 - 1)
phi_rq <- sqrt(1/0.95 - 1)
1 - 1/(phi_rq^2 + 1) # should give corr

# 4. Wave (sine): C(tau) = sin(tau/phi)/(tau/phi)
f_wave <- function(phi) sin(1/phi)/(1/phi) - 0.05
f_wave(0.1)  # should be negative
f_wave(2)
phi_wave <- uniroot(f_wave, c(0.01, 2))$root

# 5. Matern: using geoR's matern function, nu=1.5
f_mat <- function(phi) {
  geoR::matern(1, phi=phi, kappa=1.5) - 0.05
}
phi_mat <- uniroot(f_mat, c(0.01, 100))$root
phi_mat
(1 + 1/phi_mat)*exp(-1/phi_mat)

# PART 2: covvar functions
C_powexp <- function(tau, phi) exp(-(tau/phi))
C_sph    <- function(tau, phi) ifelse(tau <= phi,
                                      1 - 1.5*(tau/phi) + 0.5*(tau/phi)^3, 0)
C_rq     <- function(tau, phi) 1 - tau^2/(phi^2 + tau^2)
C_wave   <- function(tau, phi) ifelse(tau==0, 1, sin(tau/phi)/(tau/phi))
C_mat    <- function(tau, phi) geoR::matern(tau, phi=phi, kappa=1.5)


# variogram functions: gamma(tau) = sigma^2 - C(tau) 
G_powexp <- function(tau, phi) sigma2*(1 - C_powexp(tau, phi))
G_sph    <- function(tau, phi) sigma2*(1 - C_sph(tau, phi))
G_rq     <- function(tau, phi) sigma2*(1 - C_rq(tau, phi))
G_wave   <- function(tau, phi) sigma2*(1 - C_wave(tau, phi))
G_mat    <- function(tau, phi) sigma2*(1 - C_mat(tau, phi))

# plot covariograms
par(mfrow=c(1,2), mar=c(4,4,2,1))

plot(tau, sigma2*C_powexp(tau, phi_powexp), type="l", col=1, lwd=2,
     ylim=c(-0.1,1.1), xlab=expression(tau), ylab="C(tau)",
     main="Covariograms")
lines(tau, sigma2*C_sph(tau, phi_sph),  col=2, lwd=2)
lines(tau, sigma2*C_rq(tau, phi_rq),    col=3, lwd=2)
lines(tau, sigma2*C_wave(tau, phi_wave), col=4, lwd=2)
lines(tau, sigma2*C_mat(tau, phi_mat),  col=5, lwd=2)
abline(h=0.05, lty=2, col="gray")
abline(v=1,    lty=2, col="gray")
legend("topright", 
       legend=c("Pow.Exp","Spherical","Rat.Quad","Wave","Matern(1.5)"),
       col=1:5, lwd=2, cex=0.7)

# plot variograms
plot(tau, G_powexp(tau, phi_powexp), type="l", col=1, lwd=2,
     ylim=c(0,1.2), xlab=expression(tau), ylab=expression(gamma(tau)),
     main="Variograms")
lines(tau, G_sph(tau, phi_sph),  col=2, lwd=2)
lines(tau, G_rq(tau, phi_rq),    col=3, lwd=2)
lines(tau, G_wave(tau, phi_wave), col=4, lwd=2)
lines(tau, G_mat(tau, phi_mat),  col=5, lwd=2)
abline(h=sigma2, lty=2, col="gray")
legend("bottomright",
       legend=c("Pow.Exp","Spherical","Rat.Quad","Wave","Matern(1.5)"),
       col=1:5, lwd=2, cex=0.7)


# ============================================================
# Q4: SIMULATE 1D GP REALIZATIONS VIA CHOLESKY
# 100 points on [0,1], one realization per covariance model
# ============================================================

set.seed(123)
n    <- 100
locs <- seq(0, 1, length.out=n)  # 1D locations
D    <- as.matrix(dist(locs))     # n x n distance matrix

# Function: simulate GP realization given covariance matrix
sim_gp <- function(Sigma) {
  # Cholesky decomposition: Sigma = L L'
  # X = L * Z, Z ~ N(0,I) => X ~ N(0, Sigma)
  L <- chol(Sigma + 1e-10*diag(n))  # small nugget for numerical stability
  Z <- rnorm(n)
  as.vector(t(L) %*% Z)
}

# Build covariance matrices
Sigma_powexp <- sigma2 * C_powexp(D, phi_powexp)
Sigma_sph    <- sigma2 * C_sph(D, phi_sph)
Sigma_rq     <- sigma2 * C_rq(D, phi_rq)
Sigma_wave   <- sigma2 * C_wave(D, phi_wave)
Sigma_mat    <- sigma2 * C_mat(D, phi_mat)

# Simulate
X_powexp <- sim_gp(Sigma_powexp)
X_sph    <- sim_gp(Sigma_sph)
X_rq     <- sim_gp(Sigma_rq)
X_wave   <- sim_gp(Sigma_wave)
X_mat    <- sim_gp(Sigma_mat)


# Plot realizations
par(mfrow=c(3,2), mar=c(4,4,2,1))

plot(locs, X_powexp, type="l", col=1, lwd=1.5,
     main="Powered Exponential", xlab="s", ylab="X(s)")
plot(locs, X_sph, type="l", col=2, lwd=1.5,
     main="Spherical", xlab="s", ylab="X(s)")
plot(locs, X_rq, type="l", col=3, lwd=1.5,
     main="Rational Quadratic", xlab="s", ylab="X(s)")
plot(locs, X_wave, type="l", col=4, lwd=1.5,
     main="Wave", xlab="s", ylab="X(s)")
plot(locs, X_mat, type="l", col=5, lwd=1.5,
     main="Matern (nu=1.5)", xlab="s", ylab="X(s)")