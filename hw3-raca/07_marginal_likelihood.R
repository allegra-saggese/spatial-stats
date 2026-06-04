# 07_marginal_likelihood.R
# Marginal likelihood for the range parameter phi, integrating out
# trend parameters beta and variance sigma^2 analytically.
#
# Model: Y ~ N(F*beta, sigma^2 * V_tilde)
#   V_tilde = R(phi, kappa) + delta*I,  delta = tau^2/sigma^2 (plugged in)
#
# Under flat prior on beta and Jeffreys prior 1/sigma^2 on sigma^2:
#
#   l_marg(phi) = -1/2 log|V_tilde|
#               - 1/2 log|F' V_tilde^{-1} F|
#               - (n-p)/2 * log(S^2(phi))
#
# where S^2(phi) = (Y - F*beta_hat)' V_tilde^{-1} (Y - F*beta_hat)
# and beta_hat = GLS estimate at each phi.

library(readr)
library(dplyr)
library(geoR)

va   <- read_csv("VA_soil_carbon_clean.csv", show_col_types = FALSE)
va_m <- filter(va, soc_measured == "measured") |> mutate(lsoc = log(soc))

ref_lon <- mean(va_m$lon); ref_lat <- mean(va_m$lat)
va_m <- va_m |>
  mutate(x_km = (lon - ref_lon) * cos(ref_lat * pi / 180) * 111.32,
         y_km = (lat - ref_lat) * 111.32)

set.seed(42)
gd <- as.geodata(va_m, coords.col = c("x_km", "y_km"), data.col = "lsoc")
gd <- jitterDupCoords(gd, max = 0.5)

# ── Setup ──────────────────────────────────────────────────────────────────────
Y      <- gd$data
n      <- length(Y)
D      <- as.matrix(dist(gd$coords))
coords <- gd$coords

F_mat <- cbind(1,
               coords[,1], coords[,2],
               coords[,1]^2, coords[,2]^2,
               coords[,1] * coords[,2])
p <- ncol(F_mat)   # 6 trend parameters

# ML estimates from step 06 — used to plug in delta = tau^2 / sigma^2
# kappa: 0.5   sigma^2=0.556  tau^2=0.027  delta=0.049
# kappa: 1.0   sigma^2=0.608  tau^2=0.033  delta=0.055
# kappa: 1.5   sigma^2=0.325  tau^2=0.212  delta=0.652
# kappa: 2.5   sigma^2=0.335  tau^2=0.212  delta=0.633
kappas <- c(0.5,   1.0,   1.5,   2.5)
deltas <- c(0.0272/0.5564, 0.0333/0.6079, 0.212/0.3252, 0.212/0.3347)

cat("Relative nuggets (delta = tau^2/sigma^2) plugged in:\n")
for (i in seq_along(kappas))
  cat(sprintf("  kappa=%.1f  delta=%.4f\n", kappas[i], deltas[i]))

# ── Marginal log-likelihood function ──────────────────────────────────────────
marg_loglik <- function(phi, kappa, delta, D, Y, F_mat) {
  n <- nrow(D); p <- ncol(F_mat)

  R <- matrix(matern(as.vector(D), phi = phi, kappa = kappa), n, n)
  diag(R) <- 1
  V <- R + delta * diag(n)          # V_tilde

  L <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(L)) return(NA_real_)

  log_det_V <- 2 * sum(log(diag(L)))

  Li_F <- backsolve(L, F_mat, transpose = TRUE)   # L^{-T} F
  Li_Y <- backsolve(L, Y,     transpose = TRUE)   # L^{-T} Y

  FtVF     <- crossprod(Li_F)                     # F' V^{-1} F
  LF       <- chol(FtVF)
  log_det_FtVF <- 2 * sum(log(diag(LF)))

  beta_hat <- backsolve(LF, backsolve(LF, crossprod(Li_F, Li_Y), transpose = TRUE))
  resid    <- Li_Y - Li_F %*% beta_hat
  S2       <- sum(resid^2)

  if (S2 <= 0) return(NA_real_)

  -0.5 * log_det_V - 0.5 * log_det_FtVF - (n - p) / 2 * log(S2)
}

# ── Evaluate over phi grid ─────────────────────────────────────────────────────
phi_grid <- seq(2, 200, length.out = 200)
cols     <- c("#2166ac", "#4dac26", "#d6604d", "#8B1A1A")
ltys     <- 1:4

cat("\nEvaluating marginal likelihood over phi grid...\n")
marg_ll <- lapply(seq_along(kappas), function(i) {
  cat(sprintf("  kappa = %.1f\n", kappas[i]))
  vals <- sapply(phi_grid, marg_loglik,
                 kappa = kappas[i], delta = deltas[i],
                 D = D, Y = Y, F_mat = F_mat)
  vals - max(vals, na.rm = TRUE)   # normalise so peak = 0
})

# ── Plot ────────────────────────────────────────────────────────────────────────
png("figures/16_marginal_likelihood_range.png", width = 820, height = 560, res = 120)

ylim <- range(unlist(marg_ll), na.rm = TRUE)
ylim[1] <- max(ylim[1], -10)

plot(NA, xlim = range(phi_grid), ylim = ylim,
     main = "Marginal log-likelihood for range φ  (β and σ² integrated out)",
     xlab = "Range φ (km)",
     ylab = "Marginal log-likelihood (relative to maximum)")

abline(h = -1.92, lty = 3, col = "grey55")
text(195, -1.92 + 0.25, "95% CI", col = "grey55", adj = 1, cex = 0.8)
abline(h = 0, lty = 1, col = "grey85")

for (i in seq_along(kappas))
  lines(phi_grid, marg_ll[[i]], col = cols[i], lwd = 2, lty = ltys[i])

# Mark MLEs
for (i in seq_along(kappas)) {
  phi_mle <- phi_grid[which.max(marg_ll[[i]])]
  points(phi_mle, 0, pch = 4, cex = 1.2, col = cols[i], lwd = 2)
}

legend("topright", bty = "n", lwd = 2, lty = ltys, col = cols,
       legend = sprintf("κ = %.1f   δ = %.3f", kappas, deltas))
dev.off()
cat("Saved: figures/16_marginal_likelihood_range.png\n")

# ── Report MLEs and 95% CIs ────────────────────────────────────────────────────
cat("\n── Marginal likelihood summary ───────────────────────────────\n")
for (i in seq_along(kappas)) {
  ml  <- marg_ll[[i]]
  mle <- phi_grid[which.max(ml)]
  ci  <- range(phi_grid[!is.na(ml) & ml > -1.92])
  cat(sprintf("kappa=%.1f  phi MLE=%.1f km  95%%CI=[%.1f, %.1f] km\n",
              kappas[i], mle, ci[1], ci[2]))
}
