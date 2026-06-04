# 02b_trend_geoR.R
# Fit spatial trend using geoR likfit() — GLS estimates of trend parameters
# accounting for spatial correlation (vs OLS from lm() in 02_eda.R)
# Uses best-fitting covariance model: Matern kappa=1.0, trend="2nd"

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

# ── Fit null / 1st / 2nd order trend with spatial covariance ─────────────────
# kappa=1.0, exponential-like Matern — best fit from step 05
fit_cte <- likfit(gd, trend = "cte", ini.cov.pars = c(0.20, 40),
                  nugget = 0.03, cov.model = "matern", kappa = 1,
                  lik.method = "ML", messages = FALSE)

fit_1st <- likfit(gd, trend = "1st", ini.cov.pars = c(0.20, 40),
                  nugget = 0.03, cov.model = "matern", kappa = 1,
                  lik.method = "ML", messages = FALSE)

fit_2nd <- likfit(gd, trend = "2nd", ini.cov.pars = c(0.20, 40),
                  nugget = 0.03, cov.model = "matern", kappa = 1,
                  lik.method = "ML", messages = FALSE)

# ── Likelihood ratio tests ────────────────────────────────────────────────────
# LRT: -2 * (logLik_reduced - logLik_full) ~ chi-sq(df = extra params)
lrt_1st <- -2 * (fit_cte$loglik - fit_1st$loglik)  # 2 extra params (x, y)
lrt_2nd <- -2 * (fit_1st$loglik - fit_2nd$loglik)  # 3 extra params (x^2, y^2, xy)

p_1st <- pchisq(lrt_1st, df = 2, lower.tail = FALSE)
p_2nd <- pchisq(lrt_2nd, df = 3, lower.tail = FALSE)

cat("── Likelihood ratio tests (spatial GLS) ─────────────────────\n")
cat(sprintf("1st order vs null :  LRT = %.3f  df = 2  p = %.4f\n", lrt_1st, p_1st))
cat(sprintf("2nd order vs 1st  :  LRT = %.3f  df = 3  p = %.4f\n", lrt_2nd, p_2nd))

cat("\n── Model log-likelihoods and AIC ────────────────────────────\n")
models <- list("Null (constant)" = fit_cte,
               "1st order"       = fit_1st,
               "2nd order"       = fit_2nd)
for (nm in names(models)) {
  f <- models[[nm]]
  cat(sprintf("  %-18s  logLik = %7.2f   AIC = %7.2f\n",
              nm, f$loglik, -2*f$loglik + 2*f$npars))
}

# ── Beta coefficients from 2nd order fit ─────────────────────────────────────
# beta.vcov gives the GLS covariance matrix of the trend estimates
beta  <- fit_2nd$beta
se    <- sqrt(diag(fit_2nd$beta.vcov))
tstat <- beta / se
pval  <- 2 * pt(abs(tstat), df = length(gd$data) - length(beta), lower.tail = FALSE)

param_names <- c("Intercept", "x (E-W km)", "y (N-S km)",
                 "x^2", "y^2", "x*y")

cat("\n── 2nd order trend coefficients (GLS, kappa=1.0) ────────────\n")
coef_tbl <- data.frame(
  term    = param_names,
  beta    = round(beta,  5),
  SE      = round(se,    5),
  t       = round(tstat, 3),
  p_value = round(pval,  4)
)
print(coef_tbl, row.names = FALSE)

cat("\n── Covariance parameters at 2nd order fit ───────────────────\n")
cat(sprintf("  sigma^2 (partial sill) = %.4f\n", fit_2nd$cov.pars[1]))
cat(sprintf("  phi     (range, km)    = %.2f\n",  fit_2nd$cov.pars[2]))
cat(sprintf("  tau^2   (nugget)       = %.4f\n",  fit_2nd$nugget))
