## simulate_data.R
##
## Implements the data-generating process described in Section 3.1 of the
## manuscript: CVOT-representative baseline covariates, a logistic
## non-susceptibility (incidence) model, and Weibull-PH event times among
## susceptibles. Supports both the PH-S (primary) and NPH-S (robustness)
## settings.
##
## The manuscript gives the *form* of the DGP and the covariate distributions,
## but not the exact regression coefficients used to produce the published
## table (they are not reported). The coefficients below are documented,
## sensible choices calibrated to hit the three headline Su(tau) targets
## (~0.02 well-identified, ~0.80 REWIND-calibrated, ~0.95 insufficient
## follow-up) via root-finding on the incidence/latency intercepts, exactly as
## Su(tau) is described functioning as the identifiability diagnostic in the
## paper. This reproduces the qualitative findings (Cox attenuation bias
## growing with the non-susceptible fraction; PHMC incidence robustness;
## PHMC latency degradation at high Su(tau)) rather than the literal
## (unreported) numeric table.

#' Simulate one CVOT-like trial dataset under the mixture-cure DGP
#'
#' @param n total sample size (n/2 per arm)
#' @param alpha0 incidence-model intercept (controls baseline cure fraction)
#' @param alpha1 incidence-model treatment log-odds-ratio (TRUE incidence effect)
#' @param alpha2 numeric vector of length 5, incidence-model covariate coefficients
#'   for (age, diabDur, sex, BMI, sbp), on the standardized covariate scale
#' @param gamma0 latency-model (Weibull) log-scale intercept
#' @param beta1 latency-model treatment log-HR (TRUE latency effect)
#' @param beta2 numeric vector of length 5, latency-model covariate log-HRs
#' @param weibull_shape Weibull shape parameter k for the susceptible hazard
#' @param tau administrative follow-up horizon (years)
#' @param admin_dropout_rate rate of an additional independent
#'   (random, non-informative) censoring process representing loss to follow-up
#' @param nph if TRUE, use a non-proportional-hazards-among-susceptibles (NPH-S)
#'   generating mechanism: the treatment log-HR on the susceptible hazard
#'   changes over time (crossing hazards), per Section 3.1
#' @param nph_beta1_late treatment log-HR after `nph_crossing_time` under NPH-S
#' @param nph_crossing_time time (years) at which the susceptible-hazard
#'   treatment effect changes, under NPH-S
#'
#' @return data.frame with columns: trt, age, diabDur, sex, bmi, sbp,
#'   time (observed, possibly censored), event, C_true (1=non-susceptible),
#'   T_true (latent true event time, Inf if non-susceptible)
simulate_phmc_trial <- function(n,
                                 alpha0, alpha1,
                                 alpha2 = c(age = 0.15, diabDur = 0.05, sex = -0.10, bmi = 0.10, sbp = 0.05),
                                 gamma0, beta1,
                                 beta2 = c(age = 0.20, diabDur = 0.05, sex = 0.05, bmi = 0.05, sbp = 0.10),
                                 weibull_shape = 1.2,
                                 tau,
                                 admin_dropout_rate = 0.01,
                                 nph = FALSE,
                                 nph_beta1_late = NULL,
                                 nph_crossing_time = NULL) {

  ## ---- baseline covariates (Section 3.1) ----
  trt     <- rbinom(n, 1, 0.5)
  age     <- rnorm(n, 63, 8)
  diabDur <- rgamma(n, shape = 2.5, rate = 0.2)
  sex     <- rbinom(n, 1, 0.54)
  bmi     <- rnorm(n, 32, 5)
  sbp     <- rnorm(n, 137, 16)

  # standardize continuous covariates for numerically stable, interpretable
  # coefficients (coefficients above are on this standardized scale)
  z_age     <- (age - 63) / 8
  z_diabDur <- (diabDur - mean(diabDur)) / sd(diabDur)
  z_sex     <- sex - 0.54
  z_bmi     <- (bmi - 32) / 5
  z_sbp     <- (sbp - 137) / 16
  Xz <- cbind(age = z_age, diabDur = z_diabDur, sex = z_sex, bmi = z_bmi, sbp = z_sbp)

  ## ---- incidence (susceptibility) model: logit P(C=1 | trt, x) ----
  lp_inc <- alpha0 + alpha1 * trt + as.numeric(Xz %*% alpha2)
  p_cure <- plogis(lp_inc)
  C_true <- rbinom(n, 1, p_cure)     # 1 = non-susceptible ("cured")

  ## ---- latency model: Weibull PH hazard among susceptibles ----
  if (!nph) {
    lp_lat <- gamma0 + beta1 * trt + as.numeric(Xz %*% beta2)
    U <- runif(n)
    T_true <- (-log(U) / exp(lp_lat))^(1 / weibull_shape)
  } else {
    stopifnot(!is.null(nph_beta1_late), !is.null(nph_crossing_time))
    # simulate via piecewise-constant treatment log-HR (crossing hazards):
    # generate from the "early" regime hazard, then re-map draws whose
    # implied time exceeds the crossing time using the "late" regime hazard,
    # preserving a valid (continuous) piecewise Weibull-PH survival function.
    lp_lat_early <- gamma0 + beta1 * trt + as.numeric(Xz %*% beta2)
    lp_lat_late  <- gamma0 + nph_beta1_late * trt + as.numeric(Xz %*% beta2)
    tc <- nph_crossing_time
    U <- runif(n)
    # cumulative hazard to time t: H(t) = exp(lp_early) t^k  for t<=tc
    #                              H(t) = exp(lp_early) tc^k + exp(lp_late)(t^k - tc^k) for t>tc
    H_at_tc <- exp(lp_lat_early) * tc^weibull_shape
    S_at_tc <- exp(-H_at_tc)
    T_true <- numeric(n)
    # Invert the piecewise cumulative hazard: draw target cumulative hazard
    # Etar ~ Exp(1) and solve H(T_true) = Etar piecewise.
    Etar <- -log(U)
    is_early <- Etar <= H_at_tc
    T_true[is_early] <- (Etar[is_early] / exp(lp_lat_early[is_early]))^(1 / weibull_shape)
    rem <- Etar[!is_early] - H_at_tc[!is_early]
    T_true[!is_early] <- (tc^weibull_shape + rem / exp(lp_lat_late[!is_early]))^(1 / weibull_shape)
  }
  T_true[C_true == 1] <- Inf   # non-susceptibles never experience the event

  ## ---- censoring: administrative at tau + light independent dropout ----
  C_dropout <- rexp(n, rate = admin_dropout_rate)
  C_obs <- pmin(tau, C_dropout)

  time  <- pmin(T_true, C_obs)
  event <- as.integer(T_true <= C_obs)

  data.frame(trt = trt, age = age, diabDur = diabDur, sex = sex, bmi = bmi, sbp = sbp,
             time = time, event = event, C_true = C_true, T_true = T_true)
}

#' Root-find the incidence intercept alpha0 that yields a target marginal
#' cure fraction (averaged over both arms, at the design's covariate
#' distribution), holding all other DGP parameters fixed.
calibrate_alpha0 <- function(target_cure, alpha1, alpha2, n_cal = 20000, ...) {
  f <- function(a0) {
    d <- simulate_phmc_trial(n = n_cal, alpha0 = a0, alpha1 = alpha1, alpha2 = alpha2,
                              gamma0 = 0, beta1 = 0, tau = 1e6, admin_dropout_rate = 0, ...)
    mean(d$C_true) - target_cure
  }
  uniroot(f, interval = c(-6, 6), tol = 1e-3)$root
}

#' Root-find the Weibull log-scale intercept gamma0 that yields a target
#' Su(tau): P(T_true > tau | susceptible), at a given follow-up horizon tau,
#' holding shape/covariate coefficients fixed. This is the calibration
#' target used to construct the Su(tau) identifiability gradient across the
#' 13 simulation scenarios (Section 3.1, Table 4).
calibrate_gamma0 <- function(target_Su_tau, tau, beta1, beta2, weibull_shape,
                              alpha0, alpha1, alpha2, n_cal = 20000) {
  f <- function(g0) {
    d <- simulate_phmc_trial(n = n_cal, alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
                              gamma0 = g0, beta1 = beta1, beta2 = beta2,
                              weibull_shape = weibull_shape, tau = tau, admin_dropout_rate = 0)
    ds <- d[d$C_true == 0, ]
    mean(ds$T_true > tau) - target_Su_tau
  }
  uniroot(f, interval = c(-6, 6), tol = 1e-3)$root
}
