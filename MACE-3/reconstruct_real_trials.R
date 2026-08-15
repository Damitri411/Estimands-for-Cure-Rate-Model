## reconstruct_real_trials.R
##
## Calibrates a synthetic, mixture-cure-consistent individual-patient dataset
## per trial so that its OBSERVABLE aggregate quantities -- per-arm sample
## size, per-arm event rate, overall marginal Cox HR, and median follow-up --
## match the real, published summary statistics in R/real_trial_stats.R as
## closely as numerical calibration allows. The Cox/PHMC pipeline
## (R/fit_wrapper.R) is then run on this calibrated synthetic dataset to
## produce Table 5/6- and Figure 1/4-style output.
##
## IMPORTANT: this is a reconstruction, not a re-analysis. No individual
## patient-level data from SURPASS-CVOT or REWIND is used or was available;
## only the aggregate numbers reported in the publications are used as
## calibration targets. The incidence/latency SPLIT of any given trial's
## marginal HR is not identified by the reported aggregate HR alone (that
## split is exactly the estimand the manuscript's PHMC model targets) -- so
## the calibration below picks a plausible, documented incidence/latency
## split (roughly matching the manuscript's own qualitative finding that
## effects concentrate in incidence rather than latency) as a STARTING
## assumption, then numerically adjusts it so the synthetic data's fitted
## marginal Cox HR matches the reported HR. Treat the resulting PHMC
## incidence/latency decomposition as illustrative of the method applied to
## realistic aggregate targets, not as the trials' true decomposition.
##
## Depends on: R/simulate_data.R, R/fit_wrapper.R, R/real_trial_stats.R
## (source those first, e.g. via run_real_trials.R).

#' Calibrate a synthetic dataset to a trial's real aggregate summary stats
#'
#' @param stats one element of REAL_TRIAL_STATS (e.g. REAL_TRIAL_STATS[["REWIND"]])
#' @param n_multiplier simulate at n_multiplier x the real trial's sample
#'   size then subsample down, OR simulate at exactly the real N (default 1).
#'   Kept as an option in case you want a larger calibration sample for
#'   smoother root-finding; final dataset is always trimmed to the real N.
#' @param assumed_incidence_share in [0,1]: how much of the reported
#'   log(1/HR) treatment benefit to initially assign to the incidence
#'   (cure-probability) pathway vs. the latency (event-timing) pathway,
#'   before the joint calibration step corrects the marginal Cox HR to match
#'   the reported value. Default 0.7, reflecting the manuscript's own
#'   qualitative finding (effect concentrated in incidence) as a documented,
#'   labeled starting assumption -- see the file-level note above.
#' @return list(data = synthetic IPD data.frame, params = calibrated DGP
#'   parameters, achieved = achieved aggregate stats for comparison to target)
calibrate_trial <- function(stats, n_multiplier = 1, assumed_incidence_share = 0.7,
                             weibull_shape = 1.3, cure_baseline = 0.15, seed = 1) {
  set.seed(seed)
  n_total <- round((stats$n_control + stats$n_treated) * n_multiplier)
  ## NOTE: this simplified DGP has no staggered accrual (everyone "enters" at
  ## time 0, administrative censoring at a single fixed tau) -- unlike a real
  ## trial, where staggered enrollment makes median follow-up shorter than
  ## the total trial duration. With event rates under 50% here, most
  ## simulated subjects are administratively censored at tau, so the
  ## simulated median observed time lands at tau itself; setting tau equal
  ## to the reported median follow-up directly reproduces that target.
  tau <- stats$median_followup_years

  alpha2 <- c(age = 0.10, diabDur = 0.05, sex = -0.05, bmi = 0.08, sbp = 0.03)
  beta2  <- c(age = 0.15, diabDur = 0.05, sex = 0.03,  bmi = 0.05, sbp = 0.05)

  ## Step 1: calibrate the CONTROL-arm hazard (gamma0) and incidence intercept
  ## (alpha0) to reproduce the control arm's reported event rate at tau,
  ## for the assumed baseline cure fraction.
  alpha0 <- calibrate_alpha0(target_cure = cure_baseline, alpha1 = 0, alpha2 = alpha2)

  target_event_rate_control <- stats$event_rate_control
  g0_obj <- function(g0) {
    set.seed(90001)  # fixed seed -> deterministic function of g0 (removes MC noise from uniroot's search)
    d <- simulate_phmc_trial(n = 20000, alpha0 = alpha0, alpha1 = 0, alpha2 = alpha2,
                              gamma0 = g0, beta1 = 0, beta2 = beta2,
                              weibull_shape = weibull_shape, tau = tau, admin_dropout_rate = 0)
    d0 <- d[d$trt == 0, ]
    mean(d0$event) - target_event_rate_control
  }
  gamma0 <- uniroot(g0_obj, interval = c(-8, 4), tol = 1e-3)$root

  ## Step 2: jointly calibrate (alpha1, beta1) -- the incidence and latency
  ## treatment effects -- so that (a) the treated-arm event rate matches the
  ## reported value and (b) the overall marginal Cox HR (fit on the full,
  ## two-arm synthetic dataset) matches the reported HR. Two targets, two
  ## free parameters -> solved via a 2D optimizer starting from the
  ## documented incidence/latency split assumption above.
  log_total_effect <- log(stats$hr_reported)  # negative if treatment is beneficial
  a1_start <- assumed_incidence_share * (-log_total_effect)  # incidence: benefit -> MORE cure -> positive alpha1
  b1_start <- (1 - assumed_incidence_share) * log_total_effect  # latency: benefit -> lower hazard -> negative beta1

  discrepancy <- function(par) {
    a1 <- par[1]; b1 <- par[2]
    set.seed(90002)  # fixed seed -> deterministic function of (a1,b1) (removes MC noise from optim's search)
    d <- simulate_phmc_trial(n = 120000, alpha0 = alpha0, alpha1 = a1, alpha2 = alpha2,
                              gamma0 = gamma0, beta1 = b1, beta2 = beta2,
                              weibull_shape = weibull_shape, tau = tau, admin_dropout_rate = 0)
    d1 <- d[d$trt == 1, ]
    er1 <- mean(d1$event)
    cx <- tryCatch(
      survival::coxph(survival::Surv(time, event) ~ trt, data = d),
      error = function(e) NULL
    )
    hr_hat <- if (is.null(cx)) NA else exp(coef(cx)["trt"])
    if (is.na(hr_hat)) return(1e3)
    (er1 - stats$event_rate_treated)^2 * 400 + (log(hr_hat) - log_total_effect)^2 * 400
  }
  opt <- optim(par = c(a1_start, b1_start), fn = discrepancy, method = "Nelder-Mead",
               control = list(reltol = 1e-9, maxit = 300))
  alpha1 <- opt$par[1]; beta1 <- opt$par[2]

  ## Step 3: final large-sample simulation at the (approximately) real N,
  ## trimmed to exactly match the real per-arm sample sizes.
  set.seed(90003 + seed)
  d_full <- simulate_phmc_trial(n = max(n_total * 3, 20000), alpha0 = alpha0, alpha1 = alpha1,
                                 alpha2 = alpha2, gamma0 = gamma0, beta1 = beta1, beta2 = beta2,
                                 weibull_shape = weibull_shape, tau = tau, admin_dropout_rate = 0)
  d0 <- d_full[d_full$trt == 0, ][seq_len(stats$n_control), ]
  d1 <- d_full[d_full$trt == 1, ][seq_len(stats$n_treated), ]
  d <- rbind(d0, d1)
  d$arm <- factor(d$trt, levels = c(0, 1), labels = stats$arm_labels[c("control", "treated")])

  cx <- survival::coxph(survival::Surv(time, event) ~ trt, data = d)
  achieved <- list(
    n_control = nrow(d0), n_treated = nrow(d1),
    event_rate_control = mean(d0$event), event_rate_treated = mean(d1$event),
    median_followup = median(d$time),
    cox_hr = unname(exp(coef(cx)["trt"]))
  )

  list(data = d,
       params = list(alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
                      gamma0 = gamma0, beta1 = beta1, beta2 = beta2,
                      weibull_shape = weibull_shape, tau = tau),
       achieved = achieved, target = stats)
}
