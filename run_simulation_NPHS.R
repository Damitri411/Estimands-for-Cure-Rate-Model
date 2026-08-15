## run_simulation_NPHS.R
##
## Reproduces the non-proportional-hazards-among-susceptibles (NPH-S)
## robustness scenarios referenced in Section 3.1: the treatment log-HR on the
## susceptible hazard changes over time (crossing hazards), assessing whether
## the PHMC incidence estimand remains robust (it should) while the PHMC
## latency HR (fit under a working PH assumption) becomes biased.
##
## NOTE: this file depends on R/simulate_data.R, R/fit_wrapper.R and the DGP
## constants/helpers defined in R/run_simulation_PHS.R (source that file first).

## ---- NPH-S scenario grid ----
## A smaller set of scenarios than the primary PH-S grid (the manuscript
## presents NPH-S as a secondary robustness check, "Additional panels ... are
## provided in the supplementary material"), spanning mild to strong crossing.
NPHS_GRID <- data.frame(
  tau               = c(8, 5.4, 5.4, 3.5, 3.5),
  target_Su         = c(0.02, 0.80, 0.80, 0.90, 0.90),
  nph_beta1_late    = log(c(1.3, 1.3, 1.6, 1.3, 1.6)),  # effect REVERSES sign after crossing
  nph_crossing_time = c(4, 2.7, 2.7, 1.8, 1.8),
  label = c("Well-identified, mild crossing", "REWIND-like, mild crossing",
            "REWIND-like, strong crossing", "Short follow-up, mild crossing",
            "Short follow-up, strong crossing"),
  stringsAsFactors = FALSE
)

run_one_replicate_nphs <- function(tau, gamma0, nph_beta1_late, nph_crossing_time, dgp, n_per_trial) {
  d <- simulate_phmc_trial(n = n_per_trial, alpha0 = dgp$alpha0, alpha1 = DGP_ALPHA1, alpha2 = DGP_ALPHA2,
                            gamma0 = gamma0, beta1 = DGP_BETA1, beta2 = DGP_BETA2,
                            weibull_shape = DGP_SHAPE, tau = tau,
                            nph = TRUE, nph_beta1_late = nph_beta1_late,
                            nph_crossing_time = nph_crossing_time)
  covs <- c("age", "diabDur", "sex", "bmi", "sbp")

  cx <- tryCatch(fit_cox_conflated(d$time, d$event, d$trt, covariates = covs, data = d),
                 error = function(e) NULL)
  ph <- tryCatch(fit_phmc_smcure(d$time, d$event, d$trt, covariates = covs, data = d),
                 error = function(e) NULL)

  data.frame(
    tau = tau,
    event_rate = mean(d$event),
    cure_frac_true = mean(d$C_true),
    cox_log_hr  = if (is.null(cx)) NA_real_ else cx$log_hr,
    phmc_log_or = if (is.null(ph)) NA_real_ else log(ph$or_hat),
    phmc_log_hr = if (is.null(ph)) NA_real_ else log(ph$hr_hat)
  )
}

#' Run the NPH-S robustness simulation
#'
#' @param n_rep replications per scenario (manuscript default: 500)
#' @param n_per_trial per-trial sample size
#' @return list(raw, summary)
run_simulation_NPHS <- function(n_rep = 500, n_per_trial = 1000, seed = 2027,
                                 n_cores = max(1, parallel::detectCores() - 1)) {
  set.seed(seed)
  dgp <- calibrate_dgp(scenario_grid = data.frame(tau = NPHS_GRID$tau, target_Su = NPHS_GRID$target_Su))
  if (.Platform$OS.type == "windows") n_cores <- 1

  cat(sprintf("Running NPH-S simulation: %d scenarios x %d replicates\n", nrow(NPHS_GRID), n_rep))

  all_results <- list()
  for (i in seq_len(nrow(NPHS_GRID))) {
    row <- NPHS_GRID[i, ]
    gamma0 <- dgp$gamma0_by_scenario[[i]]
    cat(sprintf("  scenario %d (%s) ... ", i, row$label))
    t0 <- Sys.time()
    reps <- parallel::mclapply(seq_len(n_rep), function(r)
      run_one_replicate_nphs(row$tau, gamma0, row$nph_beta1_late, row$nph_crossing_time, dgp, n_per_trial),
      mc.cores = n_cores)
    reps_df <- do.call(rbind, reps)
    reps_df$label <- row$label
    all_results[[i]] <- reps_df
    cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))
  }
  raw <- do.call(rbind, all_results)
  rownames(raw) <- NULL

  ## bias against the EARLY-regime true log-HR (beta1) -- since the treatment
  ## effect changes sign after nph_crossing_time, the "true" summary HR is not
  ## a single number; we report bias against the early-window effect (as the
  ## manuscript does: NPH-S bias is defined relative to the assumed working
  ## PH model's target, illustrating that a single scalar HR is inadequate)
  summarize_scn <- function(df) {
    data.frame(
      label = df$label[1], tau = df$tau[1], n_rep_ok = sum(!is.na(df$phmc_log_hr)),
      mean_event_rate = mean(df$event_rate, na.rm = TRUE),
      cox_bias_log_hr  = mean(df$cox_log_hr, na.rm = TRUE) - DGP_BETA1,
      phmc_bias_log_or = mean(df$phmc_log_or, na.rm = TRUE) - DGP_ALPHA1,
      phmc_bias_log_hr = mean(df$phmc_log_hr, na.rm = TRUE) - DGP_BETA1
    )
  }
  summary_df <- do.call(rbind, lapply(split(raw, raw$label), summarize_scn))
  summary_df <- summary_df[match(NPHS_GRID$label, summary_df$label), ]
  rownames(summary_df) <- NULL

  list(raw = raw, summary = summary_df)
}
