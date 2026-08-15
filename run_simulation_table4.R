## run_simulation_table4.R
##
## Reproduces the Table 4 structure of the manuscript: bias AND empirical 95%
## CI coverage for Cox HR, PHMC incidence OR, and PHMC latency HR, for the
## three named representative scenarios (well-identified / REWIND-calibrated /
## insufficient follow-up). Coverage requires a per-replicate CI, which we
## obtain via `phmc_bootstrap_ci()` (percentile bootstrap) -- this is much
## slower than the point-estimate-only PH-S/NPH-S scans, so this script uses
## fewer replications and fewer bootstrap draws by default. Increase
## n_rep/nboot for a closer match to the manuscript's 500-replication design
## if you have the compute budget (see README for expected runtimes).
##
## NOTE: depends on R/simulate_data.R, R/fit_wrapper.R, and the DGP constants
## in R/run_simulation_PHS.R (source those first).

TABLE4_SCENARIOS <- data.frame(
  tau = c(8, 5.4, 1.5),
  target_Su = c(0.02, 0.80, 0.95),
  label = c("Well-identified (8yr)", "REWIND-calibrated (5.4yr)", "Insufficient follow-up (1.5yr)"),
  stringsAsFactors = FALSE
)

#' Run one Table-4-style replicate: point estimates + bootstrap CIs for Cox HR
#' (Wald CI) and PHMC OR/HR (percentile bootstrap CI)
run_one_table4_replicate <- function(tau, gamma0, dgp, n_per_trial, nboot) {
  d <- simulate_phmc_trial(n = n_per_trial, alpha0 = dgp$alpha0, alpha1 = DGP_ALPHA1, alpha2 = DGP_ALPHA2,
                            gamma0 = gamma0, beta1 = DGP_BETA1, beta2 = DGP_BETA2,
                            weibull_shape = DGP_SHAPE, tau = tau)
  covs <- c("age", "diabDur", "sex", "bmi", "sbp")

  cx <- tryCatch(fit_cox_conflated(d$time, d$event, d$trt, covariates = covs, data = d), error = function(e) NULL)
  bs <- tryCatch(phmc_bootstrap_ci(d$time, d$event, d$trt, covariates = covs, data = d, nboot = nboot),
                  error = function(e) NULL)

  data.frame(
    tau = tau,
    event_rate = mean(d$event),
    cure_frac_true = mean(d$C_true),
    cox_hr = if (is.null(cx)) NA else cx$hr,
    cox_ci_l = if (is.null(cx)) NA else cx$ci_lower,
    cox_ci_u = if (is.null(cx)) NA else cx$ci_upper,
    phmc_or = if (is.null(bs)) NA else bs$point$or_hat,
    phmc_or_ci_l = if (is.null(bs)) NA else bs$or_ci[1],
    phmc_or_ci_u = if (is.null(bs)) NA else bs$or_ci[2],
    phmc_hr = if (is.null(bs)) NA else bs$point$hr_hat,
    phmc_hr_ci_l = if (is.null(bs)) NA else bs$hr_ci[1],
    phmc_hr_ci_u = if (is.null(bs)) NA else bs$hr_ci[2]
  )
}

#' Run the Table-4-style bias/coverage simulation
#'
#' @param n_rep replications per scenario (manuscript: 500; default here: 100
#'   for tractable runtime -- see README for timing guidance)
#' @param nboot bootstrap draws per replicate for the PHMC percentile CI
#'   (manuscript/smcure default: 100; default here: 40)
run_simulation_table4 <- function(n_rep = 100, n_per_trial = 1000, nboot = 40, seed = 2028,
                                   n_cores = max(1, parallel::detectCores() - 1)) {
  set.seed(seed)
  dgp <- calibrate_dgp(scenario_grid = data.frame(tau = TABLE4_SCENARIOS$tau,
                                                    target_Su = TABLE4_SCENARIOS$target_Su))
  if (.Platform$OS.type == "windows") n_cores <- 1

  cat(sprintf("Running Table-4-style simulation: %d scenarios x %d replicates x %d bootstrap draws\n",
              nrow(TABLE4_SCENARIOS), n_rep, nboot))
  cat("This is the slow path (bootstrap CIs per replicate) -- expect several minutes per scenario.\n")

  all_results <- list()
  for (i in seq_len(nrow(TABLE4_SCENARIOS))) {
    row <- TABLE4_SCENARIOS[i, ]
    gamma0 <- dgp$gamma0_by_scenario[[i]]
    cat(sprintf("  scenario %d (%s) ... ", i, row$label))
    t0 <- Sys.time()
    reps <- parallel::mclapply(seq_len(n_rep), function(r)
      run_one_table4_replicate(row$tau, gamma0, dgp, n_per_trial, nboot), mc.cores = n_cores)
    reps_df <- do.call(rbind, reps)
    reps_df$label <- row$label
    all_results[[i]] <- reps_df
    cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))
  }
  raw <- do.call(rbind, all_results)
  rownames(raw) <- NULL

  TRUE_hr <- exp(DGP_BETA1)   # true latency HR
  TRUE_or <- exp(DGP_ALPHA1)  # true incidence OR

  summarize_scn <- function(df) {
    n_ok <- sum(!is.na(df$cox_hr))
    data.frame(
      label = df$label[1], tau = df$tau[1],
      n_rep_ok = n_ok,
      mean_event_rate = mean(df$event_rate, na.rm = TRUE),
      mean_cure_frac_true = mean(df$cure_frac_true, na.rm = TRUE),
      cox_hr_mean = mean(df$cox_hr, na.rm = TRUE),
      cox_hr_bias = mean(df$cox_hr, na.rm = TRUE) - TRUE_hr,
      cox_hr_coverage = mean(df$cox_ci_l <= TRUE_hr & TRUE_hr <= df$cox_ci_u, na.rm = TRUE),
      phmc_or_mean = mean(df$phmc_or, na.rm = TRUE),
      phmc_or_bias = mean(df$phmc_or, na.rm = TRUE) - TRUE_or,
      phmc_or_coverage = mean(df$phmc_or_ci_l <= TRUE_or & TRUE_or <= df$phmc_or_ci_u, na.rm = TRUE),
      phmc_hr_mean = mean(df$phmc_hr, na.rm = TRUE),
      phmc_hr_bias = mean(df$phmc_hr, na.rm = TRUE) - TRUE_hr,
      phmc_hr_coverage = mean(df$phmc_hr_ci_l <= TRUE_hr & TRUE_hr <= df$phmc_hr_ci_u, na.rm = TRUE)
    )
  }
  summary_df <- do.call(rbind, lapply(split(raw, raw$label), summarize_scn))
  summary_df <- summary_df[match(TABLE4_SCENARIOS$label, summary_df$label), ]
  rownames(summary_df) <- NULL

  list(raw = raw, summary = summary_df, truth = list(hr = TRUE_hr, or = TRUE_or))
}
