## generate_table5_table6.R
##
## Produces Table 5 (Cox vs. PHMC mixture-cure estimates, treatment-policy
## estimand) and Table 6 (treatment-policy vs. while-on-treatment estimand
## comparison) in the manuscript's layout, fit to the calibrated synthetic
## datasets from reconstruct_real_trials.R (see that file's header for what
## "calibrated to the real trial" does and does not mean).
##
## The while-on-treatment estimand in the original manuscript censors
## follow-up at intercurrent events (treatment discontinuation, rescue
## therapy). No real intercurrent-event timing data is available here, so
## Table 6 below simulates a plausible, documented synthetic discontinuation
## process (independent exponential "intercurrent event" time) and censors
## at it -- illustrating the treatment-policy vs. while-on-treatment
## MECHANICS described in the manuscript (Section 2.4.5.2 / 3.2), not a
## real intercurrent-event analysis of either trial.
##
## Depends on: R/simulate_data.R, R/fit_wrapper.R, R/real_trial_stats.R,
## R/reconstruct_real_trials.R (source those first).

#' Build Table 5 for one trial: Cox HR (population) + PHMC incidence
#' (cure rates, risk difference) + PHMC latency (HR among susceptibles)
build_table5_row <- function(trial_name, calib) {
  d <- calib$data
  covs <- NULL  # unadjusted, matching the manuscript's Table 5 (population/unadjusted HR)

  cx <- fit_cox_conflated(d$time, d$event, d$trt, covariates = covs, data = d)
  ph <- fit_phmc_smcure(d$time, d$event, d$trt, covariates = covs, data = d)

  arm_ctrl <- calib$target$arm_labels[["control"]]
  arm_trt  <- calib$target$arm_labels[["treated"]]

  data.frame(
    Trial = trial_name,
    Component = c("Standard Cox Proportional Hazards",
                   "Mixture Cure Model: Incidence (Non-Susceptible)",
                   "Mixture Cure Model: Incidence (Non-Susceptible)",
                   "Mixture Cure Model: Incidence (Non-Susceptible)",
                   "Mixture Cure Model: Latency (Among Susceptible)"),
    Estimand = c("Hazard Ratio", "Cure Rate", "Cure Rate", "Risk Difference", "Hazard Ratio"),
    Group = c("Population", arm_ctrl, arm_trt, paste0(arm_trt, " - ", arm_ctrl), paste0(arm_trt, " / ", arm_ctrl)),
    Estimate = c(round(cx$hr, 3), round(ph$pi0_hat, 3), round(ph$pi1_hat, 3),
                 round(ph$rd_hat, 3), round(ph$hr_hat, 3)),
    CI_lower = c(round(cx$ci_lower, 3), NA, NA, NA, NA),
    CI_upper = c(round(cx$ci_upper, 3), NA, NA, NA, NA),
    stringsAsFactors = FALSE
  )
}

#' Simulate a synthetic while-on-treatment (WoT) version of a calibrated
#' dataset: an independent "intercurrent event" (discontinuation/rescue)
#' time is drawn per subject; follow-up is censored at that time if earlier
#' than the observed time.
apply_synthetic_wot_censoring <- function(d, ice_rate = 0.06, seed = 2) {
  set.seed(seed)
  ice_time <- rexp(nrow(d), rate = ice_rate)  # ~ mean 1/rate years to ICE
  d$time_wot  <- pmin(d$time, ice_time)
  d$event_wot <- ifelse(ice_time < d$time, 0L, d$event)
  d
}

#' Build Table 6 for one trial: treatment-policy vs. while-on-treatment
#' cure rate + latency HR comparison
build_table6_rows <- function(trial_name, calib) {
  d <- calib$data
  arm_ctrl <- calib$target$arm_labels[["control"]]
  arm_trt  <- calib$target$arm_labels[["treated"]]

  ph_tp <- fit_phmc_smcure(d$time, d$event, d$trt, data = d)

  d_wot <- apply_synthetic_wot_censoring(d)
  ph_wot <- fit_phmc_smcure(d_wot$time_wot, d_wot$event_wot, d_wot$trt, data = d_wot)

  data.frame(
    Trial = trial_name,
    Estimand = c("Treatment Policy", "Treatment Policy",
                 "While-on-Treatment", "While-on-Treatment"),
    Treatment = c(arm_ctrl, arm_trt, arm_ctrl, arm_trt),
    Cure_Rate = round(c(ph_tp$pi0_hat, ph_tp$pi1_hat, ph_wot$pi0_hat, ph_wot$pi1_hat), 3),
    Latency_HR = round(c(ph_tp$hr_hat, NA, ph_wot$hr_hat, NA), 3),
    stringsAsFactors = FALSE
  )
}

#' Run both trials and print/return Table 5 and Table 6
run_table5_table6 <- function() {
  results <- list()
  for (nm in names(REAL_TRIAL_STATS)) {
    cat("Calibrating synthetic dataset for", nm, "...\n")
    results[[nm]] <- calibrate_trial(REAL_TRIAL_STATS[[nm]], seed = which(names(REAL_TRIAL_STATS) == nm))
  }

  table5 <- do.call(rbind, lapply(names(results), function(nm) build_table5_row(nm, results[[nm]])))
  table6 <- do.call(rbind, lapply(names(results), function(nm) build_table6_rows(nm, results[[nm]])))

  list(table5 = table5, table6 = table6, calibrations = results)
}
