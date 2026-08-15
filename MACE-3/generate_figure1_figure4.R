## generate_figure1_figure4.R
##
## Produces Figure 1 (Kaplan-Meier curves by trial/arm) and Figure 4
## (forest plot of adjusted Cox hazard ratios) in the manuscript's layout,
## using the calibrated synthetic datasets from reconstruct_real_trials.R.
##
## FIGURE 4 CAVEAT: the source publications report SUBGROUP hazard-ratio
## forest plots (treatment effect within age/sex/region strata etc.), not a
## single multivariable Cox model's per-covariate coefficients. The
## manuscript's own Figure 4, which this reproduces the *structure* of, uses
## the latter (a multivariable model's fitted covariate coefficients). What
## is plotted below is therefore the calibrated synthetic dataset's own
## fitted multivariable Cox model -- illustrative of the display format, not
## an extraction of either trial's actual published subgroup analysis.
##
## Depends on: survival, R/simulate_data.R, R/reconstruct_real_trials.R
## (source those first).

#' Figure 1: Kaplan-Meier curves for both trials, side by side
#' @param calibrations named list as returned by run_table5_table6()$calibrations
plot_figure1_km <- function(calibrations) {
  op <- par(mfrow = c(1, length(calibrations)), mar = c(4.5, 4.5, 3, 1))
  on.exit(par(op))
  for (nm in names(calibrations)) {
    calib <- calibrations[[nm]]
    d <- calib$data
    sf <- survival::survfit(survival::Surv(time, event) ~ trt, data = d)
    plot(sf, col = c("grey40", "steelblue"), lwd = 2, conf.int = FALSE,
         xlab = "Years", ylab = "Survival probability (MACE-3 free)",
         main = paste0(nm, "\nKaplan-Meier (synthetic, calibrated to reported event rates)"),
         ylim = c(min(0.75, min(summary(sf)$surv) - 0.02), 1))
    legend("bottomleft", legend = calib$target$arm_labels[c("control", "treated")],
           col = c("grey40", "steelblue"), lwd = 2, bty = "n", cex = 0.85)
  }
}

#' Figure 4: forest plot of adjusted Cox hazard ratios, one panel per trial
#' @param calibrations named list as returned by run_table5_table6()$calibrations
plot_figure4_forest <- function(calibrations) {
  op <- par(mfrow = c(1, length(calibrations)), mar = c(4.5, 9, 3, 1))
  on.exit(par(op))
  for (nm in names(calibrations)) {
    calib <- calibrations[[nm]]
    d <- calib$data
    fit <- survival::coxph(survival::Surv(time, event) ~ trt + age + diabDur + sex + bmi + sbp, data = d)
    s <- summary(fit)
    hr <- s$conf.int[, "exp(coef)"]
    lo <- s$conf.int[, "lower .95"]
    hi <- s$conf.int[, "upper .95"]
    labs <- rownames(s$conf.int)
    labs[labs == "trt"] <- "Treatment"

    ord <- rev(seq_along(hr))
    xlim <- c(min(0.5, lo), max(2, hi))
    plot(hr[ord], seq_along(hr), xlim = xlim, pch = 16, log = "x",
         yaxt = "n", xlab = "Hazard Ratio (95% CI)", ylab = "",
         main = paste0(nm, "\nAdjusted Cox HRs (synthetic, calibrated)"))
    axis(2, at = seq_along(hr), labels = labs[ord], las = 1, cex.axis = 0.85)
    segments(lo[ord], seq_along(hr), hi[ord], seq_along(hr))
    abline(v = 1, lty = 2, col = "grey50")
  }
}
