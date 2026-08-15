## make_figure3.R
##
## Reproduces the manuscript's Figure 3: estimation bias (log scale) versus
## the Su(tau) identifiability diagnostic, for the conflated Cox HR and the
## PHMC incidence/latency components. Uses base R graphics only (no ggplot2
## dependency, since ggplot2 is not part of R's default/recommended package
## set and may not be installed in every environment).
##
## Usage:
##   source("R/make_figure3.R")
##   plot_figure3(phs_result)   # phs_result <- run_simulation_PHS(...)
##   # or, to save to file:
##   png("figure3_bias_vs_identifiability.png", width=1000, height=700, res=120)
##   plot_figure3(phs_result)
##   dev.off()

#' @param phs_result the list returned by run_simulation_PHS()
plot_figure3 <- function(phs_result) {
  s <- phs_result$summary
  s <- s[order(s$Su_tau), ]

  ylim <- range(c(s$cox_bias_log_hr, s$phmc_bias_log_or, s$phmc_bias_log_hr), na.rm = TRUE)
  ylim <- ylim + c(-0.05, 0.05) * diff(ylim)

  plot(s$Su_tau, s$cox_bias_log_hr, type = "n", #ylim = ylim,
       xlab = expression("S"[u]*"("*tau*"): conditional survival among susceptibles at end of follow-up"),
       ylab = "Bias (log scale)",
       main = "Estimation Bias vs Identifiability Diagnostic Across Simulation Scenarios")
  abline(h = 0, lty = 2, col = "grey40")
  grid(col = "grey90")

  lines(s$Su_tau, s$cox_bias_log_hr, col = "black", lwd = 2)
  points(s$Su_tau, s$cox_bias_log_hr, pch = 18, col = "black", cex = 1.4)

  lines(s$Su_tau, s$phmc_bias_log_or, col = "steelblue", lwd = 2)
  points(s$Su_tau, s$phmc_bias_log_or, pch = 16, col = "steelblue", cex = 1.3)

  lines(s$Su_tau, s$phmc_bias_log_hr, col = "darkorange", lwd = 2)
  points(s$Su_tau, s$phmc_bias_log_hr, pch = 17, col = "darkorange", cex = 1.3)

  legend("bottomleft",
         legend = c("Cox HR (conflated)", "PHMC Incidence (OR)", "PHMC Latency (HR)"),
         col = c("black", "steelblue", "darkorange"), pch = c(18, 16, 17), lwd = 2, bty = "n")

  mtext(expression("Low S"[u]*"("*tau*") "%<->% " well-identified;  High S"[u]*"("*tau*") "%<->% " poorly identified"),
        side = 3, line = 0.3, cex = 0.85, col = "grey30")
  invisible(s)
}

#' Optional ggplot2 version, used automatically if ggplot2 is installed.
#' Falls back to plot_figure3() (base graphics) otherwise.
plot_figure3_gg <- function(phs_result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; using base-R plot_figure3() instead.")
    return(plot_figure3(phs_result))
  }
  s <- phs_result$summary
  long <- data.frame(
    Su_tau = rep(s$Su_tau, 3),
    bias = c(s$cox_bias_log_hr, s$phmc_bias_log_or, s$phmc_bias_log_hr),
    Estimand = rep(c("Cox HR (conflated)", "PHMC Incidence (OR)", "PHMC Latency (HR)"), each = nrow(s))
  )
  ggplot2::ggplot(long, ggplot2::aes(x = Su_tau, y = bias, color = Estimand, shape = Estimand)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(x = expression("S"[u]*"("*tau*")"), y = "Bias (log scale)",
                  title = "Estimation Bias vs Identifiability Diagnostic Across Simulation Scenarios") +
    ggplot2::theme_minimal()
}
